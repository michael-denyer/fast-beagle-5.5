#!/usr/bin/env python3
"""No call that can reach util_exit may run while a pthread mutex is held.

util_exit inside util_try longjmps to the try frame, so a mutex locked at the
call stays locked. This check reads src/ as text. It finds each function, the
calls it makes and the mutexes held at each call, and fails on a call under a
mutex that can reach util_exit. util_oom is not such a call: it exits without
unwinding.

A mutex is named by the argument text of its lock call. An unlock in a block
nested inside the lock's block is taken as an early exit, so the mutex is held
again after that block. An unlock of a mutex the function did not lock releases
the caller's mutex until the block ends, which is how a callee may unlock
before its own util_exit.

The check errs towards failing: functions that share a name count as one, and
a call through a struct member or a parameter counts as reaching util_exit. It
does not see a call through a local function-pointer variable, and it takes an
unlock that is the whole body of an if without braces as a release.
"""

import re
import unittest
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UNWINDS = "util_exit"
NEVER_UNWINDS = {"util_oom"}
POINTER = "a function pointer"
LOCK, UNLOCK = "pthread_mutex_lock", "pthread_mutex_unlock"
KEYWORDS = {"if", "while", "for", "switch", "return", "sizeof"}
# Comments, string and character literals, and preprocessor lines.
NOT_CODE = re.compile(
    r"/\*.*?\*/|//[^\n]*|\"(?:\\.|[^\"\\\n])*\"|'(?:\\.|[^'\\\n])*'|^[ \t]*#(?:[^\n]*\\\n)*[^\n]*",
    re.DOTALL | re.MULTILINE,
)
TOKEN = re.compile(r"[A-Za-z_]\w*|->|\S")


@dataclass
class Function:
    name: str
    path: str
    # Callees called while a mutex of this function's caller would still be held.
    reaches: set[str] = field(default_factory=set)
    # (line, callee, mutex) for each call made while this function holds mutex.
    under_lock: list[tuple[int, str, str]] = field(default_factory=list)


def is_name(token):
    return token[0].isalpha() or token[0] == "_"


def closing(tokens, i, open_, close):
    """Return the index of the token that closes the bracket at tokens[i]."""
    depth = 0
    for j in range(i, len(tokens)):
        depth += (tokens[j][0] == open_) - (tokens[j][0] == close)
        if depth == 0:
            return j
    raise ValueError(f"unclosed {open_}")


def scan_body(fn, params, tokens, start, end):
    depth = 1
    held, released = {}, {}  # mutex -> block depth of its lock, and of an unlock that lasts until that block ends
    i = start
    while i < end:
        token, line = tokens[i]
        if token == "{":
            depth += 1
        elif token == "}":
            released = {m: d for m, d in released.items() if d < depth}
            depth -= 1
            held = {m: min(d, depth) for m, d in held.items()}
        elif is_name(token) and token not in KEYWORDS and tokens[i + 1][0] == "(":
            if token in (LOCK, UNLOCK):
                last = closing(tokens, i + 1, "(", ")")
                mutex = "".join(t for t, _ in tokens[i + 2 : last])
                if token == LOCK:
                    released.pop(mutex, None)
                    held.setdefault(mutex, depth)
                elif held.get(mutex) == depth:
                    del held[mutex]
                else:
                    released[mutex] = depth
                i = last
            else:
                callee = POINTER if tokens[i - 1][0] in ("->", ".") or token in params else token
                locked = [m for m in held if m not in released]
                if locked:
                    fn.under_lock.append((line, callee, locked[0]))
                if all(m in held for m in released):
                    fn.reaches.add(callee)
        i += 1


def functions(path, text):
    text = NOT_CODE.sub(lambda m: re.sub(r"[^\n]", " ", m.group()), text)
    # Each token with its line number. No token spans a line.
    tokens = [(m.group(), n) for n, line in enumerate(text.split("\n"), 1) for m in TOKEN.finditer(line)]

    found = []
    head = 0  # the first token of the current top-level declaration
    i = 0
    while i < len(tokens):
        token = tokens[i][0]
        if token == "{":
            end = closing(tokens, i, "{", "}")
            if i > head and tokens[i - 1][0] == ")":
                open_ = max(j for j in range(head, i) if tokens[j][0] == "(" and closing(tokens, j, "(", ")") == i - 1)
                fn = Function(tokens[open_ - 1][0], path)
                params = {t for t, _ in tokens[open_ + 1 : i - 1] if is_name(t)}
                scan_body(fn, params, tokens, i + 1, end)
                found.append(fn)
            else:
                for t, line in tokens[i:end]:
                    if t in (LOCK, UNLOCK, UNWINDS):
                        raise ValueError(f"{path}:{line}: {t} in a block that is not a function body")
            i = end
        if token in ("{", ";"):
            head = i + 1
        i += 1
    return found


def violations(sources):
    """Return one line per call under a mutex that can reach util_exit, for {path: text}."""
    fns = [fn for path, text in sorted(sources.items()) for fn in functions(path, text)]
    via = {UNWINDS: None, POINTER: None}  # function -> the callee through which it reaches util_exit
    grew = True
    while grew:
        grew = False
        for fn in fns:
            exits = sorted(fn.reaches & via.keys())
            if exits and fn.name not in via and fn.name not in NEVER_UNWINDS:
                via[fn.name] = exits[0]
                grew = True

    def chain(name):
        return name if via[name] is None else f"{name} > {chain(via[name])}"

    return [
        f"{fn.path}:{line}: {fn.name} holds {mutex} and calls {chain(callee)}"
        for fn in fns
        for line, callee, mutex in fn.under_lock
        if callee in via
    ]


class LockExit(unittest.TestCase):
    def check(self, code, want):
        self.assertEqual(violations({"a.c": code}), want)

    def test_source_tree(self):
        sources = {str(p.relative_to(ROOT)): p.read_text(encoding="utf-8") for p in (ROOT / "src").rglob("*.[ch]")}
        found = violations(sources)
        if found:
            self.fail("\n" + "\n".join(found))

    def test_exit_under_lock(self):
        self.check(
            """void f(int n) {
    pthread_mutex_lock(&m);
    if (n == 0) util_exit("empty");
    pthread_mutex_unlock(&m);
}""",
            ["a.c:3: f holds &m and calls util_exit"],
        )

    def test_exit_through_callees(self):
        self.check(
            """static void parse(int n) { if (n == 0) util_exit("empty"); }
static void add(int n) { parse(n); }
void f(int n) {
    pthread_mutex_lock(&m);
    add(n);
    pthread_mutex_unlock(&m);
    add(n);
}""",
            ["a.c:5: f holds &m and calls add > parse > util_exit"],
        )

    def test_early_exit_unlocks_for_its_block_only(self):
        self.check(
            """void f(reader *r) {
    for (;;) {
        pthread_mutex_lock(&r->mutex);
        if (r->stop) {
            pthread_mutex_unlock(&r->mutex);
            util_exit("stopped");
        }
        util_exit("still locked");
        pthread_mutex_unlock(&r->mutex);
        util_exit("released");
    }
}""",
            ["a.c:8: f holds &r->mutex and calls util_exit"],
        )

    def test_callee_unlocks_before_its_exit(self):
        code = """static FILE *open_seam(const char *path) {
    FILE *out = fopen(path, "w");
    if (out == NULL) {
        pthread_mutex_unlock(&lock);
        util_exit("Error opening %s", path);
    }
    CHECK_WRITE
    return out;
}
void trace_line(const char *path) {
    pthread_mutex_lock(&lock);
    open_seam(path);
    pthread_mutex_unlock(&lock);
}"""
        self.check(code.replace("CHECK_WRITE", ""), [])
        self.check(
            code.replace("CHECK_WRITE", 'if (ferror(out)) util_exit("Error writing %s", path);'),
            ["a.c:12: trace_line holds &lock and calls open_seam > util_exit"],
        )

    def test_another_mutex_stays_locked(self):
        self.check(
            """void f(void) {
    pthread_mutex_lock(&a);
    pthread_mutex_lock(&b);
    pthread_mutex_unlock(&b);
    util_exit("a is locked");
    pthread_mutex_unlock(&a);
}""",
            ["a.c:5: f holds &a and calls util_exit"],
        )

    def test_allocation_under_lock(self):
        self.check(
            """void util_oom(void) { util_exit("out of memory"); }
void *util_malloc(int n) { void *p = malloc(n); if (p == 0) util_oom(); return p; }
void f(void) {
    pthread_mutex_lock(&m);
    util_malloc(8);
    pthread_mutex_unlock(&m);
}""",
            [],
        )

    def test_call_through_a_pointer(self):
        self.check(
            """static void run(shared *sh, int item) { sh->fn(item); }
void f(shared *sh, callback fn) {
    pthread_mutex_lock(&sh->mu);
    run(sh, 0);
    fn(0);
    pthread_mutex_unlock(&sh->mu);
}""",
            [
                "a.c:4: f holds &sh->mu and calls run > a function pointer",
                "a.c:5: f holds &sh->mu and calls a function pointer",
            ],
        )

    def test_comments_and_strings_are_not_calls(self):
        self.check(
            """void f(void) {
    pthread_mutex_lock(&m);
    /* util_exit("no") would longjmp past the unlock */
    puts("util_exit(1)");  // util_exit(2)
    pthread_mutex_unlock(&m);
}""",
            [],
        )

    def test_lock_outside_a_function_body(self):
        with self.assertRaisesRegex(ValueError, r"a\.c:2: pthread_mutex_lock in a block that is not a function body"):
            violations({"a.c": "static hook hooks[] = {\n    pthread_mutex_lock(&m),\n};"})


if __name__ == "__main__":
    unittest.main()
