#!/bin/bash
# Checks the groups CI builds its job matrix from, and that each check's
# declared tier and macos attribute decide where it runs.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
fail=0

expect() {  # label got want
  if [ "$2" = "$3" ]; then
    echo "PASS $1"
  else
    echo "FAIL $1: got '$2', want '$3'"; fail=1
  fi
}
gate() {  # root VAR=value...
  local root=$1; shift
  env "$@" "$root/tests/gate-steps.sh" "$root"
}

expect "full-tier groups" "$(gate "$ROOT" GATE_LIST_GROUPS=1 GATE_TIER=full)" \
  $'core\ncases\njava-jar\njava-source\ntrace\nbgen\nsanitizers\ntsan'
expect "C-tier groups" "$(gate "$ROOT" GATE_LIST_GROUPS=1 GATE_TIER=c)" \
  $'core\nbgen\nsanitizers\ntsan'
expect "undeclared group rejected" \
  "$(gate "$ROOT" GATE_GROUP=not-a-group GATE_LIST=1 > /dev/null 2>&1; echo $?)" 2
expect "full-only group is a C-tier no-op" "$(gate "$ROOT" GATE_TIER=c GATE_GROUP=java-jar GATE_LIST=1)" \
  $'setup fixtures\nsetup c-build'
expect "setup and selected checks keep declaration order" \
  "$(gate "$ROOT" GATE_TIER=full GATE_GROUP=bgen GATE_LIST=1)" $'setup fixtures\nsetup c-build\nbgen bgen'

mkdir -p "$tmp/read-only/tests"
ln -s "$ROOT/tests/gate-steps.sh" "$tmp/read-only/tests/gate-steps.sh"
expect "group listing succeeds and creates no files" \
  "$(gate "$tmp/read-only" GATE_LIST_GROUPS=1 > /dev/null; echo $?) $(ls -A "$tmp/read-only")" "0 tests"

# A fixture with a C-only group and a macOS-only check outside the tsan group.
mkdir -p "$tmp/fixture/tests" "$tmp/linux" "$tmp/darwin"
sed -e 's/^step c core fuzz-regressions /step c c-only fuzz-regressions /' \
    -e 's/^step macos tsan tsan /step macos core tsan /' \
  "$ROOT/tests/gate-steps.sh" > "$tmp/fixture/tests/gate-steps.sh"
chmod +x "$tmp/fixture/tests/gate-steps.sh"
printf '#!/bin/sh\necho Linux\n' > "$tmp/linux/uname"
printf '#!/bin/sh\necho Darwin\n' > "$tmp/darwin/uname"
chmod +x "$tmp/linux/uname" "$tmp/darwin/uname"

expect "C-only group is planned in the C tier" "$(gate "$tmp/fixture" GATE_LIST_GROUPS=1 GATE_TIER=c)" \
  $'core\nbgen\nsanitizers\nc-only'
expect "C-only group is not planned in the full tier" "$(gate "$tmp/fixture" GATE_LIST_GROUPS=1 GATE_TIER=full)" \
  $'core\ncases\njava-jar\njava-source\ntrace\nbgen\nsanitizers'
expect "C-only group runs in the C tier" "$(gate "$tmp/fixture" GATE_TIER=c GATE_GROUP=c-only GATE_LIST=1)" \
  $'setup fixtures\nsetup c-build\nc-only fuzz-regressions'
expect "C-only group is a full-tier no-op" "$(gate "$tmp/fixture" GATE_TIER=full GATE_GROUP=c-only GATE_LIST=1)" \
  $'setup fixtures\nsetup c-build'
expect "macos check skipped off macOS" \
  "$(gate "$tmp/fixture" PATH="$tmp/linux:$PATH" GATE_GROUP=core GATE_LIST=1 | awk '$2 == "tsan"')" ""
expect "macos check runs on macOS" \
  "$(gate "$tmp/fixture" PATH="$tmp/darwin:$PATH" GATE_GROUP=core GATE_LIST=1 | awk '$2 == "tsan"')" "core tsan"
exit "$fail"
