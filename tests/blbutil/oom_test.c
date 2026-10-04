/* A failed allocation inside util_try ends the process instead of returning to
 * the try frame, so it cannot unwind past a lock the caller holds. */
/* fork and waitpid are POSIX, which glibc hides under -std=c11. */
#ifndef __APPLE__
#define _POSIX_C_SOURCE 200809L
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "blbutil/utilities.h"

static void allocate_too_much(void *arg) {
    (void)arg;
    util_malloc(SIZE_MAX / 2);
}

static void fail_input(void *arg) {
    (void)arg;
    util_exit("bad input");
}

int main(void) {
    char *error = util_try(fail_input, NULL);
    assert(error != NULL && strcmp(error, "bad input") == 0);
    puts("PASS util_try returns an input error");

    int fds[2];
    assert(pipe(fds) == 0);
    fflush(stdout);
    pid_t child = fork();
    assert(child >= 0);
    if (child == 0) {
        dup2(fds[1], STDERR_FILENO);
        util_try(allocate_too_much, NULL);
        _exit(0);
    }
    close(fds[1]);
    char message[64] = {0};
    size_t n = 0;
    ssize_t got;
    while (n < sizeof message - 1 && (got = read(fds[0], message + n, sizeof message - 1 - n)) > 0) n += (size_t)got;
    int status;
    assert(waitpid(child, &status, 0) == child);
    assert(WIFEXITED(status) && WEXITSTATUS(status) == 1);
    assert(n > 0 && strcmp(message, "ERROR: out of memory\n") == 0);
    puts("PASS a failed allocation inside util_try exits 1 with its message");
    return 0;
}
