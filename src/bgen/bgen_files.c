/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "bgen/bgen_files.h"

#include <assert.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdlib.h>

#include "blbutil/utilities.h"

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static const char *created[3];
static int n_created;
static bool finished;

static void remove_partial(void) {
    pthread_mutex_lock(&lock);
    if (!finished) {
        for (int j = 0; j < n_created; ++j) remove(created[j]);
        finished = true;
    }
    pthread_mutex_unlock(&lock);
}

void bgen_files_register_cleanup(void) {
    atexit(remove_partial);
}

FILE *bgen_files_open(const run_outputs *out, run_output_file file) {
    const char *path = run_outputs_path(out, file), *conflict = NULL;
    pthread_mutex_lock(&lock);
    assert(n_created < (int)(sizeof created / sizeof *created));
    FILE *stream = finished ? NULL : run_outputs_open_stream(out, file, &conflict);
    if (stream != NULL) created[n_created++] = path;
    pthread_mutex_unlock(&lock);
    if (conflict) util_exit(PROGRAM ": output file %s aliases %s", path, conflict);
    return stream;
}

void bgen_files_complete(void) {
    pthread_mutex_lock(&lock);
    finished = true;
    pthread_mutex_unlock(&lock);
}
