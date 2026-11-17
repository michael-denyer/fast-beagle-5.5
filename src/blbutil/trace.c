/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "blbutil/trace.h"

#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "main/run_outputs.h"

static const char *trace_dir;
static const run_outputs *trace_outputs;
static struct {
    char *name;
    FILE *out;
} *seams;
static int n_seams, seams_cap;
/* parallel_for workers trace too (rev_pbwt_phaser_init), so the seam table
 * and each line are written under one lock. */
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

void trace_init(const char *dir, const run_outputs *outputs) {
    trace_outputs = outputs;
    trace_dir = dir;
}

bool trace_on(void) {
    return trace_dir != NULL;
}

static FILE *seam_file(const char *seam) {
    for (int j = 0; j < n_seams; ++j) {
        if (strcmp(seams[j].name, seam) == 0) return seams[j].out;
    }
    if (n_seams == seams_cap) {
        seams_cap = seams_cap == 0 ? 32 : 2 * seams_cap;
        seams = util_realloc(seams, (size_t)seams_cap * sizeof *seams);
    }
    size_t size = strlen(trace_dir) + strlen(seam) + sizeof "/.txt";
    char *path = util_malloc(size);
    snprintf(path, size, "%s/%s.txt", trace_dir, seam);
    const char *conflict;
    FILE *out = run_outputs_open_trace(trace_outputs, path, &conflict);
    if (out == NULL) {
        pthread_mutex_unlock(&lock);   /* util_exit can return to a util_try */
        if (conflict) util_exit(PROGRAM ": trace file %s aliases %s", path, conflict);
        util_exit("Error opening %s", path);
    }
    free(path);
    seams[n_seams].name = util_strndup(seam, strlen(seam));
    seams[n_seams].out = out;
    ++n_seams;
    return out;
}

void trace_line(const char *seam, const char *fmt, ...) {
    pthread_mutex_lock(&lock);
    FILE *out = seam_file(seam);
    va_list ap;
    va_start(ap, fmt);
    vfprintf(out, fmt, ap);
    va_end(ap);
    fputc('\n', out);
    pthread_mutex_unlock(&lock);
}

void trace_close(void) {
    for (int j = 0; j < n_seams; ++j) {
        if (fclose(seams[j].out) != 0) util_exit("Error writing trace seam %s", seams[j].name);
        free(seams[j].name);
    }
    n_seams = 0;
    free(seams);
    seams = NULL;
    seams_cap = 0;
}
