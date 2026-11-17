/* Exercise process-exit cleanup with real files and a concurrent late open. */
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "bgen/bgen_files.h"
#include "blbutil/utilities.h"
#include "main/run_outputs.h"

static run_outputs *outputs;
static pthread_t late_thread;
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t changed = PTHREAD_COND_INITIALIZER;
static bool resume, reopened;

static void fail(void *arg) {
    (void)arg;
    util_exit("caught reader error");
}

static void write_member(run_output_file file, const char *text) {
    FILE *stream = bgen_files_open(outputs, file);
    if (stream == NULL) util_exit("fixture open failed");
    if (fputs(text, stream) == EOF || fclose(stream) != 0) util_exit("fixture write failed");
}

static void *late_open(void *arg) {
    (void)arg;
    pthread_mutex_lock(&mutex);
    while (!resume) pthread_cond_wait(&changed, &mutex);
    pthread_mutex_unlock(&mutex);
    /* A completion notification during abort must not enable a new open. */
    bgen_files_complete();
    FILE *stream = bgen_files_open(outputs, RUN_OUTPUT_SAMPLE);
    reopened = stream != NULL;
    if (stream != NULL) fclose(stream);
    return NULL;
}

/* Registered before the cleanup, so it runs after the cleanup. */
static void after_cleanup(void) {
    pthread_mutex_lock(&mutex);
    resume = true;
    pthread_cond_broadcast(&changed);
    pthread_mutex_unlock(&mutex);
    pthread_join(late_thread, NULL);
    if (reopened) {
        fputs("terminal cleanup allowed a later open\n", stderr);
        _Exit(2);
    }
    run_outputs_free(outputs);
}

int main(int argc, char **argv) {
    if (argc != 3) return 2;
    const char *mode = argv[1];
    bool aborting = strcmp(mode, "abort") == 0;
    if (aborting) atexit(after_cleanup);
    par p = {.out = argv[2], .bgen = BGEN_PHASED};
    outputs = run_outputs_new(&p);
    bgen_files_register_cleanup();
    if (aborting && pthread_create(&late_thread, NULL, late_open, NULL) != 0) return 2;

    FILE *stream = bgen_files_open(outputs, RUN_OUTPUT_BGEN);
    if (stream == NULL) return 2;
    if (strcmp(mode, "caught") == 0) {
        char *error = util_try(fail, NULL);
        if (error == NULL || strcmp(error, "caught reader error") != 0) return 2;
        free(error);
    }
    if (fputs("data after reader error\n", stream) == EOF || fclose(stream) != 0) return 2;
    write_member(RUN_OUTPUT_INFO, "info\n");
    if (strcmp(mode, "partial") == 0 || aborting) util_exit("fixture failure");

    write_member(RUN_OUTPUT_SAMPLE, "sample\n");
    bgen_files_complete();
    if (strcmp(mode, "complete") == 0) util_exit("fixture failure");
    run_outputs_free(outputs);
    return 0;
}
