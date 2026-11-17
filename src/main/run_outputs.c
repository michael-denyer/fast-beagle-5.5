/* SPDX-License-Identifier: GPL-3.0-or-later */
#define _POSIX_C_SOURCE 200809L
#include "main/run_outputs.h"

#include <errno.h>
#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "blbutil/utilities.h"

/* Also the preflight order, which determines the first collision reported. */
static const char *const suffixes[RUN_OUTPUT_COUNT] = {
    ".vcf.gz", ".log", ".bgen", ".info", ".sample", ".vcf.gz.tbi"
};

struct run_outputs {
    char *paths[RUN_OUTPUT_COUNT];
    const char *inputs[7];
};

static bool enabled(const par *p, run_output_file file) {
    if (file == RUN_OUTPUT_TBI) return p->tbi;
    if (file == RUN_OUTPUT_BGEN || file == RUN_OUTPUT_INFO || file == RUN_OUTPUT_SAMPLE)
        return p->bgen != BGEN_NONE;
    return true;
}

static char *destination(const par *p, run_output_file file) {
    size_t size = strlen(p->out) + strlen(suffixes[file]) + 1;
    char *path = util_malloc(size);
    snprintf(path, size, "%s%s", p->out, suffixes[file]);
    return path;
}

/* java.io.File normalization applies only to comparison copies. Writers use
 * the original spelling, including repeated and trailing slashes. */
static char *file_path(char *s) {
    size_t k = 0;
    for (size_t j = 0; s[j] != '\0'; ++j) {
        if (s[j] != '/' || k == 0 || s[k - 1] != '/') s[k++] = s[j];
    }
    if (k > 1 && s[k - 1] == '/') --k;
    s[k] = '\0';
    return s;
}

static bool same_inode(const char *path, dev_t dev, ino_t ino) {
    struct stat other;
    return path && stat(path, &other) == 0 && dev == other.st_dev && ino == other.st_ino;
}

static void check_file(const run_outputs *out, run_output_file file) {
    const char *path = out->paths[file];
    char *output = file_path(util_strndup(path, strlen(path)));
    struct stat st;
    bool exists = stat(path, &st) == 0;
    for (size_t j = 0; j < sizeof out->inputs / sizeof *out->inputs; ++j) {
        if (out->inputs[j] == NULL) continue;
        char *input = file_path(util_strndup(out->inputs[j], strlen(out->inputs[j])));
        bool same_path = strcmp(output, input) == 0;
        if (file == RUN_OUTPUT_VCF && j < 2 && same_path)
            util_exit("ERROR: VCF output file equals input file: %s", input);
        if (same_path || (exists && same_inode(input, st.st_dev, st.st_ino)))
            util_exit(PROGRAM ": output file %s equals input file %s", output, input);
        free(input);
    }
    free(output);
}

/* A disabled index path may hold an unrelated file (even a declared input).
 * Only a recognizable tabix index needs invalidating on a VCF-only rerun. */
static bool existing_tbi(const char *path) {
    struct stat st;
    if (stat(path, &st) != 0) {
        if (errno == ENOENT || errno == ENOTDIR) return false;
        util_exit("Error checking %s", path);
    }
    if (!S_ISREG(st.st_mode)) return false;
    hFILE *fp = hopen(path, "r");
    if (fp == NULL) util_exit("Error reading %s", path);
    htsFormat fmt;
    int rc = hts_detect_format(fp, &fmt);
    if (hclose(fp) != 0 || rc < 0) util_exit("Error reading %s", path);
    return fmt.format == tbi;
}

run_outputs *run_outputs_new(const par *p) {
    run_outputs *out = util_malloc(sizeof *out);
    *out = (run_outputs){.inputs = {p->ref, p->gt, p->map, p->excludesamples, p->excludemarkers, p->ped, p->truth}};
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file) {
        if (enabled(p, file)) out->paths[file] = destination(p, file);
        else if (file == RUN_OUTPUT_TBI) {
            char *path = destination(p, file);
            if (existing_tbi(path)) out->paths[file] = path;
            else free(path);
        }
    }
    return out;
}

void run_outputs_check(const par *p) {
    struct stat st;
    if (stat(p->out, &st) == 0 && S_ISDIR(st.st_mode))
        util_exit("ERROR: \"out\" parameter cannot be a directory: \"%s\"", p->out);
    run_outputs *out = run_outputs_new(p);
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file)
        if (out->paths[file]) check_file(out, file);
    for (run_output_file a = 0; a < RUN_OUTPUT_COUNT; ++a) {
        if (!out->paths[a] || stat(out->paths[a], &st) != 0 || !S_ISREG(st.st_mode)) continue;
        for (run_output_file b = a + 1; b < RUN_OUTPUT_COUNT; ++b)
            if (same_inode(out->paths[b], st.st_dev, st.st_ino))
                util_exit(PROGRAM ": output file %s equals output file %s", out->paths[a], out->paths[b]);
    }
    run_outputs_free(out);
}

const char *run_outputs_path(const run_outputs *out, run_output_file file) {
    return out->paths[file];
}

/* Opening without O_TRUNC gives a previously dangling link an inode. Check
 * that inode before either writer can change its contents. A trace seam
 * has no fixed output slot and must be checked against every destination. */
static int open_fd(const run_outputs *out, const char *path, run_output_file self, const char **conflict, struct stat *st) {
    *conflict = NULL;
    int fd = open(path, O_WRONLY | O_CREAT, 0666);
    if (fd < 0) return -1;
    if (fstat(fd, st) != 0) goto fail;
    for (size_t j = 0; j < sizeof out->inputs / sizeof *out->inputs; ++j) {
        if (same_inode(out->inputs[j], st->st_dev, st->st_ino)) {
            *conflict = out->inputs[j];
            goto fail;
        }
    }
    if (S_ISREG(st->st_mode)) {
        for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file) {
            if (file != self && same_inode(out->paths[file], st->st_dev, st->st_ino)) {
                *conflict = out->paths[file];
                goto fail;
            }
        }
    }
    return fd;

fail:
    close(fd);
    return -1;
}

static int truncate_output(int fd, const struct stat *st) {
    if (fd >= 0 && S_ISREG(st->st_mode) && ftruncate(fd, 0) != 0) {
        close(fd);
        return -1;
    }
    return fd;
}

int run_outputs_open_vcf(const run_outputs *out) {
    const char *path = out->paths[RUN_OUTPUT_VCF], *conflict;
    struct stat st;
    int fd = open_fd(out, path, RUN_OUTPUT_VCF, &conflict, &st);
    if (conflict) util_exit(PROGRAM ": output file %s aliases %s", path, conflict);
    if (fd < 0) util_exit("Error opening %s", path);
    /* Invalidate the old index just before truncation, including when a
     * later error prevents the new VCF from being completed. */
    const char *index = out->paths[RUN_OUTPUT_TBI];
    if (index && unlink(index) != 0 && errno != ENOENT) {
        close(fd);
        util_exit("Error removing %s", index);
    }
    if (truncate_output(fd, &st) < 0) util_exit("Error opening %s", path);
    return fd;
}

static FILE *open_stream(const run_outputs *out, const char *path, run_output_file self, const char **conflict) {
    struct stat st;
    int fd = open_fd(out, path, self, conflict, &st);
    if (truncate_output(fd, &st) < 0) return NULL;
    FILE *fp = fdopen(fd, "wb");
    if (fp == NULL) close(fd);
    return fp;
}

FILE *run_outputs_open_stream(const run_outputs *out, run_output_file file, const char **conflict) {
    return open_stream(out, out->paths[file], file, conflict);
}

FILE *run_outputs_open_trace(const run_outputs *out, const char *path, const char **conflict) {
    return open_stream(out, path, RUN_OUTPUT_COUNT, conflict);
}

void run_outputs_free(run_outputs *out) {
    for (run_output_file file = 0; file < RUN_OUTPUT_COUNT; ++file) free(out->paths[file]);
    free(out);
}
