/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) blbutil/InputIt.java and
 * blbutil/BGZipIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "blbutil/line_reader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <htslib/bgzf.h>

#include "blbutil/utilities.h"
#include "jcompat/jutf8.h"

#define BUF_SIZE (1 << 16)

struct line_reader {
    char *path;
    const char *name;
    BGZF *bgzf;           /* compressed input */
    FILE *file;           /* uncompressed input */
    bool bgzip_lines;     /* BGZipIt line splitting, else BufferedReader */
    bool skip_lf;         /* BufferedReader: a '\r' ended the last line */
    const char *buf;      /* the bytes not yet split into lines: a bgzf block, or file_buf */
    char file_buf[BUF_SIZE];
    size_t pos, len;
    bool eof;
};

static bool ends_with(const char *s, const char *suffix) {
    size_t n = strlen(s), m = strlen(suffix);
    return n >= m && strcmp(s + n - m, suffix) == 0;
}

/* BGZipIt.beginsWithBgzipBlock */
static bool begins_with_bgzip_block(const char *path) {
    unsigned char b[16];
    FILE *f = fopen(path, "rb");
    if (f == NULL) return false;
    size_t n = fread(b, 1, sizeof b, f);
    fclose(f);
    return n == sizeof b && b[0] == 0x1f && b[1] == 0x8b && b[2] == 8 && (b[3] & 4) != 0
            && b[10] == 6 && b[11] == 0 && b[12] == 'B' && b[13] == 'C' && b[14] == 2 && b[15] == 0;
}

line_reader *line_reader_open(const char *path, int n_threads) {
    line_reader *r = util_malloc(sizeof *r);
    memset(r, 0, sizeof *r);
    r->path = util_strndup(path, strlen(path));
    const char *slash = strrchr(r->path, '/');
    r->name = slash == NULL ? r->path : slash + 1;
    if (ends_with(r->name, ".gz") || ends_with(r->name, ".bgz")) {
        /* GZIPInputStream rejects input without the gzip magic bytes, where
         * htslib would read it as plain text. */
        unsigned char magic[2] = {0, 0};
        FILE *f = fopen(path, "rb");
        if (f == NULL) util_exit("Error opening %s", path);
        size_t n = fread(magic, 1, 2, f);
        fclose(f);
        if (n < 2) util_exit("java.io.EOFException: %s", path);
        if (magic[0] != 0x1f || magic[1] != 0x8b) util_exit("java.util.zip.ZipException: Not in GZIP format: %s", path);
        r->bgzip_lines = begins_with_bgzip_block(path);
        r->bgzf = bgzf_open(path, "r");
        if (r->bgzf == NULL) util_exit("Error opening %s", path);
        if (r->bgzip_lines && n_threads > 1 && bgzf_mt(r->bgzf, n_threads, 256) != 0) util_exit("Error opening %s", path);
    } else {
        r->file = fopen(path, "rb");
        if (r->file == NULL) util_exit("Error opening %s", path);
    }
    return r;
}

/* Points buf at the next bytes. A bgzf block is split in place, as
 * bgzf_getline does, instead of being copied out with bgzf_read; the block
 * stays valid until the next bgzf_read_block. */
static bool fill(line_reader *r) {
    if (r->eof) return false;
    ssize_t n;
    if (r->bgzf != NULL) {
        BGZF *fp = r->bgzf;
        if (bgzf_read_block(fp) != 0) util_exit("Error reading %s", r->path);
        n = fp->block_length - fp->block_offset;
        r->buf = (const char *)fp->uncompressed_block + fp->block_offset;
        fp->block_offset = 0;
        fp->block_length = 0;
    } else {
        n = (ssize_t)fread(r->file_buf, 1, BUF_SIZE, r->file);
        if (n == 0 && ferror(r->file)) n = -1;
        r->buf = r->file_buf;
    }
    if (n < 0) util_exit("Error reading %s", r->path);
    r->pos = 0;
    r->len = (size_t)n;
    r->eof = n == 0;
    return n > 0;
}

bool line_reader_next_raw(line_reader *r, kstring_t *line) {
    line->l = 0;
    if (kputsn("", 0, line) < 0) util_oom();  /* line->s is never NULL */
    bool any = false;
    for (;;) {
        if (r->pos == r->len && !fill(r)) break;
        if (r->skip_lf) {
            r->skip_lf = false;
            if (r->buf[r->pos] == '\n') {
                ++r->pos;
                continue;
            }
        }
        any = true;
        const char *start = r->buf + r->pos;
        size_t avail = r->len - r->pos;
        const char *lf = memchr(start, '\n', avail);
        size_t i = lf == NULL ? avail : (size_t)(lf - start);
        if (!r->bgzip_lines) {
            const char *cr = memchr(start, '\r', i);
            if (cr != NULL) i = (size_t)(cr - start);
        }
        if (kputsn(start, i, line) < 0) util_oom();
        if (i == avail) {
            r->pos = r->len;
            continue;
        }
        r->skip_lf = !r->bgzip_lines && start[i] == '\r';
        r->pos += i + 1;
        if (r->bgzip_lines && line->l > 0 && line->s[line->l - 1] == '\r') {
            line->s[--line->l] = '\0';
        }
        return true;
    }
    /* BGZipIt drops bytes after the last '\n'; BufferedReader returns them. */
    return any && !r->bgzip_lines;
}

bool line_reader_next(line_reader *r, kstring_t *line) {
    if (!line_reader_next_raw(r, line)) return false;
    jutf8_sanitize(line);
    return true;
}

const char *line_reader_name(const line_reader *r) {
    return r->name;
}

void line_reader_close(line_reader *r) {
    if (r->bgzf != NULL) bgzf_close(r->bgzf);
    if (r->file != NULL) fclose(r->file);
    free(r->path);
    free(r);
}
