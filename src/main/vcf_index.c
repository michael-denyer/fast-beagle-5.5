/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "main/vcf_index.h"

#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include "blbutil/str_set.h"
#include "blbutil/utilities.h"

typedef struct {
    int tid;
    hts_pos_t beg, end;
    uint64_t u_end;
} index_rec;

struct vcf_index {
    char *vcf_path, *tbi_path;
    uint64_t u_header_end;
    str_set *chroms;   /* tid order, as tabix numbers them */
    index_rec *recs;
    size_t n_recs, cap;
};

vcf_index *vcf_index_new(const char *vcf_path, const char *tbi_path, uint64_t u_header_end) {
    remove(tbi_path);
    vcf_index *x = util_malloc(sizeof *x);
    *x = (vcf_index){.vcf_path = util_strndup(vcf_path, strlen(vcf_path)),
        .tbi_path = util_strndup(tbi_path, strlen(tbi_path)), .u_header_end = u_header_end,
        .chroms = str_set_new()};
    return x;
}

/* tbx_parse1 (htslib tbx.c) for tbx_conf_vcf, which htslib does not export:
 * 0-based POS, and POS + length(REF) or else an INFO END= past POS. Returns
 * false where tabix skips the line. */
static bool parse(char *line, size_t len, char **chrom_end, hts_pos_t *beg, hts_pos_t *end) {
    *chrom_end = NULL;
    *beg = *end = -1;
    size_t b = 0;
    int id = 1;
    for (size_t i = 0; i <= len; ++i) {
        if (i < len && line[i] != '\t') continue;
        if (id == 1) {
            *chrom_end = line + i;
        } else if (id == 2) {
            char *s;
            *beg = strtoll(line + b, &s, 0);
            if (s == line + b) return false;
            if (--*beg < 0) *beg = 0;
            if (*end < 1) *end = 1;
        } else if (id == 4) {
            if (b < i) *end = *beg + (hts_pos_t)(i - b);
        } else if (id == 8) {
            char c = line[i];
            line[i] = '\0';
            char *s = strstr(line + b, "END=");
            if (s == line + b) s += 4;
            else if (s != NULL && (s = strstr(line + b, ";END=")) != NULL) s += 5;
            if (s != NULL && *s != '.') {
                long long info_end = strtoll(s, NULL, 0);
                if (info_end > *beg) *end = info_end;
            }
            line[i] = c;
        }
        b = i + 1;
        ++id;
    }
    return *chrom_end != NULL && *beg >= 0 && *end >= 0;
}

void vcf_index_add(vcf_index *x, char *line, size_t info_end, uint64_t u_end) {
    /* tabix reads such a line as a header line */
    if (line[0] == '#') {
        if (x->n_recs == 0) x->u_header_end = u_end;
        return;
    }
    char *chrom_end;
    hts_pos_t beg, end;
    if (!parse(line, info_end, &chrom_end, &beg, &end)) return;
    if (x->n_recs == x->cap) {
        x->cap = x->cap == 0 ? 1024 : 2 * x->cap;
        x->recs = util_realloc(x->recs, x->cap * sizeof *x->recs);
    }
    int tid = str_set_index(x->chroms, line, (size_t)(chrom_end - line));
    x->recs[x->n_recs++] = (index_rec){tid, beg, end, u_end};
}

/* The BGZF blocks of a closed file, read from their headers and ISIZE
 * trailers without decompressing. */
typedef struct {
    hFILE *f;
    const char *path;
    uint64_t addr, next;     /* the current block's address and the next block's */
    uint64_t ustart, uend;   /* the current block's uncompressed range */
} bgzf_blocks;

static uint32_t le32(const uint8_t *p) {
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

static bool next_block(bgzf_blocks *bl) {
    uint8_t h[18], isize[4];
    if (hseek(bl->f, (off_t)bl->next, SEEK_SET) < 0) util_exit("Error reading %s", bl->path);
    ssize_t n = hread(bl->f, h, sizeof h);
    if (n == 0) return false;
    if (n != (ssize_t)sizeof h || h[0] != 31 || h[1] != 139 || h[12] != 'B' || h[13] != 'C')
        util_exit("Error reading %s: not a BGZF block at offset %" PRIu64, bl->path, bl->next);
    uint64_t size = (uint64_t)(h[16] | h[17] << 8) + 1;
    if (hseek(bl->f, (off_t)(bl->next + size - 4), SEEK_SET) < 0 || hread(bl->f, isize, 4) != 4)
        util_exit("Error reading %s", bl->path);
    bl->addr = bl->next;
    bl->next += size;
    bl->ustart = bl->uend;
    bl->uend += le32(isize);
    return true;
}

/* bgzf_tell after bgzf_getline has read through uncompressed offset u > 0,
 * for u not decreasing between calls: an offset in the block holding byte
 * u - 1, or the start of the following block once the line has used all of
 * it. */
static uint64_t voffset(bgzf_blocks *bl, uint64_t u) {
    while (u > bl->uend) {
        if (!next_block(bl)) util_exit("Error reading %s: it is shorter than was written", bl->path);
    }
    return u < bl->uend ? bl->addr << 16 | (u - bl->ustart) : bl->next << 16;
}

static void put_le32(uint8_t *p, uint32_t v) {
    for (int j = 0; j < 4; ++j) p[j] = (uint8_t)(v >> 8 * j);
}

/* tbx_set_meta: the tabix configuration, then the NUL-terminated chromosome
 * names in tid order. */
static void set_meta(hts_idx_t *idx, const str_set *chroms) {
    int n = str_set_size(chroms);
    size_t l_nm = 0;
    for (int t = 0; t < n; ++t) l_nm += strlen(str_set_get(chroms, t)) + 1;
    uint8_t *meta = util_malloc(28 + l_nm);
    const tbx_conf_t *c = &tbx_conf_vcf;
    int32_t head[7] = {c->preset, c->sc, c->bc, c->ec, c->meta_char, c->line_skip, (int32_t)l_nm};
    for (int j = 0; j < 7; ++j) put_le32(meta + 4 * j, (uint32_t)head[j]);
    size_t l = 28;
    for (int t = 0; t < n; ++t) {
        size_t k = strlen(str_set_get(chroms, t)) + 1;
        memcpy(meta + l, str_set_get(chroms, t), k);
        l += k;
    }
    if (hts_idx_set_meta(idx, (uint32_t)l, meta, 0) != 0) util_oom();
}

void vcf_index_write(vcf_index *x, uint64_t u_total) {
    const char *vcf_path = x->vcf_path;
    bgzf_blocks bl = {.f = hopen(vcf_path, "r"), .path = vcf_path};
    if (bl.f == NULL) util_exit("Error opening %s", vcf_path);
    hts_idx_t *idx = hts_idx_init(0, HTS_FMT_TBI, voffset(&bl, x->u_header_end), 14, 5);
    if (idx == NULL) util_oom();
    for (size_t j = 0; j < x->n_recs; ++j) {
        const index_rec *r = &x->recs[j];
        if (hts_idx_push(idx, r->tid, r->beg, r->end, voffset(&bl, r->u_end), 1) < 0)
            util_exit("Error indexing %s: records are not sorted", vcf_path);
    }
    if (hts_idx_finish(idx, voffset(&bl, u_total)) != 0) util_exit("Error indexing %s", vcf_path);
    while (next_block(&bl)) {}
    if (bl.uend != u_total) util_exit("Error reading %s: it is longer than was written", vcf_path);
    if (hclose(bl.f) != 0) util_exit("Error reading %s", vcf_path);
    set_meta(idx, x->chroms);
    if (hts_idx_save_as(idx, vcf_path, x->tbi_path, HTS_FMT_TBI) != 0) util_exit("Error writing %s", x->tbi_path);
    hts_idx_destroy(idx);
    str_set_free(x->chroms);
    free(x->recs);
    free(x->vcf_path);
    free(x->tbi_path);
    free(x);
}
