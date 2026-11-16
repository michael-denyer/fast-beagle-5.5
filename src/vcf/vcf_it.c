/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) vcf/VcfIt.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/vcf_it.h"

#include <stdlib.h>

#include "blbutil/line_reader.h"
#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jutf8.h"
#include "vcf/filter_util.h"
#include "vcf/gt_rec.h"
#include "vcf/vcf_header.h"

/* VcfIt.DEFAULT_BUFFER_SIZE: lines parsed together on the worker threads. */
#define VCF_IT_BATCH 1024

typedef struct {
    line_reader *reader;
    vcf_header header;
    const str_set *exclude;
    int n_threads;
    kstring_t line;
    bool have_line;     /* line holds the next unread data line */
    /* The batch: recs[next..n) are not yet returned. A parse error, and an
     * error reading the lines after lines[n - 1], are raised when the record
     * they follow is returned, as when each line was parsed on demand. The
     * read error is raised in place of the batch's parse errors, because
     * VcfIt.fillStringBuffer reads a batch's lines before any is parsed. */
    kstring_t lines[VCF_IT_BATCH];
    gt_rec *recs[VCF_IT_BATCH];
    char *errors[VCF_IT_BATCH];
    char *read_error;
    int n, next;
} vcf_it;

void vcf_it_trace_marker(const char *seam, const marker *m) {
    span id = marker_id(m), al = marker_alleles(m), end = marker_end_value(m);
    trace_line(seam, "%s\t%d\t%.*s\t%.*s\t%d\t%.*s", marker_chrom(m), m->pos, id.n, id.s, al.n, al.s,
            marker_n_alleles(m), end.n, end.s);
}

/* VcfIt.readLine skips blank lines, except a final one. */
static bool is_blank(const kstring_t *s) {
    for (size_t j = 0; j < s->l; ++j) {
        if ((unsigned char)s->s[j] > ' ') return false;
    }
    return true;
}

/* Lines are sanitised on the parse workers; a blank line is blank either way. */
static void read_line(vcf_it *it) {
    it->have_line = line_reader_next_raw(it->reader, &it->line);
    kstring_t probe = {0, 0, NULL};
    while (it->have_line && is_blank(&it->line)) {
        if (!line_reader_next_raw(it->reader, &probe)) break;
        kstring_t tmp = it->line;
        it->line = probe;
        probe = tmp;
    }
    free(probe.s);
}

static void trace_rec(const gt_rec *rec) {
    vcf_it_trace_marker("T1a-target", &rec->marker);
    kstring_t s = {0, 0, NULL};
    ksprintf(&s, "%s\t%s\t", gt_rec_class_name(rec), rec->is_phased ? "true" : "false");
    for (int h = 0; h < rec->n_haps; ++h) {
        int a = gt_rec_get(rec, h);
        if (h > 0) kputc(',', &s);
        if (a < 0) kputc('.', &s);
        else kputw(a, &s);
    }
    trace_line("T1c-target", "%s", s.s);
    free(s.s);
}

static void read_lines(void *self) {
    vcf_it *it = self;
    while (it->n < VCF_IT_BATCH && it->have_line) {
        kstring_t t = it->lines[it->n];
        it->lines[it->n++] = it->line;
        it->line = t;
        read_line(it);
    }
}

typedef struct {
    vcf_it *it;
    int j;
} parse_ctx;

static void parse_line(void *arg) {
    const parse_ctx *c = arg;
    kstring_t *line = &c->it->lines[c->j];
    jutf8_sanitize(line);
    gt_rec *rec = util_malloc(sizeof *rec);
    gt_rec_parse(rec, line->s, line->l, &c->it->header);
    rec->refs = 1;
    c->it->recs[c->j] = rec;
}

static void parse_task(void *worker, int j) {
    vcf_it *it = worker;
    it->recs[j] = NULL;
    it->errors[j] = util_try(parse_line, &(parse_ctx){it, j});
}

/* VcfIt.fillEmissionBuffer: reads a batch of lines and parses them in parallel. */
static void fill_batch(vcf_it *it) {
    it->n = it->next = 0;
    it->read_error = util_try(read_lines, it);
    parallel_for(parallel_threads(it->n_threads, it->n), it->n, it, 0, parse_task);
}

/* VcfIt.next: the next record that passes the marker filter. */
static void *vcf_it_next(void *self) {
    vcf_it *it = self;
    for (;;) {
        if (it->next == it->n) {
            if (!it->have_line) return NULL;
            fill_batch(it);
        }
        int j = it->next++;
        if (it->read_error != NULL && (it->errors[j] != NULL || it->next == it->n)) util_exit("%s", it->read_error);
        if (it->errors[j] != NULL) util_exit("%s", it->errors[j]);
        gt_rec *rec = it->recs[j];
        if (filter_accept_marker(it->exclude, &rec->marker)) {
            if (trace_on()) trace_rec(rec);
            return rec;
        }
        gt_rec_release(rec);
    }
}

static const samples *vcf_it_samples(const void *self) {
    const vcf_it *it = self;
    return &it->header.samples;
}

static void vcf_it_close(void *self) {
    vcf_it *it = self;
    for (int j = it->next; j < it->n; ++j) gt_rec_release(it->recs[j]);
    for (int j = 0; j < VCF_IT_BATCH; ++j) {
        free(it->lines[j].s);
        free(it->errors[j]);
    }
    free(it->read_error);
    vcf_header_free(&it->header);
    free(it->line.s);
    line_reader_close(it->reader);
    free(it);
}

static const marker *rec_marker(const void *rec) {
    return &((const gt_rec *)rec)->marker;
}

static void rec_release(void *rec) {
    gt_rec_release(rec);
}

static const sample_file_it_ops vcf_it_ops = {.samples = vcf_it_samples, .next = vcf_it_next, .close = vcf_it_close,
        .marker = rec_marker, .release = rec_release};

sample_file_it vcf_it_open(const char *path, const str_set *exclude_samples, const str_set *exclude_markers,
        int n_threads) {
    vcf_it *it = util_malloc(sizeof *it);
    *it = (vcf_it){0};
    it->reader = line_reader_open(path, n_threads);
    it->exclude = exclude_markers;
    it->n_threads = n_threads;
    vcf_header_read(&it->header, it->reader, &it->line, exclude_samples);
    it->have_line = true;
    if (trace_on()) vcf_header_trace(&it->header, "T1b-target");
    return (sample_file_it){&vcf_it_ops, it};
}
