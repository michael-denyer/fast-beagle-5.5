/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) vcf/SlidingWindow.java,
 * vcf/RefTargSlidingWindow.java, vcf/TargSlidingWindow.java and
 * vcf/Window.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/sliding_window.h"

#include <inttypes.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "bref/bref3_it.h"
#include "vcf/filter_util.h"
#include "vcf/interval_it.h"
#include "vcf/markers.h"
#include "vcf/ref_it.h"
#include "vcf/vcf_it.h"

/* Growable arrays of owned record pointers. */
typedef struct {
    gt_rec **v;
    int n, cap;
} targ_list;

typedef struct {
    ref_gt_rec **v;
    bool *in_targ;
    int n, cap;
} ref_list;

static void targ_add(targ_list *l, gt_rec *rec) {
    if (l->n == l->cap) {
        l->cap = l->cap == 0 ? 1024 : 2 * l->cap;
        l->v = util_realloc(l->v, (size_t)l->cap * sizeof *l->v);
    }
    l->v[l->n++] = rec;
}

static void ref_add(ref_list *l, ref_gt_rec *rec, bool in_targ) {
    if (l->n == l->cap) {
        l->cap = l->cap == 0 ? 1024 : 2 * l->cap;
        l->v = util_realloc(l->v, (size_t)l->cap * sizeof *l->v);
        l->in_targ = util_realloc(l->in_targ, (size_t)l->cap * sizeof *l->in_targ);
    }
    l->in_targ[l->n] = in_targ;
    l->v[l->n++] = rec;
}

/* RefTargSlidingWindow.refIt: the reference reader for the file's extension. */
static sample_file_it ref_file_it_open(const char *path, const str_set *exclude_samples,
        const str_set *exclude_markers, int n_threads) {
    size_t n = strlen(path);
    if (n >= 5 && strcmp(path + n - 5, ".bref") == 0) {
        util_exit("\nERROR: bref format (.bref) is not supported\n       Reference files should be in bref3 format (.bref3)");
    }
    bool bref3 = n >= 6 && strcmp(path + n - 6, ".bref3") == 0;
    bool vcf_name = (n >= 4 && strcmp(path + n - 4, ".vcf") == 0) || (n >= 7 && strcmp(path + n - 7, ".vcf.gz") == 0)
            || (n >= 8 && strcmp(path + n - 8, ".vcf.bgz") == 0);
    if (!bref3 && !vcf_name) {
        fprintf(stderr, "\nERROR: unrecognized reference filename extension: \n"
                "       expected \".bref3\", \".bref4\", \".vcf\", \".vcf.gz\", or \".vcf.bgz\"\n\n");
    }
    return bref3 ? bref3_it_open(path, exclude_samples, exclude_markers)
            : ref_it_open(path, exclude_samples, exclude_markers, n_threads);
}

struct sliding_window {
    str_set *exclude_samples, *exclude_markers;
    sample_file_it targ_it;  /* gt_rec records */
    sample_file_it ref_it;   /* ref_gt_rec records; ops is NULL without a reference */
    genetic_map *gen_map;
    float window_cm, overlap_cm;
    int window_markers, overlap_markers;
    bool impute;
    ref_gt_rec *ref_peek;    /* a record read by a hasNext() check */
    gt_rec *next_targ;
    ref_gt_rec *next_ref;
    targ_list targ_overlap;
    ref_list ref_overlap;
    int window_index;
    int64_t n_read;          /* records taken from both iterators */
    int64_t cum_targ_markers, cum_markers;
    bool done;               /* the reader has read the last window */
    /* The reader thread reads each window while the caller holds the one
     * before, as Java's reader thread does. */
    pthread_t reader;
    pthread_mutex_t mutex;
    pthread_cond_t changed;
    bool ready;              /* ahead, or error, awaits sliding_window_next */
    bool stop;
    window *ahead;           /* the next window, or NULL after the last */
    char *error;             /* the reader's util_exit message */
};

static void *read_ahead(void *arg);

/* Java's `it.hasNext() ? it.next() : null`. A true hasNext() is always
 * followed by next(), so reading the record here consumes the same records. */
static gt_rec *targ_next(sliding_window *sw) {
    ++sw->n_read;
    return sample_file_it_next(sw->targ_it);
}

static ref_gt_rec *ref_next(sliding_window *sw) {
    ++sw->n_read;
    if (sw->ref_peek != NULL) {
        ref_gt_rec *rec = sw->ref_peek;
        sw->ref_peek = NULL;
        return rec;
    }
    return sample_file_it_next(sw->ref_it);
}

static bool ref_has_next(sliding_window *sw) {
    if (sw->ref_peek == NULL) sw->ref_peek = sample_file_it_next(sw->ref_it);
    return sw->ref_peek != NULL;
}

sliding_window *sliding_window_open(const par *p) {
    if (p->ped != NULL) util_exit(PROGRAM ": the ped= parameter is not supported");
    sliding_window *sw = util_malloc(sizeof *sw);
    *sw = (sliding_window){0};
    const chrom_interval *interval = p->has_chrom_int ? &p->chrom_int : NULL;
    sw->exclude_samples = filter_id_set(p->excludesamples);
    sw->exclude_markers = filter_id_set(p->excludemarkers);
    sw->targ_it = interval_it_open(vcf_it_open(p->gt, sw->exclude_samples, sw->exclude_markers, p->nthreads), interval);
    if (p->ref != NULL) {
        sw->ref_it = interval_it_open(ref_file_it_open(p->ref, sw->exclude_samples, sw->exclude_markers, p->nthreads),
                interval);
    }
    sw->gen_map = genetic_map_open(p->map, interval);
    sw->window_cm = p->window;
    sw->window_markers = p->window_markers;
    sw->overlap_cm = p->overlap;
    sw->overlap_markers = p->window_markers >> 2;
    sw->impute = p->impute;
    pthread_mutex_init(&sw->mutex, NULL);
    pthread_cond_init(&sw->changed, NULL);
    if (pthread_create(&sw->reader, NULL, read_ahead, sw) != 0) util_exit(PROGRAM ": cannot create thread");
    return sw;
}

const samples *sliding_window_targ_samples(const sliding_window *sw) {
    return sample_file_it_samples(sw->targ_it);
}

int sliding_window_n_ref_samples(const sliding_window *sw) {
    return sw->ref_it.ops == NULL ? 0 : sample_file_it_samples(sw->ref_it)->n;
}

const genetic_map *sliding_window_gen_map(const sliding_window *sw) {
    return sw->gen_map;
}

static void trace_window(const window *w, double end_cm, int32_t end_pos) {
    const marker_indices *ix = &w->indices;
    uint64_t bits;
    memcpy(&bits, &end_cm, sizeof bits);
    trace_line("T2", "window\t%d\t%s\t%s\t%" PRIx64 "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", w->index,
            chrom_ids_id(w->chrom_index), w->last_window ? "true" : "false", bits, end_pos, ix->n_markers,
            ix->n_targ_markers, ix->overlap_end, ix->overlap_start, ix->prev_splice, ix->next_splice,
            ix->targ_overlap_end, ix->targ_overlap_start, ix->targ_prev_splice, ix->targ_next_splice);
    kstring_t s = {0, 0, NULL};
    kputs("positions\t", &s);
    int n = w->n_ref > 0 ? w->n_ref : w->n_targ;
    for (int j = 0; j < n; ++j) {
        if (j > 0) kputc(',', &s);
        kputw(w->n_ref > 0 ? w->ref[j]->marker.pos : w->targ[j]->marker.pos, &s);
    }
    trace_line("T2", "%s", s.s);
    s.l = 0;
    kputs("targMarkerToMarker\t", &s);
    for (int j = 0; j < ix->n_targ_markers; ++j) {
        if (j > 0) kputc(',', &s);
        kputw(ix->targ_marker_to_marker[j], &s);
    }
    trace_line("T2", "%s", s.s);
    free(s.s);
}

static int first_index_with_pos(const marker *const *m, int index) {
    int32_t pos = m[index]->pos;
    while (index > 0 && m[index - 1]->pos == pos) --index;
    return index;
}

/* The overlap start of a window whose markers are m[0..n): the first marker at
 * or after the base position overlap= cM before end_gen_pos, searching only
 * the last overlap_markers markers. */
static int overlap_start(const sliding_window *sw, const marker *const *m, int n, double end_gen_pos) {
    int chrom = m[n - 1]->chrom_index;
    double start_gen_pos = end_gen_pos - sw->overlap_cm;
    int32_t key = genetic_map_base_pos(sw->gen_map, chrom, start_gen_pos);
    int low = n - sw->overlap_markers > 0 ? n - sw->overlap_markers : 0;
    int high = n - 1;
    while (low <= high) {
        int mid = (int)((unsigned)(low + high) >> 1);
        int32_t mid_pos = m[mid]->pos;
        if (mid_pos < key) low = mid + 1;
        else if (mid_pos > key) high = mid - 1;
        else return first_index_with_pos(m, mid);
    }
    return first_index_with_pos(m, high > 0 ? high : 0);
}

/* Window.chromIndex is the chromosome of the first target marker, which can be
 * an overlap marker from the previous chromosome. */
static window *make_window(int index, bool last, gt_rec **targ, int n_targ, ref_gt_rec **ref, int n_ref) {
    window *w = util_malloc(sizeof *w);
    *w = (window){.index = index, .last_window = last, .chrom_index = targ[0]->marker.chrom_index,
            .n_targ = n_targ, .n_ref = n_ref};
    w->targ = util_malloc((size_t)(n_targ > 0 ? n_targ : 1) * sizeof *w->targ);
    for (int j = 0; j < n_targ; ++j) w->targ[j] = gt_rec_retain(targ[j]);
    if (n_ref > 0) {
        w->ref = util_malloc((size_t)n_ref * sizeof *w->ref);
        for (int j = 0; j < n_ref; ++j) w->ref[j] = ref_gt_rec_retain(ref[j]);
    }
    return w;
}

static const marker **targ_markers(gt_rec **recs, int n) {
    const marker **m = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *m);
    for (int j = 0; j < n; ++j) m[j] = &recs[j]->marker;
    return m;
}

static const marker **ref_markers(ref_gt_rec **recs, int n) {
    const marker **m = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *m);
    for (int j = 0; j < n; ++j) m[j] = &recs[j]->marker;
    return m;
}

/* A window that read no new record and carries all of its records into the
 * next window leaves the input and the overlap unchanged, so every later window
 * would be the same one. Java Beagle repeats it without end; exit instead. */
static void check_advances(const sliding_window *sw, int64_t n_read_before, int ov_start, bool last, int chrom,
        int32_t end_pos) {
    if (!last && ov_start == 0 && sw->n_read == n_read_before) {
        util_exit("The window ending at %s:%d does not advance: its overlap with the next window is the whole\n"
                "window. Increase the window parameter or decrease the overlap parameter.", chrom_ids_id(chrom), end_pos);
    }
}

/* TargSlidingWindow.Reader: one window of target records. */
static window *next_targ_window(sliding_window *sw) {
    int chrom = sw->next_targ->marker.chrom_index;
    double end_cm = genetic_map_gen_pos(sw->gen_map, chrom, sw->next_targ->marker.pos);
    end_cm += sw->targ_overlap.n == 0 ? (double)sw->window_cm : (double)(sw->window_cm - sw->overlap_cm);
    int32_t end_pos = genetic_map_base_pos(sw->gen_map, chrom, end_cm);

    int64_t n_read_before = sw->n_read;
    int overlap_end = sw->targ_overlap.n;
    targ_list recs = sw->targ_overlap;
    sw->targ_overlap = (targ_list){0};
    while (sw->next_targ != NULL && sw->next_targ->marker.chrom_index == chrom && sw->next_targ->marker.pos < end_pos
            && recs.n < sw->window_markers) {
        targ_add(&recs, sw->next_targ);
        sw->next_targ = targ_next(sw);
    }
    /* BasicGT rejects an empty window before its markers are inspected. */
    if (recs.n == 0) util_exit("java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0");
    const marker **m = targ_markers(recs.v, recs.n);
    markers_check(m, recs.n);
    bool last = sw->next_targ == NULL;
    bool chrom_end = last || sw->next_targ->marker.chrom_index != chrom;
    int ov_start = chrom_end ? recs.n
            : overlap_start(sw, m, recs.n, genetic_map_gen_pos(sw->gen_map, m[recs.n - 1]->chrom_index, m[recs.n - 1]->pos));
    free(m);
    check_advances(sw, n_read_before, ov_start, last, chrom, end_pos);
    window *w = make_window(++sw->window_index, last, recs.v, recs.n, NULL, 0);
    marker_indices_init_targ(&w->indices, overlap_end, ov_start, recs.n);
    if (trace_on()) trace_window(w, end_cm, end_pos);
    for (int j = ov_start; j < recs.n; ++j) targ_add(&sw->targ_overlap, gt_rec_retain(recs.v[j]));
    for (int j = 0; j < recs.n; ++j) gt_rec_release(recs.v[j]);
    free(recs.v);
    sw->done = last;
    return w;
}

static _Noreturn void empty_window_error(int chrom, int32_t end_pos, const ref_list *ref) {
    const char *id = chrom_ids_id(chrom);
    if (ref->n == 0) {
        util_exit("The window ending at %s:%d\ncontains no reference markers\n"
                "Do the reference and target VCF files contain the same\nchromosomes in the same order?\n", id, end_pos);
    }
    util_exit("The reference and target VCF files contain no markers in common in the window: \n%s:%d-%d\n"
            "Do both VCF files share any markers in this window?\n"
            "Do both VCF files contain the same chromosomes in the same order?\n", id, ref->v[0]->marker.pos, end_pos);
}

/* RefTargSlidingWindow.Reader: one window of matched target and reference records. */
static window *next_ref_window(sliding_window *sw) {
    int chrom = sw->next_targ->marker.chrom_index;
    /* advanceRefItToChrom */
    while (sw->next_ref != NULL && sw->next_ref->marker.chrom_index != chrom && ref_has_next(sw)) {
        ref_gt_rec_release(sw->next_ref);
        sw->next_ref = ref_next(sw);
    }
    double end_cm = genetic_map_gen_pos(sw->gen_map, sw->next_ref->marker.chrom_index, sw->next_ref->marker.pos);
    end_cm += sw->ref_overlap.n == 0 ? (double)sw->window_cm : (double)(sw->window_cm - sw->overlap_cm);
    int32_t end_pos = genetic_map_base_pos(sw->gen_map, chrom, end_cm);

    /* readWindow */
    int64_t n_read_before = sw->n_read;
    int ref_overlap_end = sw->ref_overlap.n;
    targ_list targ = sw->targ_overlap;
    ref_list ref = sw->ref_overlap;
    sw->targ_overlap = (targ_list){0};
    sw->ref_overlap = (ref_list){0};
    while (sw->next_targ != NULL && sw->next_targ->marker.chrom_index == chrom && sw->next_targ->marker.pos < end_pos
            && ref.n < sw->window_markers) {
        const marker *tm = &sw->next_targ->marker;
        while (sw->next_ref != NULL && sw->next_ref->marker.chrom_index == chrom
                && (sw->next_ref->marker.pos < tm->pos
                    || (sw->next_ref->marker.pos == tm->pos && !marker_equals(tm, &sw->next_ref->marker)))) {
            if (sw->impute) ref_add(&ref, sw->next_ref, false);
            else ref_gt_rec_release(sw->next_ref);
            sw->next_ref = ref_next(sw);
        }
        if (sw->next_ref != NULL && marker_equals(&sw->next_ref->marker, tm)) {
            targ_add(&targ, sw->next_targ);
            ref_add(&ref, sw->next_ref, true);
            sw->next_ref = ref_next(sw);
        } else {
            gt_rec_release(sw->next_targ);
        }
        sw->next_targ = targ_next(sw);
    }
    if (sw->impute) {
        while (sw->next_ref != NULL && sw->next_ref->marker.chrom_index == chrom && sw->next_ref->marker.pos < end_pos
                && ref.n < sw->window_markers) {
            ref_add(&ref, sw->next_ref, false);
            sw->next_ref = ref_next(sw);
        }
    }

    /* window */
    if (targ.n == 0 || ref.n == 0) empty_window_error(chrom, end_pos, &ref);
    const marker **rm = ref_markers(ref.v, ref.n);
    markers_check(rm, ref.n);
    const marker **tm = targ_markers(targ.v, targ.n);
    markers_check(tm, targ.n);
    free(tm);
    bool last = sw->next_targ == NULL || sw->next_ref == NULL;
    bool chrom_end = last || ref.v[0]->marker.chrom_index != sw->next_ref->marker.chrom_index;
    int ov_start = chrom_end ? ref.n
            : overlap_start(sw, rm, ref.n, genetic_map_gen_pos(sw->gen_map, rm[ref.n - 1]->chrom_index, end_pos - 1));
    free(rm);
    check_advances(sw, n_read_before, ov_start, last, chrom, end_pos);
    window *w = make_window(++sw->window_index, last, targ.v, targ.n, ref.v, ref.n);
    marker_indices_init(&w->indices, ref.in_targ, ref.n, ref_overlap_end, ov_start);
    if (trace_on()) trace_window(w, end_cm, end_pos);

    /* The overlap carried into the next window. */
    int targ_ov_start = w->indices.targ_overlap_start;
    for (int j = targ_ov_start; j < targ.n; ++j) targ_add(&sw->targ_overlap, gt_rec_retain(targ.v[j]));
    for (int j = ov_start; j < ref.n; ++j) ref_add(&sw->ref_overlap, ref_gt_rec_retain(ref.v[j]), ref.in_targ[j]);
    for (int j = 0; j < targ.n; ++j) gt_rec_release(targ.v[j]);
    for (int j = 0; j < ref.n; ++j) ref_gt_rec_release(ref.v[j]);
    free(targ.v);
    free(ref.v);
    free(ref.in_targ);
    sw->done = last;
    return w;
}

/* Gives the caller w, or NULL after the last window, and waits until the
 * caller takes it. Returns false once sliding_window_close stops the reader. */
static bool hand_over(sliding_window *sw, window *w) {
    pthread_mutex_lock(&sw->mutex);
    sw->ahead = w;
    sw->ready = true;
    pthread_cond_broadcast(&sw->changed);
    while (sw->ready && !sw->stop) pthread_cond_wait(&sw->changed, &sw->mutex);
    bool stop = sw->stop;
    pthread_mutex_unlock(&sw->mutex);
    return !stop;
}

static void read_windows(void *arg) {
    sliding_window *sw = arg;
    sw->next_targ = targ_next(sw);
    if (sw->ref_it.ops == NULL) {
        if (sw->next_targ == NULL) util_exit("Error: no genotype data");
    } else {
        sw->next_ref = ref_next(sw);
        if (sw->next_targ == NULL || sw->next_ref == NULL) util_exit("no genotype data");
    }
    for (;;) {
        window *w = NULL;
        if (!sw->done) w = sw->ref_it.ops == NULL ? next_targ_window(sw) : next_ref_window(sw);
        if (!hand_over(sw, w) || w == NULL) return;
    }
}

/* The reader thread reads window k+1 only once the caller has taken window k,
 * so at most one window is read ahead. An error reading window k+1 waits
 * until the caller asks for window k+1, so an error the caller raises for
 * window k is reported first. (Java's reader runs ahead of its queue of one
 * window, so which of the two errors Java reports depends on thread timing.) */
static void *read_ahead(void *arg) {
    sliding_window *sw = arg;
    char *error = util_try(read_windows, sw);
    if (error != NULL) {
        pthread_mutex_lock(&sw->mutex);
        sw->error = error;
        sw->ready = true;
        pthread_cond_broadcast(&sw->changed);
        pthread_mutex_unlock(&sw->mutex);
    }
    return NULL;
}

window *sliding_window_next(sliding_window *sw) {
    pthread_mutex_lock(&sw->mutex);
    while (!sw->ready) pthread_cond_wait(&sw->changed, &sw->mutex);
    char *error = sw->error;
    window *w = sw->ahead;
    if (error == NULL && w != NULL) {
        sw->ahead = NULL;
        sw->ready = false;
        pthread_cond_broadcast(&sw->changed);
    }
    pthread_mutex_unlock(&sw->mutex);
    if (error != NULL) util_exit("%s", error);
    if (w == NULL) return NULL;
    sw->cum_targ_markers += w->indices.n_targ_markers - w->indices.targ_overlap_end;
    sw->cum_markers += w->indices.n_markers - w->indices.overlap_end;
    return w;
}

int64_t sliding_window_cum_targ_markers(const sliding_window *sw) {
    return sw->cum_targ_markers;
}

int64_t sliding_window_cum_markers(const sliding_window *sw) {
    return sw->cum_markers;
}

void window_free(window *w) {
    for (int j = 0; j < w->n_targ; ++j) gt_rec_release(w->targ[j]);
    for (int j = 0; j < w->n_ref; ++j) ref_gt_rec_release(w->ref[j]);
    free(w->targ);
    free(w->ref);
    marker_indices_free(&w->indices);
    free(w);
}

void sliding_window_close(sliding_window *sw) {
    pthread_mutex_lock(&sw->mutex);
    sw->stop = true;
    pthread_cond_broadcast(&sw->changed);
    pthread_mutex_unlock(&sw->mutex);
    pthread_join(sw->reader, NULL);
    pthread_mutex_destroy(&sw->mutex);
    pthread_cond_destroy(&sw->changed);
    if (sw->ahead != NULL) window_free(sw->ahead);
    free(sw->error);
    gt_rec_release(sw->next_targ);
    ref_gt_rec_release(sw->ref_peek);
    ref_gt_rec_release(sw->next_ref);
    for (int j = 0; j < sw->targ_overlap.n; ++j) gt_rec_release(sw->targ_overlap.v[j]);
    for (int j = 0; j < sw->ref_overlap.n; ++j) ref_gt_rec_release(sw->ref_overlap.v[j]);
    free(sw->targ_overlap.v);
    free(sw->ref_overlap.v);
    free(sw->ref_overlap.in_targ);
    sample_file_it_close(sw->targ_it);
    if (sw->ref_it.ops != NULL) sample_file_it_close(sw->ref_it);
    genetic_map_free(sw->gen_map);
    str_set_free(sw->exclude_samples);
    str_set_free(sw->exclude_markers);
    free(sw);
}
