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
#ifndef VCF_SLIDING_WINDOW_H
#define VCF_SLIDING_WINDOW_H

#include <stdbool.h>

#include "blbutil/sample_file_it.h"
#include "main/par.h"
#include "vcf/genetic_map.h"
#include "vcf/gt_rec.h"
#include "vcf/marker_indices.h"
#include "vcf/ref_gt_rec.h"
#include "vcf/samples.h"

/* One analysis window. The target records are those also in the reference
 * (or all target records without a reference); the reference records include
 * the reference-only markers to impute. Records at the end of a window are
 * shared with the next window's overlap. */
typedef struct {
    int index;             /* from 1 */
    bool last_window;
    int chrom_index;       /* of the first target marker */
    marker_indices indices;
    int n_targ;
    gt_rec **targ;
    int n_ref;             /* 0 without a reference */
    ref_gt_rec **ref;
} window;

/* Target haplotype hap's phased allele at target marker m of a window, from
 * the phasing result ctx. */
typedef int (*phased_allele_fn)(const void *ctx, int m, int hap);

/* Reads windows of window= cM (at most window-markers= markers) that overlap
 * by overlap= cM, as Beagle's SlidingWindow implementations do. */
typedef struct sliding_window sliding_window;

sliding_window *sliding_window_open(const par *p);
/* fast-beagle only: the windows of open record sources, so a test can cut
 * windows from in-memory records. targ returns gt_rec records; ref returns
 * ref_gt_rec records, or has NULL ops without a reference. Takes ownership of
 * both and of gm. p gives the window parameters. */
sliding_window *sliding_window_open_its(const par *p, sample_file_it targ, sample_file_it ref, genetic_map *gm);
const samples *sliding_window_targ_samples(const sliding_window *sw);
/* The reference samples after excludesamples=, 0 without a reference. */
int sliding_window_n_ref_samples(const sliding_window *sw);
const genetic_map *sliding_window_gen_map(const sliding_window *sw);
/* The next window, or NULL after the last. The caller frees it with window_free.
 * A reader thread reads it while the caller works on the previous window, and
 * an error reading it exits here. */
window *sliding_window_next(sliding_window *sw);
/* SlidingWindow.cumTargMarkers and cumMarkers: the target and all markers of
 * the windows returned so far, each counted once. */
int64_t sliding_window_cum_targ_markers(const sliding_window *sw);
int64_t sliding_window_cum_markers(const sliding_window *sw);
void window_free(window *w);
void sliding_window_close(sliding_window *sw);

#endif
