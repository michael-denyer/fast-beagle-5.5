/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) phase/PhaseBaum2.java and
 * phase/SwapRate.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/phase_baum2.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "phase/hmm_updater.h"
#include "phase/marker_cluster.h"

void phase_baum2_init(phase_baum2 *pb, const pbwt_phase_ibs *ibs) {
    phase_data *pd = ibs->pd;
    const par *p = pd->par;
    int n_markers = pd->fpd->n_stage1;
    pb->pd = pd;
    pb->burnin = pd->it < p->burnin;
    pb->lr_threshold = pd->lr_threshold;
    pb->mask_trailing_hets = pb->lr_threshold < 50.0f;
    pb->max_states = p->phase_states;
    basic_phase_states_init(&pb->states, ibs, pb->max_states);
    pb->n_states = 0;
    pb->rows = NULL;
    pb->rows_cap = 0;
    pb->zero_row = util_malloc((size_t)pb->max_states);
    memset(pb->zero_row, 0, (size_t)pb->max_states);
    for (int i = 0; i < 3; ++i) {
        pb->mismatch[i] = util_malloc((size_t)n_markers * sizeof *pb->mismatch[i]);
        pb->fwd[i] = util_malloc((size_t)pb->max_states * sizeof *pb->fwd[i]);
        pb->bwd[i] = util_malloc((size_t)pb->max_states * sizeof *pb->bwd[i]);
    }
    pb->p_mismatch = pd->p_mismatch;
    pb->em_probs[0] = 1.0f - pb->p_mismatch;
    pb->em_probs[1] = pb->p_mismatch;
    pb->n_miss_cap = pb->n_het_cap = 0;
    pb->ref_alleles = NULL;
    pb->bwd_miss1 = pb->bwd_miss2 = pb->bwd_het1 = pb->bwd_het2 = NULL;
    pb->swap_haps = false;
    pb->n_swaps = 0;
}

void phase_baum2_free(phase_baum2 *pb) {
    for (int i = 0; i < 3; ++i) {
        free(pb->mismatch[i]);
        free(pb->fwd[i]);
        free(pb->bwd[i]);
    }
    free(pb->rows);
    free(pb->zero_row);
    for (int j = 0; j < pb->n_miss_cap; ++j) {
        free(pb->ref_alleles[j]);
        free(pb->bwd_miss1[j]);
        free(pb->bwd_miss2[j]);
    }
    for (int j = 0; j < pb->n_het_cap; ++j) {
        free(pb->bwd_het1[j]);
        free(pb->bwd_het2[j]);
    }
    free(pb->ref_alleles);
    free(pb->bwd_miss1);
    free(pb->bwd_miss2);
    free(pb->bwd_het1);
    free(pb->bwd_het2);
    basic_phase_states_free(&pb->states);
}

static void ensure_capacity(phase_baum2 *pb, int n_clusters, int n_unph, int n_miss) {
    size_t row = (size_t)pb->max_states;
    size_t n_rows_bytes = 2 * (size_t)n_clusters * row;
    if (pb->rows_cap < n_rows_bytes) {
        free(pb->rows);
        pb->rows = util_malloc(n_rows_bytes);
        pb->rows_cap = n_rows_bytes;
    }
    if (pb->n_miss_cap < n_miss) {
        pb->ref_alleles = util_realloc(pb->ref_alleles, (size_t)n_miss * sizeof *pb->ref_alleles);
        pb->bwd_miss1 = util_realloc(pb->bwd_miss1, (size_t)n_miss * sizeof *pb->bwd_miss1);
        pb->bwd_miss2 = util_realloc(pb->bwd_miss2, (size_t)n_miss * sizeof *pb->bwd_miss2);
        for (int j = pb->n_miss_cap; j < n_miss; ++j) {
            pb->ref_alleles[j] = util_malloc(row * sizeof **pb->ref_alleles);
            pb->bwd_miss1[j] = util_malloc(row * sizeof **pb->bwd_miss1);
            pb->bwd_miss2[j] = util_malloc(row * sizeof **pb->bwd_miss2);
        }
        pb->n_miss_cap = n_miss;
    }
    if (pb->n_het_cap < n_unph) {
        pb->bwd_het1 = util_realloc(pb->bwd_het1, (size_t)n_unph * sizeof *pb->bwd_het1);
        pb->bwd_het2 = util_realloc(pb->bwd_het2, (size_t)n_unph * sizeof *pb->bwd_het2);
        for (int j = pb->n_het_cap; j < n_unph; ++j) {
            pb->bwd_het1[j] = util_malloc(row * sizeof **pb->bwd_het1);
            pb->bwd_het2[j] = util_malloc(row * sizeof **pb->bwd_het2);
        }
        pb->n_het_cap = n_unph;
    }
}

/* The emission probabilities for a cluster of n markers: a mismatch anywhere
 * in it, capped at 0.5. */
static void set_cluster_em_probs(phase_baum2 *pb, int n_markers_in_cluster) {
    float clust_em = n_markers_in_cluster * pb->p_mismatch;
    if (clust_em >= 0.5f) clust_em = 0.5f;
    pb->em_probs[1] = clust_em;
    pb->em_probs[0] = 1.0f - clust_em;
}

static void copy_states(const phase_baum2 *pb, float *dst, const float *src) {
    memcpy(dst, src, (size_t)pb->n_states * sizeof *dst);
}

/* The backward pass over a sample's clusters, held as a position so that
 * two samples can take their steps together. */
typedef struct {
    phase_baum2 *pb;
    const marker_cluster *mc;
    int c;   /* the cluster the next step reaches, from n_clusters - 2 down to 0 */
    int miss_index;
    int unph_het_index;
} bwd_pass;

static void bwd_begin(bwd_pass *bp, phase_baum2 *pb, const marker_cluster *mc) {
    const sample_phase *sp = mc->sp;
    bp->pb = pb;
    bp->mc = mc;
    bp->miss_index = sp->clust_type_cnt[CLUST_MISSING_GT] + sp->clust_type_cnt[CLUST_MASKED_HET] - 1;
    bp->unph_het_index = sp->clust_type_cnt[CLUST_UNPHASED_HET] - 1;
    for (int k = 0; k < pb->n_states; ++k) pb->bwd[0][k] = 1.0f / pb->n_states;
    copy_states(pb, pb->bwd[1], pb->bwd[0]);
    copy_states(pb, pb->bwd[2], pb->bwd[0]);
    int last_cluster = mc->n_clusters - 1;
    if (marker_cluster_is_missing_or_masked(mc, last_cluster)) {
        copy_states(pb, pb->bwd_miss1[bp->miss_index], pb->bwd[0]);
        copy_states(pb, pb->bwd_miss2[bp->miss_index], pb->bwd[0]);
        --bp->miss_index;
    }
    bp->c = last_cluster - 1;
}

/* The update that takes the backward values from cluster c + 1 to cluster c,
 * as lanes[0, 3). */
static void bwd_step(bwd_pass *bp, hmm_lane *lanes) {
    phase_baum2 *pb = bp->pb;
    const marker_cluster *mc = bp->mc;
    int c_p1 = bp->c + 1;
    set_cluster_em_probs(pb, mc->ends[c_p1] - marker_cluster_start(mc, c_p1));
    for (int i = 0; i < 3; ++i) {
        lanes[i] = (hmm_lane){NULL, pb->bwd[i], pb->mismatch[i][c_p1], pb->em_probs, NULL, mc->p_recomb[c_p1], pb->n_states};
    }
}

static void bwd_end_step(bwd_pass *bp) {
    phase_baum2 *pb = bp->pb;
    int c = bp->c;
    if (marker_cluster_is_missing_or_masked(bp->mc, c)) {
        copy_states(pb, pb->bwd_miss1[bp->miss_index], pb->bwd[1]);
        copy_states(pb, pb->bwd_miss2[bp->miss_index], pb->bwd[2]);
        --bp->miss_index;
    }
    if (bp->mc->sp->clust_type[c + 1] == CLUST_UNPHASED_HET) {
        copy_states(pb, pb->bwd_het1[bp->unph_het_index], pb->bwd[1]);
        copy_states(pb, pb->bwd_het2[bp->unph_het_index], pb->bwd[2]);
        copy_states(pb, pb->bwd[1], pb->bwd[0]);
        copy_states(pb, pb->bwd[2], pb->bwd[0]);
        --bp->unph_het_index;
    }
    --bp->c;
}

static void bwd_finish(bwd_pass *bp) {
    while (bp->c >= 0) {
        hmm_lane lanes[3];
        bwd_step(bp, lanes);
        hmm_bwd_lanes(lanes, 3);
        bwd_end_step(bp);
    }
}

static void swap_haps(phase_baum2 *pb, const marker_cluster *mc, int start_clust, int end_clust) {
    for (int c = start_clust; c < end_clust; ++c) {
        const uint8_t *tmp = pb->mismatch[1][c];
        pb->mismatch[1][c] = pb->mismatch[2][c];
        pb->mismatch[2][c] = tmp;
    }
    sample_phase_swap_haps(mc->sp, marker_cluster_start(mc, start_clust), mc->ends[end_clust - 1]);
}

static void phase_het(phase_baum2 *pb, sample_phase *sp, int unph_het_index, int cluster) {
    const float *b1 = pb->bwd_het1[unph_het_index];
    const float *b2 = pb->bwd_het2[unph_het_index];
    float p11 = 0.0f, p12 = 0.0f, p21 = 0.0f, p22 = 0.0f;
    for (int k = 0; k < pb->n_states; ++k) {
        p11 += pb->fwd[1][k] * b1[k];
        p12 += pb->fwd[1][k] * b2[k];
        p21 += pb->fwd[2][k] * b1[k];
        p22 += pb->fwd[2][k] * b2[k];
    }
    float num = p11 * p22;
    float den = p12 * p21;
    bool last_swap_haps = pb->swap_haps;
    pb->swap_haps = num < den;
    if (pb->swap_haps != last_swap_haps) ++pb->n_swaps;
    if (!pb->burnin && (num >= den * pb->lr_threshold || (pb->swap_haps && den >= num * pb->lr_threshold))) {
        sample_phase_mark_unphased_het_as_phased(sp, cluster);
    }
}

static void impute_missing_gt(sample_phase *sp, int marker, const float *al_freq1, const float *al_freq2, int n_alleles) {
    int a1 = 0;
    int a2 = 0;
    for (int j = 1; j < n_alleles; ++j) {
        if (al_freq1[j] > al_freq1[a1]) a1 = j;
        if (al_freq2[j] > al_freq2[a2]) a2 = j;
    }
    sample_phase_set_allele1(sp, marker, a1);
    sample_phase_set_allele2(sp, marker, a2);
}

static void impute_masked_het(const phase_baum2 *pb, sample_phase *sp, int cluster, int marker, const float *al_freq1, const float *al_freq2) {
    int a1 = sample_phase_allele1(sp, marker);
    int a2 = sample_phase_allele2(sp, marker);
    float p_no_switch = al_freq1[a1] * al_freq2[a2];
    float p_switch = al_freq1[a2] * al_freq2[a1];
    if (p_switch > p_no_switch) {
        sample_phase_set_allele1(sp, marker, a2);
        sample_phase_set_allele2(sp, marker, a1);
        if (p_switch >= pb->lr_threshold * p_no_switch) sample_phase_mark_masked_het_as_phased(sp, cluster);
    } else if (p_no_switch >= pb->lr_threshold * p_switch) {
        sample_phase_mark_masked_het_as_phased(sp, cluster);
    }
}

static void impute_alleles(phase_baum2 *pb, const marker_cluster *mc, int cluster, int miss_index) {
    float *state_probs1 = pb->bwd_miss1[miss_index];
    float *state_probs2 = pb->bwd_miss2[miss_index];
    if (pb->swap_haps) {
        float *tmp = state_probs1;
        state_probs1 = state_probs2;
        state_probs2 = tmp;
    }
    const int *ref_al = pb->ref_alleles[miss_index];
    for (int k = 0; k < pb->n_states; ++k) {
        state_probs1[k] *= pb->fwd[1][k];
        state_probs2[k] *= pb->fwd[2][k];
    }
    int marker = marker_cluster_start(mc, cluster);
    int n_alleles = fpd_n_alleles(pb->pd->fpd, marker);
    float *al_freq1 = util_malloc((size_t)n_alleles * sizeof *al_freq1);
    float *al_freq2 = util_malloc((size_t)n_alleles * sizeof *al_freq2);
    for (int a = 0; a < n_alleles; ++a) al_freq1[a] = al_freq2[a] = 0.0f;
    for (int k = 0; k < pb->n_states; ++k) {
        al_freq1[ref_al[k]] += state_probs1[k];
        al_freq2[ref_al[k]] += state_probs2[k];
    }
    sample_phase *sp = mc->sp;
    if (sp->clust_type[cluster] == CLUST_MISSING_GT) {
        impute_missing_gt(sp, marker, al_freq1, al_freq2, n_alleles);
    } else if (sp->clust_type[cluster] == CLUST_MASKED_HET) {
        impute_masked_het(pb, sp, cluster, marker, al_freq1, al_freq2);
    }
    free(al_freq1);
    free(al_freq2);
}

/* The forward pass, held as a position in the same way. */
typedef struct {
    phase_baum2 *pb;
    const marker_cluster *mc;
    int c;   /* the cluster the next step reaches */
    int miss_index;
    int unph_het_index;
} fwd_pass;

static void fwd_begin(fwd_pass *fp, phase_baum2 *pb, const marker_cluster *mc) {
    *fp = (fwd_pass){pb, mc, 0, 0, 0};
    for (int k = 0; k < pb->n_states; ++k) pb->fwd[0][k] = 1.0f / pb->n_states;
    copy_states(pb, pb->fwd[1], pb->fwd[0]);
    copy_states(pb, pb->fwd[2], pb->fwd[0]);
    pb->fwd_sums[2] = pb->fwd_sums[1] = pb->fwd_sums[0] = 1.0f;
}

/* Phases cluster c if it is an unphased heterozygote, then gives the update
 * that takes the forward values to cluster c, as lanes[0, 3). */
static void fwd_step(fwd_pass *fp, hmm_lane *lanes) {
    phase_baum2 *pb = fp->pb;
    const marker_cluster *mc = fp->mc;
    int c = fp->c;
    if (mc->sp->clust_type[c] == CLUST_UNPHASED_HET) {
        phase_het(pb, mc->sp, fp->unph_het_index, c);
        ++fp->unph_het_index;
        if (pb->swap_haps) {
            int swap_end = fp->unph_het_index < mc->n_unph_het ? mc->unph_het_clusters[fp->unph_het_index] : mc->n_clusters;
            swap_haps(pb, mc, c, swap_end);
        }
        copy_states(pb, pb->fwd[1], pb->fwd[0]);
        copy_states(pb, pb->fwd[2], pb->fwd[0]);
        pb->fwd_sums[1] = pb->fwd_sums[2] = pb->fwd_sums[0];
    }
    set_cluster_em_probs(pb, mc->ends[c] - marker_cluster_start(mc, c));
    for (int i = 0; i < 3; ++i) {
        lanes[i] = (hmm_lane){pb->fwd[i], pb->fwd[i], pb->mismatch[i][c], pb->em_probs, &pb->fwd_sums[i], mc->p_recomb[c], pb->n_states};
    }
}

static void fwd_end_step(fwd_pass *fp) {
    if (marker_cluster_is_missing_or_masked(fp->mc, fp->c)) impute_alleles(fp->pb, fp->mc, fp->c, fp->miss_index++);
    ++fp->c;
}

static void fwd_finish(fwd_pass *fp) {
    while (fp->c < fp->mc->n_clusters) {
        hmm_lane lanes[3];
        fwd_step(fp, lanes);
        hmm_fwd_lanes(lanes, 3);
        fwd_end_step(fp);
    }
}

/* Masks, clusters and builds the states of a sample. False if the sample has
 * nothing left to phase or impute. */
static bool begin_sample(phase_baum2 *pb, int sample, marker_cluster *mc) {
    sample_phase *sp = &pb->pd->phase[sample];
    if (pb->mask_trailing_hets) sample_phase_mask_trailing_unphased_hets(sp);
    int n_unph_hets = sp->clust_type_cnt[CLUST_UNPHASED_HET];
    int n_missing_or_masked = sp->clust_type_cnt[CLUST_MISSING_GT] + sp->clust_type_cnt[CLUST_MASKED_HET];
    if (n_missing_or_masked == 0 && n_unph_hets == 0) return false;
    pb->n_swaps = 0;
    pb->swap_haps = false;
    marker_cluster_init(mc, pb->pd, sample);
    ensure_capacity(pb, mc->n_clusters, n_unph_hets, n_missing_or_masked);
    pb->n_states = basic_phase_states_cluster_states(&pb->states, mc, pb->ref_alleles, pb->mismatch, pb->rows, pb->zero_row);
    return true;
}

static void end_sample(const phase_baum2 *pb, marker_cluster *mc, swap_rate *rate) {
    rate->n_swaps += pb->n_swaps;
    rate->n_unph_hets += mc->n_unph_het;
    marker_cluster_free(mc);
}

/* Each sample is phased as PhaseBaum2.phase phases it. The samples are
 * independent, so their backward passes, and then their forward passes, take
 * their steps together until the shorter one ends. */
void phase_baum2_phase_pair(phase_baum2 *pb0, int sample0, phase_baum2 *pb1, int sample1, swap_rate *rate) {
    marker_cluster mc0 = {0}, mc1 = {0};
    bool run0 = begin_sample(pb0, sample0, &mc0);
    bool run1 = sample1 >= 0 && begin_sample(pb1, sample1, &mc1);
    bwd_pass b0 = {0}, b1 = {0};
    if (run0) bwd_begin(&b0, pb0, &mc0);
    if (run1) bwd_begin(&b1, pb1, &mc1);
    while (run0 && run1 && b0.c >= 0 && b1.c >= 0) {
        hmm_lane lanes[6];
        bwd_step(&b0, lanes);
        bwd_step(&b1, lanes + 3);
        hmm_bwd_lanes(lanes, 6);
        bwd_end_step(&b0);
        bwd_end_step(&b1);
    }
    if (run0) bwd_finish(&b0);
    if (run1) bwd_finish(&b1);
    fwd_pass f0 = {0}, f1 = {0};
    if (run0) fwd_begin(&f0, pb0, &mc0);
    if (run1) fwd_begin(&f1, pb1, &mc1);
    while (run0 && run1 && f0.c < mc0.n_clusters && f1.c < mc1.n_clusters) {
        hmm_lane lanes[6];
        fwd_step(&f0, lanes);
        fwd_step(&f1, lanes + 3);
        hmm_fwd_lanes(lanes, 6);
        fwd_end_step(&f0);
        fwd_end_step(&f1);
    }
    if (run0) fwd_finish(&f0);
    if (run1) fwd_finish(&f1);
    if (run0) end_sample(pb0, &mc0, rate);
    if (run1) end_sample(pb1, &mc1, rate);
}
