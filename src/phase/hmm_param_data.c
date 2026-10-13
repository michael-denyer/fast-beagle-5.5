/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) phase/HmmParamData.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/hmm_param_data.h"

#include <stdlib.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "phase/hmm_updater.h"

/* Markers per block of backward values. */
enum { BWD_BLOCK = 128 };

void hmm_param_data_init(hmm_param_data *hpd, const pbwt_phase_ibs *ibs) {
    const phase_data *pd = ibs->pd;
    int max_states = pd->par->phase_states;
    int n = pd->fpd->n_stage1;
    size_t row_bytes = (size_t)max_states * sizeof *hpd->fwd[0];
    size_t n_blocks = ((size_t)n + BWD_BLOCK - 1) / BWD_BLOCK;
    hpd->pd = pd;
    hpd->n_markers = n;
    hpd->max_states = max_states;
    basic_phase_states_init(&hpd->states, ibs, max_states);
    hpd->al_match = util_malloc(2 * sizeof *hpd->al_match);
    for (int i = 0; i < 2; ++i) {
        hpd->al_match[i] = util_malloc((size_t)n * sizeof **hpd->al_match);
        for (int m = 0; m < n; ++m) hpd->al_match[i][m] = util_malloc((size_t)max_states);
        hpd->fwd[i] = util_malloc(row_bytes);
        hpd->bwd[i] = util_malloc(row_bytes);
        hpd->block_end_bwd[i] = util_malloc(n_blocks * row_bytes);
        hpd->block_bwd[i] = util_malloc(BWD_BLOCK * row_bytes);
        hpd->joint[i] = util_malloc(row_bytes);
        hpd->state[i] = util_malloc(row_bytes);
        hpd->discord[i] = util_malloc(row_bytes);
    }
    hpd->hap2_mismatch = util_malloc((size_t)n * sizeof *hpd->hap2_mismatch);
    hpd->hap2_switch = util_malloc((size_t)n * sizeof *hpd->hap2_switch);
    hpd->em_probs[0] = 1.0f - pd->p_mismatch;
    hpd->em_probs[1] = pd->p_mismatch;
    hpd->mismatch_cnt = 0;
    hpd->sum_mismatch_prob = 0.0;
    hpd->sum_gen_dist = 0.0;
    hpd->sum_switch_prob = 0.0;
}

void hmm_param_data_free(hmm_param_data *hpd) {
    for (int i = 0; i < 2; ++i) {
        for (int m = 0; m < hpd->n_markers; ++m) free(hpd->al_match[i][m]);
        free(hpd->al_match[i]);
        free(hpd->fwd[i]);
        free(hpd->bwd[i]);
        free(hpd->block_end_bwd[i]);
        free(hpd->block_bwd[i]);
        free(hpd->joint[i]);
        free(hpd->state[i]);
        free(hpd->discord[i]);
    }
    free(hpd->al_match);
    free(hpd->hap2_mismatch);
    free(hpd->hap2_switch);
    basic_phase_states_free(&hpd->states);
}

void hmm_param_data_add_estimation_data(hmm_param_data *hpd, param_estimates *pe) {
    param_estimates_add_mismatch_data(pe, hpd->mismatch_cnt, hpd->sum_mismatch_prob);
    param_estimates_add_switch_data(pe, hpd->sum_gen_dist, hpd->sum_switch_prob);
    hpd->mismatch_cnt = 0;
    hpd->sum_mismatch_prob = 0.0;
    hpd->sum_gen_dist = 0.0;
    hpd->sum_switch_prob = 0.0;
}

/* One marker's contribution to the sums, as HmmParamData.fwdUpdate adds it. */
static void add_marker(hmm_param_data *hpd, int m, float mismatch_prob, double switch_prob) {
    ++hpd->mismatch_cnt;
    hpd->sum_mismatch_prob += mismatch_prob;
    if (switch_prob > 0.0) {
        hpd->sum_gen_dist += hpd->pd->fpd->stage1_map.gen_dist[m];
        hpd->sum_switch_prob += switch_prob;
    }
}

/* The per-state terms of haplotype i's forward step at a marker: the new
 * forward value, and the three terms that HmmParamData.fwdUpdate adds into
 * sums. A state without a mismatch adds 0 to the mismatch sum, which leaves
 * it unchanged. The loop has no sum, and `omp simd` states that its
 * iterations are independent, so the compiler can vectorise it. */
static void fwd_values(hmm_param_data *hpd, int i, float p_switch, const float *bwd_m, const uint8_t *al_discord, int n_states, float last_sum) {
    float *restrict fwd = hpd->fwd[i], *restrict joint = hpd->joint[i], *restrict state = hpd->state[i], *restrict discord = hpd->discord[i];
    const float *restrict bwd = bwd_m;
    const uint8_t *restrict d = al_discord;
    float match = hpd->em_probs[0], differ = hpd->em_probs[1];
    float shift = p_switch / n_states;
    float scale = (1.0f - p_switch) / last_sum;
    float no_switch_scale = ((1.0f - p_switch) + shift) / last_sum;
    #pragma omp simd
    for (int k = 0; k < n_states; ++k) {
        float em = d[k] ? differ : match;
        joint[k] = bwd[k] * em * no_switch_scale * fwd[k];
        float f = em * (scale * fwd[k] + shift);
        fwd[k] = f;
        float state_prob = f * bwd[k];
        state[k] = state_prob;
        discord[k] = d[k] ? state_prob : 0.0f;
    }
}

/* HmmParamData.fwdUpdate for both haplotypes at marker m. Each sum adds its
 * states in order; the eight sums advance together. The first haplotype's
 * result goes into the totals, and the second's waits in hap2_mismatch and
 * hap2_switch, because Java adds every marker of the first haplotype before
 * any of the second. */
static void fwd_update(hmm_param_data *hpd, int m, float *const bwd_m[2], int n_states, float sum[2], float h_factor) {
    float p_switch = hpd->pd->p_recomb[m];
    fwd_values(hpd, 0, p_switch, bwd_m[0], hpd->al_match[0][m], n_states, sum[0]);
    fwd_values(hpd, 1, p_switch, bwd_m[1], hpd->al_match[1][m], n_states, sum[1]);
    const float *restrict f0 = hpd->fwd[0], *restrict j0 = hpd->joint[0], *restrict s0 = hpd->state[0], *restrict d0 = hpd->discord[0];
    const float *restrict f1 = hpd->fwd[1], *restrict j1 = hpd->joint[1], *restrict s1 = hpd->state[1], *restrict d1 = hpd->discord[1];
    float fwd_sum0 = 0.0f, joint_sum0 = 0.0f, state_sum0 = 0.0f, mismatch_sum0 = 0.0f;
    float fwd_sum1 = 0.0f, joint_sum1 = 0.0f, state_sum1 = 0.0f, mismatch_sum1 = 0.0f;
    for (int k = 0; k < n_states; ++k) {
        joint_sum0 += j0[k];
        fwd_sum0 += f0[k];
        state_sum0 += s0[k];
        mismatch_sum0 += d0[k];
        joint_sum1 += j1[k];
        fwd_sum1 += f1[k];
        state_sum1 += s1[k];
        mismatch_sum1 += d1[k];
    }
    add_marker(hpd, m, mismatch_sum0 / state_sum0, h_factor * (1.0f - joint_sum0 / state_sum0));
    hpd->hap2_mismatch[m] = mismatch_sum1 / state_sum1;
    hpd->hap2_switch[m] = h_factor * (1.0f - joint_sum1 / state_sum1);
    sum[0] = fwd_sum0;
    sum[1] = fwd_sum1;
}

/* Java saves the backward values of every marker and then runs the forward
 * pass, one haplotype after the other. Here the backward pass saves only each
 * block's last row, and the forward pass recomputes a block's rows from that
 * row with the same updates, so every value is the one Java computes and the
 * rows in use stay in cache. */
static void get_param_data(hmm_param_data *hpd, int n_states) {
    int n = hpd->n_markers;
    const float *p_recomb = hpd->pd->p_recomb;
    uint8_t **al_match0 = hpd->al_match[0], **al_match1 = hpd->al_match[1];
    size_t row = (size_t)hpd->max_states;
    size_t row_bytes = (size_t)n_states * sizeof *hpd->bwd[0];
    for (int i = 0; i < 2; ++i) {
        for (int k = 0; k < n_states; ++k) hpd->bwd[i][k] = 1.0f;
    }
    for (int m = n - 1; m >= 0; --m) {
        if (m < n - 1) {
            hmm_lane lanes[2] = {{NULL, hpd->bwd[0], al_match0[m + 1], hpd->em_probs, NULL, p_recomb[m + 1], n_states},
                    {NULL, hpd->bwd[1], al_match1[m + 1], hpd->em_probs, NULL, p_recomb[m + 1], n_states}};
            hmm_bwd_lanes(lanes, 2);
        }
        if (m == n - 1 || (m + 1) % BWD_BLOCK == 0) {
            for (int i = 0; i < 2; ++i) memcpy(hpd->block_end_bwd[i] + (size_t)(m / BWD_BLOCK) * row, hpd->bwd[i], row_bytes);
        }
    }
    float h_factor = n_states / (n_states - 1.0f);
    for (int i = 0; i < 2; ++i) {
        for (int k = 0; k < n_states; ++k) hpd->fwd[i][k] = 1.0f / n_states;
    }
    float sum[2] = {1.0f, 1.0f};
    for (int start = 0; start < n; start += BWD_BLOCK) {
        int end = start + BWD_BLOCK < n ? start + BWD_BLOCK : n;
        for (int i = 0; i < 2; ++i) {
            memcpy(hpd->block_bwd[i] + (size_t)(end - 1 - start) * row, hpd->block_end_bwd[i] + (size_t)(start / BWD_BLOCK) * row, row_bytes);
        }
        for (int m = end - 2; m >= start; --m) {
            float *bwd_m[2] = {hpd->block_bwd[0] + (size_t)(m - start) * row, hpd->block_bwd[1] + (size_t)(m - start) * row};
            memcpy(bwd_m[0], bwd_m[0] + row, row_bytes);
            memcpy(bwd_m[1], bwd_m[1] + row, row_bytes);
            hmm_lane lanes[2] = {{NULL, bwd_m[0], al_match0[m + 1], hpd->em_probs, NULL, p_recomb[m + 1], n_states},
                    {NULL, bwd_m[1], al_match1[m + 1], hpd->em_probs, NULL, p_recomb[m + 1], n_states}};
            hmm_bwd_lanes(lanes, 2);
        }
        for (int m = start; m < end; ++m) {
            float *bwd_m[2] = {hpd->block_bwd[0] + (size_t)(m - start) * row, hpd->block_bwd[1] + (size_t)(m - start) * row};
            fwd_update(hpd, m, bwd_m, n_states, sum, h_factor);
        }
    }
    for (int m = 0; m < n; ++m) add_marker(hpd, m, hpd->hap2_mismatch[m], hpd->hap2_switch[m]);
}

void hmm_param_data_update(hmm_param_data *hpd, int sample) {
    int n_states = basic_phase_states_ibs_states(&hpd->states, sample, hpd->al_match);
    if (n_states > 1) get_param_data(hpd, n_states);
}
