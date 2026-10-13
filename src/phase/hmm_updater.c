/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) phase/HmmUpdater.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/hmm_updater.h"

/* The value loops below have no sum, so the compiler can vectorise them.
 * `omp simd` states that a loop's iterations are independent, so GCC at -O2
 * vectorises it without a runtime alias check. A mismatch byte is 0 or 1, so
 * the select equals p_mismatch[mismatch[s]].
 *
 * A sum is a dependent chain of float adds in state order, which
 * docs/byte-identity.md forbids vectorising. Interleaving the k lanes' chains
 * gives the CPU k chains to overlap, but only when k is a compile-time
 * constant, so that the k accumulators stay in registers: with a runtime k
 * each add goes through memory and a chr20 phasing run is 13% slower. */

/* Each lane's states in order, interleaving the lanes up to the shortest
 * lane and then finishing each longer lane on its own. */
static inline void sum_lanes_k(const hmm_lane *lanes, const int k, float *restrict sums) {
    int n = lanes[0].n_states;
    for (int l = 1; l < k; ++l) {
        if (lanes[l].n_states < n) n = lanes[l].n_states;
    }
    float acc[HMM_MAX_LANES];
    for (int l = 0; l < k; ++l) acc[l] = 0.0f;
    for (int s = 0; s < n; ++s) {
        for (int l = 0; l < k; ++l) acc[l] += lanes[l].val[s];
    }
    for (int l = 0; l < k; ++l) {
        const float *restrict v = lanes[l].val;
        float a = acc[l];
        for (int s = n; s < lanes[l].n_states; ++s) a += v[s];
        sums[l] = a;
    }
}

/* The lane counts the callers use, each inlined with its k constant. */
static void sum_lanes(const hmm_lane *lanes, int k, float *restrict sums) {
    switch (k) {
    case 1: sum_lanes_k(lanes, 1, sums); break;
    case 2: sum_lanes_k(lanes, 2, sums); break;
    case 3: sum_lanes_k(lanes, 3, sums); break;
    case 6: sum_lanes_k(lanes, 6, sums); break;
    default: sum_lanes_k(lanes, k, sums); break;
    }
}

void hmm_fwd_lanes(const hmm_lane *lanes, int k) {
    for (int l = 0; l < k; ++l) {
        const hmm_lane *ln = &lanes[l];
        const uint8_t *restrict d = ln->mismatch;
        float match = ln->p_mismatch[0], differ = ln->p_mismatch[1];
        int n_states = ln->n_states;
        float shift = ln->p_switch / n_states;
        float scale = (1.0f - ln->p_switch) / *ln->sum;
        if (ln->prev == ln->val) {
            float *restrict f = ln->val;
            #pragma omp simd
            for (int s = 0; s < n_states; ++s) f[s] = (d[s] ? differ : match) * (scale * f[s] + shift);
        } else {
            const float *restrict p = ln->prev;
            float *restrict f = ln->val;
            #pragma omp simd
            for (int s = 0; s < n_states; ++s) f[s] = (d[s] ? differ : match) * (scale * p[s] + shift);
        }
    }
    float sums[HMM_MAX_LANES];
    sum_lanes(lanes, k, sums);
    for (int l = 0; l < k; ++l) *lanes[l].sum = sums[l];
}

void hmm_bwd_lanes(const hmm_lane *lanes, int k) {
    for (int l = 0; l < k; ++l) {
        float *restrict b = lanes[l].val;
        const uint8_t *restrict d = lanes[l].mismatch;
        float match = lanes[l].p_mismatch[0], differ = lanes[l].p_mismatch[1];
        int n_states = lanes[l].n_states;
        #pragma omp simd
        for (int s = 0; s < n_states; ++s) b[s] *= d[s] ? differ : match;
    }
    float sums[HMM_MAX_LANES];
    sum_lanes(lanes, k, sums);
    for (int l = 0; l < k; ++l) {
        float *restrict b = lanes[l].val;
        int n_states = lanes[l].n_states;
        float shift = lanes[l].p_switch / n_states;
        float scale = (1.0f - lanes[l].p_switch) / sums[l];
        #pragma omp simd
        for (int s = 0; s < n_states; ++s) b[s] = scale * b[s] + shift;
    }
}
