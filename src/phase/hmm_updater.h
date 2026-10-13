/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) phase/HmmUpdater.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#ifndef PHASE_HMM_UPDATER_H
#define PHASE_HMM_UPDATER_H

#include <stdint.h>

/* One step of the scaled forward Li and Stephens HMM recursion, in float
 * (HmmUpdater.fwdUpdate). p_mismatch[mismatch[k]] is state k's emission
 * probability. Writes the step from prev into fwd and returns its sum. */
float hmm_fwd_update(const float *prev, float *fwd, float fwd_sum, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n_states);

/* The forward step, and the backward step (HmmUpdater.bwdUpdate), for two
 * independent HMMs that share p_switch. Each HMM's values are summed in
 * state order, and the two sums advance together. */
void hmm_fwd_update2(const float *const prev[2], float *const fwd[2], float fwd_sum[2], float p_switch, const float p_mismatch[2],
        const uint8_t *const mismatch[2], const int n_states[2]);
void hmm_bwd_update2(float *const bwd[2], float p_switch, const float p_mismatch[2], const uint8_t *const mismatch[2], const int n_states[2]);

/* The same steps for three HMMs sharing p_switch and n_states. Each HMM's
 * sum adds its states in order. */
void hmm_fwd_update3(float *fwd[3], float fwd_sums[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states);
void hmm_bwd_update3(float *bwd[3], float p_switch, const float p_mismatch[2], const uint8_t *m0, const uint8_t *m1, const uint8_t *m2, int n_states);

/* The arguments of one three-HMM step. sums is the forward sums, read and
 * then set; a backward step does not use it. */
typedef struct {
    float **val;   /* [3][state] */
    float *sums;
    float p_switch;
    const float *p_mismatch;
    const uint8_t *m[3];
    int n_states;
} hmm3_step;

/* hmm_fwd_update3 and hmm_bwd_update3 for two independent samples. Each
 * sum adds its states in order, and the six sums advance together. */
void hmm_fwd_update3x2(const hmm3_step *a, const hmm3_step *b);
void hmm_bwd_update3x2(const hmm3_step *a, const hmm3_step *b);

/* One HMM in a step of the scaled Li and Stephens recursion, in float, over
 * k HMMs at once (HmmUpdater.fwdUpdate and bwdUpdate). p_mismatch[mismatch[s]]
 * is state s's emission probability. The forward step reads prev, writes val
 * (prev may be val) and takes *sum from the previous step's sum to this step's.
 * The backward step reads and writes val and ignores prev and sum. */
typedef struct {
    const float *prev;
    float *val;
    const uint8_t *mismatch;
    const float *p_mismatch;   /* [2]: match, mismatch */
    float *sum;
    float p_switch;
    int n_states;
} hmm_lane;

/* The most lanes in one step. */
#define HMM_MAX_LANES 6

/* One step for lanes[0, k). Each lane's values are computed on their own and
 * its sum adds its states in order, so a lane's result does not depend on
 * which lanes it steps with; the k sums advance together. */
void hmm_fwd_lanes(const hmm_lane *lanes, int k);
void hmm_bwd_lanes(const hmm_lane *lanes, int k);

#endif
