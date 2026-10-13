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
