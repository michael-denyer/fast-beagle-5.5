/* The k-lane HMM step against the scalar form of HmmUpdater.fwdUpdate and
 * bwdUpdate, bit for bit, for every lane count the callers use. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "phase/hmm_updater.h"

#define N_LANES 6
#define MAX_STATES 1000
/* Distinct lengths, with one lane much shorter than the rest, so a step over
 * all six interleaves the sums up to the shortest lane and then finishes
 * each longer lane on its own. */
static const int n_states[N_LANES] = {7, MAX_STATES, 333, 64, 2, 517};

typedef struct {
    float prev[MAX_STATES];
    float val[MAX_STATES];
    uint8_t mismatch[MAX_STATES];
    float p_mismatch[2];
    float p_switch;
    float sum;
} lane_data;

static uint64_t rng = 0x9E3779B97F4A7C15ull;

static uint32_t next_u32(void) {
    rng ^= rng >> 12;
    rng ^= rng << 25;
    rng ^= rng >> 27;
    return (uint32_t)((rng * 0x2545F4914F6CDD1Dull) >> 32);
}

/* Spread over 26 binary orders of magnitude, so a changed summation order
 * changes the sum's bits. */
static float next_value(void) {
    float mantissa = (float)(next_u32() % 1000 + 1) / 1000.0f;
    return ldexpf(mantissa, (int)(next_u32() % 26) - 20);
}

static void fill(lane_data *d, int n) {
    for (int s = 0; s < n; ++s) {
        d->prev[s] = next_value();
        d->val[s] = d->prev[s];
        d->mismatch[s] = (uint8_t)(next_u32() & 1);
    }
    float err = (float)(next_u32() % 100 + 1) / 10000.0f;
    d->p_mismatch[0] = 1.0f - err;
    d->p_mismatch[1] = err;
    d->p_switch = (float)(next_u32() % 1000 + 1) / 10000.0f;
    d->sum = 0.5f + (float)(next_u32() % 1000) / 500.0f;
}

/* HmmUpdater.fwdUpdate as Java writes it: one loop over the states. */
static float ref_fwd(const float *prev, float *fwd, float fwd_sum, float p_switch, const float p_mismatch[2],
        const uint8_t *mismatch, int n) {
    float shift = p_switch / n;
    float scale = (1.0f - p_switch) / fwd_sum;
    float sum = 0.0f;
    for (int s = 0; s < n; ++s) {
        fwd[s] = p_mismatch[mismatch[s]] * (scale * prev[s] + shift);
        sum += fwd[s];
    }
    return sum;
}

/* HmmUpdater.bwdUpdate as Java writes it. */
static void ref_bwd(float *bwd, float p_switch, const float p_mismatch[2], const uint8_t *mismatch, int n) {
    float sum = 0.0f;
    for (int s = 0; s < n; ++s) {
        bwd[s] *= p_mismatch[mismatch[s]];
        sum += bwd[s];
    }
    float scale = (1.0f - p_switch) / sum;
    float shift = p_switch / n;
    for (int s = 0; s < n; ++s) bwd[s] = scale * bwd[s] + shift;
}

static int failures = 0;

static void expect_same(const char *what, int k, int lane, const float *got, const float *want, int n, float got_sum,
        float want_sum) {
    for (int s = 0; s < n; ++s) {
        if (memcmp(&got[s], &want[s], sizeof got[s]) != 0) {
            printf("FAIL %s k=%d lane %d state %d: %.9g, want %.9g\n", what, k, lane, s, (double)got[s], (double)want[s]);
            ++failures;
            return;
        }
    }
    if (memcmp(&got_sum, &want_sum, sizeof got_sum) != 0) {
        printf("FAIL %s k=%d lane %d sum: %.9g, want %.9g\n", what, k, lane, (double)got_sum, (double)want_sum);
        ++failures;
    }
}

int main(void) {
    static lane_data base[N_LANES], ref[N_LANES], out[N_LANES], in_place[N_LANES];
    for (int l = 0; l < N_LANES; ++l) fill(&base[l], n_states[l]);
    for (int l = 0; l < N_LANES; ++l) {
        ref[l] = base[l];
        ref[l].sum = ref_fwd(ref[l].prev, ref[l].val, ref[l].sum, ref[l].p_switch, ref[l].p_mismatch, ref[l].mismatch,
                n_states[l]);
    }
    const int lane_counts[] = {1, 2, 3, 6};
    for (int c = 0; c < 4; ++c) {
        int k = lane_counts[c];
        hmm_lane lanes[N_LANES];
        for (int l = 0; l < k; ++l) {
            out[l] = base[l];
            lanes[l] = (hmm_lane){out[l].prev, out[l].val, out[l].mismatch, out[l].p_mismatch, &out[l].sum, out[l].p_switch,
                    n_states[l]};
        }
        hmm_fwd_lanes(lanes, k);
        for (int l = 0; l < k; ++l) expect_same("forward", k, l, out[l].val, ref[l].val, n_states[l], out[l].sum, ref[l].sum);

        for (int l = 0; l < k; ++l) {
            in_place[l] = base[l];
            lanes[l] = (hmm_lane){in_place[l].val, in_place[l].val, in_place[l].mismatch, in_place[l].p_mismatch,
                    &in_place[l].sum, in_place[l].p_switch, n_states[l]};
        }
        hmm_fwd_lanes(lanes, k);
        for (int l = 0; l < k; ++l) {
            expect_same("in-place forward", k, l, in_place[l].val, ref[l].val, n_states[l], in_place[l].sum, ref[l].sum);
        }
    }
    if (failures == 0) printf("PASS forward lanes match HmmUpdater.fwdUpdate bit for bit at k = 1, 2, 3 and 6, in and out of place\n");

    for (int l = 0; l < N_LANES; ++l) {
        ref[l] = base[l];
        ref_bwd(ref[l].val, ref[l].p_switch, ref[l].p_mismatch, ref[l].mismatch, n_states[l]);
    }
    int bwd_failures_before = failures;
    for (int c = 0; c < 4; ++c) {
        int k = lane_counts[c];
        hmm_lane lanes[N_LANES];
        for (int l = 0; l < k; ++l) {
            out[l] = base[l];
            lanes[l] = (hmm_lane){NULL, out[l].val, out[l].mismatch, out[l].p_mismatch, NULL, out[l].p_switch, n_states[l]};
        }
        hmm_bwd_lanes(lanes, k);
        for (int l = 0; l < k; ++l) expect_same("backward", k, l, out[l].val, ref[l].val, n_states[l], 0.0f, 0.0f);
    }
    if (failures == bwd_failures_before) printf("PASS backward lanes match HmmUpdater.bwdUpdate bit for bit at k = 1, 2, 3 and 6\n");
    return failures == 0 ? 0 : 1;
}
