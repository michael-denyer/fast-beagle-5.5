/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) phase/HmmStateProbs.java,
 * phase/Stage2Baum.java, phase/Stage2Haps.java and PhaseLS.runStage2;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "phase/stage2.h"

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>

#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"
#include "phase/hmm_updater.h"
#include "phase/low_freq_phase_ibs.h"
#include "phase/low_freq_phase_states.h"

/* HmmStateProbs: forward-backward state probabilities of a sample's two
 * haplotypes at stage-1 markers. Java runs one haplotype at a time and keeps
 * a row for every marker. Only the markers with slot[m] >= 0 are read
 * afterwards, so only they get row slot[m] of states and probs; the forward
 * and backward values at every other marker are computed and dropped. The
 * two haplotypes' HMMs are independent and step through the markers
 * together. The last marker keeps its forward values, as in Java. */
typedef struct {
    const phase_data *pd;
    low_freq_phase_states states;
    int n_markers;
    int max_states;
    uint8_t *mismatch[2];   /* per haplotype: [marker][state] */
    float *fwd[2][2];       /* per haplotype: the forward values at a marker without a row, by marker parity */
    float *bwd[2];
    float p_mismatch[2];
} hmm_state_probs;

static void hmm_state_probs_init(hmm_state_probs *hp, const low_freq_phase_ibs *ibs) {
    const phase_data *pd = ibs->pd;
    int max_states = pd->par->phase_states / 2;
    hp->pd = pd;
    hp->n_markers = pd->fpd->n_stage1;
    hp->max_states = max_states;
    low_freq_phase_states_init(&hp->states, ibs, max_states);
    for (int i = 0; i < 2; ++i) {
        hp->mismatch[i] = util_malloc((size_t)hp->n_markers * (size_t)max_states);
        hp->fwd[i][0] = util_malloc((size_t)max_states * sizeof *hp->fwd[i][0]);
        hp->fwd[i][1] = util_malloc((size_t)max_states * sizeof *hp->fwd[i][1]);
        hp->bwd[i] = util_malloc((size_t)max_states * sizeof *hp->bwd[i]);
    }
    hp->p_mismatch[0] = 1.0f - pd->p_mismatch;
    hp->p_mismatch[1] = pd->p_mismatch;
}

static void hmm_state_probs_free(hmm_state_probs *hp) {
    for (int i = 0; i < 2; ++i) {
        free(hp->mismatch[i]);
        free(hp->fwd[i][0]);
        free(hp->fwd[i][1]);
        free(hp->bwd[i]);
    }
    low_freq_phase_states_free(&hp->states);
}

static const uint8_t *mismatch_row(const hmm_state_probs *hp, int i, int m) {
    return hp->mismatch[i] + (size_t)m * (size_t)hp->max_states;
}

static float *fwd_row(const hmm_state_probs *hp, int i, const int *slot, float *probs, int m) {
    return slot[m] >= 0 ? probs + (size_t)slot[m] * (size_t)hp->max_states : hp->fwd[i][m & 1];
}

static void run_fwd(const hmm_state_probs *hp, const int *slot, float *const probs[2], const int n_states[2]) {
    const float *p_recomb = hp->pd->p_recomb;
    const float *prev[2] = {NULL, NULL};
    float sum[2] = {0.0f, 0.0f};
    for (int i = 0; i < 2; ++i) {
        float *row = fwd_row(hp, i, slot, probs[i], 0);
        for (int j = 0; j < n_states[i]; ++j) {
            row[j] = hp->p_mismatch[hp->mismatch[i][j]];
            sum[i] += row[j];
        }
        prev[i] = row;
    }
    for (int m = 1; m < hp->n_markers; ++m) {
        float *fwd[2] = {fwd_row(hp, 0, slot, probs[0], m), fwd_row(hp, 1, slot, probs[1], m)};
        const uint8_t *mismatch[2] = {mismatch_row(hp, 0, m), mismatch_row(hp, 1, m)};
        hmm_lane lanes[2];
        for (int i = 0; i < 2; ++i) lanes[i] = (hmm_lane){prev[i], fwd[i], mismatch[i], hp->p_mismatch, &sum[i], p_recomb[m], n_states[i]};
        hmm_fwd_lanes(lanes, 2);
        prev[0] = fwd[0];
        prev[1] = fwd[1];
    }
}

static void run_bwd(hmm_state_probs *hp, const int *slot, float *const probs[2], const int n_states[2]) {
    const float *p_recomb = hp->pd->p_recomb;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < n_states[i]; ++j) hp->bwd[i][j] = 1.0f / n_states[i];
    }
    for (int m = hp->n_markers - 2; m >= 0; --m) {
        hmm_lane lanes[2];
        for (int i = 0; i < 2; ++i) {
            lanes[i] = (hmm_lane){NULL, hp->bwd[i], mismatch_row(hp, i, m + 1), hp->p_mismatch, NULL, p_recomb[m + 1], n_states[i]};
        }
        hmm_bwd_lanes(lanes, 2);
        if (slot[m] < 0) continue;
        for (int i = 0; i < 2; ++i) {
            float *row = probs[i] + (size_t)slot[m] * (size_t)hp->max_states;
            const float *bwd = hp->bwd[i];
            float sum = 0.0f;
            for (int j = 0; j < n_states[i]; ++j) {
                row[j] *= bwd[j];
                sum += row[j];
            }
            for (int j = 0; j < n_states[i]; ++j) row[j] /= sum;
        }
    }
}

static void hmm_state_probs_run(hmm_state_probs *hp, int sample, const int *slot, int *const states[2], float *const probs[2], int n_states[2]) {
    for (int i = 0; i < 2; ++i) {
        n_states[i] = low_freq_phase_states_ibs_states(&hp->states, (sample << 1) | i, slot, states[i], hp->mismatch[i]);
    }
    run_fwd(hp, slot, probs, n_states);
    run_bwd(hp, slot, probs, n_states);
}

/* Stage2Baum */
typedef struct {
    const phase_data *pd;
    const fixed_phase_data *fpd;
    stage2_haps *s2;
    hmm_state_probs hp;
    int n_states[2];
    int *slot;          /* per stage-1 marker: its row in states and probs, or -1 */
    int_list kept;      /* the stage-1 markers that have a row */
    int_list todo;      /* the sample's stage-2 markers with a heterozygous or missing genotype */
    int rows_cap;
    int *states[2];     /* [row][state] */
    float *probs[2];
    int n_targ_haps;
    uint64_t digest[2]; /* trace seam T4b */
    uint64_t decisions; /* trace seam T4d: the sample's allele probabilities at stage-2 markers */
} stage2_baum;

static void stage2_baum_init(stage2_baum *b, const low_freq_phase_ibs *ibs, stage2_haps *s2) {
    const phase_data *pd = ibs->pd;
    b->pd = pd;
    b->fpd = pd->fpd;
    b->s2 = s2;
    hmm_state_probs_init(&b->hp, ibs);
    int n = b->fpd->n_stage1;
    b->slot = util_malloc((size_t)n * sizeof *b->slot);
    for (int m = 0; m < n; ++m) b->slot[m] = -1;
    b->kept = b->todo = (int_list){0};
    b->rows_cap = 0;
    for (int i = 0; i < 2; ++i) {
        b->states[i] = NULL;
        b->probs[i] = NULL;
    }
    b->n_targ_haps = fpd_n_targ_haps(b->fpd);
}

static void stage2_baum_free(stage2_baum *b) {
    for (int i = 0; i < 2; ++i) {
        free(b->states[i]);
        free(b->probs[i]);
    }
    free(b->slot);
    free(b->kept.v);
    free(b->todo.v);
    hmm_state_probs_free(&b->hp);
}

static void keep(stage2_baum *b, int stage1_marker) {
    if (b->slot[stage1_marker] >= 0) return;
    b->slot[stage1_marker] = b->kept.n;
    int_list_add(&b->kept, stage1_marker);
}

/* The two stage-1 markers whose state probabilities are interpolated at
 * target marker m. */
static int marker_a(const fixed_phase_data *fpd, int m) {
    return fpd->prev_stage1_marker[m];
}

static int marker_b(const fixed_phase_data *fpd, int m) {
    int mkr_a = fpd->prev_stage1_marker[m];
    return mkr_a + 1 < fpd->n_stage1 - 1 ? mkr_a + 1 : fpd->n_stage1 - 1;
}

static void ensure_rows(stage2_baum *b) {
    if (b->kept.n <= b->rows_cap) return;
    b->rows_cap = b->kept.n + b->kept.n / 2;
    size_t n = (size_t)b->rows_cap * (size_t)b->hp.max_states;
    for (int i = 0; i < 2; ++i) {
        free(b->states[i]);
        free(b->probs[i]);
        b->states[i] = util_malloc(n * sizeof *b->states[i]);
        b->probs[i] = util_malloc(n * sizeof *b->probs[i]);
    }
}

static const int *states_row(const stage2_baum *b, int hap_bit, int stage1_marker) {
    return b->states[hap_bit] + (size_t)b->slot[stage1_marker] * (size_t)b->hp.max_states;
}

static const float *probs_row(const stage2_baum *b, int hap_bit, int stage1_marker) {
    return b->probs[hap_bit] + (size_t)b->slot[stage1_marker] * (size_t)b->hp.max_states;
}

/* The spliced allele of a target haplotype, or a reference haplotype's allele,
 * at target marker m. */
static int allele(const stage2_baum *b, int m, int hap) {
    const window *w = b->fpd->win;
    if (hap < b->n_targ_haps) return fpd_spliced_allele(b->fpd, m, hap);
    return ref_gt_rec_get(w->ref[w->indices.targ_marker_to_marker[m]], hap - b->n_targ_haps);
}

static bool is_low_freq(const fixed_phase_data *fpd, int m, int allele) {
    return fpd->carriers[m].allele[allele].kind != CARRIERS_HIGH_FREQ;
}

/* Probabilities that the haplotype carries a1 or a2 at a heterozygote: only
 * states whose own genotype is homozygous, or shares exactly one of the
 * heterozygote's rare alleles, contribute. */
static void unscaled_al_probs(const stage2_baum *b, int m, int hap_bit, int a1, int a2, float *al_probs, int n_alleles) {
    const fixed_phase_data *fpd = b->fpd;
    for (int a = 0; a < n_alleles; ++a) al_probs[a] = 0.0f;
    bool rare1 = is_low_freq(fpd, m, a1);
    bool rare2 = is_low_freq(fpd, m, a2);
    const int *states_a = states_row(b, hap_bit, marker_a(fpd, m));
    const float *probs_a = probs_row(b, hap_bit, marker_a(fpd, m));
    const float *probs_b = probs_row(b, hap_bit, marker_b(fpd, m));
    for (int j = 0; j < b->n_states[hap_bit]; ++j) {
        int hap = states_a[j];
        int b1 = allele(b, m, hap);
        int b2 = allele(b, m, hap ^ 1);
        if (b1 < 0 || b2 < 0) continue;
        float wt = fpd->prev_stage1_wt[m];
        float prob = wt * probs_a[j] + (1.0f - wt) * probs_b[j];
        if (b1 == b2) {
            al_probs[b1] += prob;
        } else {
            bool match1 = rare1 && (a1 == b1 || a1 == b2);
            bool match2 = rare2 && (a2 == b1 || a2 == b2);
            if (match1 ^ match2) al_probs[match1 ? a1 : a2] += prob;
        }
    }
}

/* The most probable allele at a missing genotype. The 0.55 and 0.45 weights
 * are applied in double and the sums rounded to float, as in Java. */
static int impute_allele(const stage2_baum *b, int m, int hap_bit, float *al_probs, int n_alleles) {
    const fixed_phase_data *fpd = b->fpd;
    for (int a = 0; a < n_alleles; ++a) al_probs[a] = 0.0f;
    const int *states_a = states_row(b, hap_bit, marker_a(fpd, m));
    const float *probs_a = probs_row(b, hap_bit, marker_a(fpd, m));
    const float *probs_b = probs_row(b, hap_bit, marker_b(fpd, m));
    for (int j = 0; j < b->n_states[hap_bit]; ++j) {
        float wt = fpd->prev_stage1_wt[m];
        float prob = wt * probs_a[j] + (1.0f - wt) * probs_b[j];
        int hap = states_a[j];
        int b1 = allele(b, m, hap);
        int b2 = allele(b, m, hap ^ 1);
        if (b1 < 0 || b2 < 0) continue;
        if (b1 == b2 || hap >= b->n_targ_haps) {
            al_probs[b1] += prob;
        } else {
            bool rare_b1 = is_low_freq(fpd, m, b1);
            bool rare_b2 = is_low_freq(fpd, m, b2);
            if (rare_b1 ^ rare_b2) {
                al_probs[b1] += (rare_b1 ? 0.55 : 0.45) * prob;
                al_probs[b2] += (rare_b1 ? 0.45 : 0.55) * prob;
            } else {
                al_probs[b1] += 0.5 * prob;
                al_probs[b2] += 0.5 * prob;
            }
        }
    }
    int max_index = 0;
    for (int a = 1; a < n_alleles; ++a) {
        if (al_probs[a] > al_probs[max_index]) max_index = a;
    }
    return max_index;
}

/* Stage2Haps.setPhasedGT */
static void set_phased_gt(stage2_haps *s2, int m, int sample, int a1, int a2) {
    const fixed_phase_data *fpd = s2->pd->fpd;
    bool low1 = is_low_freq(fpd, m, a1);
    bool low2 = is_low_freq(fpd, m, a2);
    if (!low1 && !low2) return;
    pthread_mutex_lock(&s2->locks[m]);
    if (low1) int_list_add(&s2->carriers[m][a1], sample << 1);
    if (low2) int_list_add(&s2->carriers[m][a2], (sample << 1) | 1);
    pthread_mutex_unlock(&s2->locks[m]);
}

static void fold_probs(stage2_baum *b, const float *al_probs, int n_alleles) {
    for (int a = 0; a < n_alleles; ++a) b->decisions = trace_fold(b->decisions, jnum_float_bits(al_probs[a]));
}

/* Stage-2 markers [start, end): a homozygous genotype needs no state
 * probabilities and is recorded at once. A heterozygous or missing genotype
 * waits in todo, and the two stage-1 markers it reads get a row. */
static void scan_interval(stage2_baum *b, int sample, int start, int end) {
    const fixed_phase_data *fpd = b->fpd;
    int hap1 = sample << 1;
    int hap2 = hap1 | 1;
    for (int m = start; m < end; ++m) {
        int a1 = fpd_spliced_allele(fpd, m, hap1);
        int a2 = fpd_spliced_allele(fpd, m, hap2);
        if (a1 >= 0 && a1 == a2) {
            set_phased_gt(b->s2, m, sample, a1, a2);
        } else {
            int_list_add(&b->todo, m);
            keep(b, marker_a(fpd, m));
            keep(b, marker_b(fpd, m));
        }
    }
}

static void impute_marker(stage2_baum *b, jrandom *r, int sample, int m) {
    int n_alleles = marker_n_alleles(&b->fpd->win->targ[m]->marker);
    float *al1 = util_malloc((size_t)n_alleles * sizeof *al1);
    float *al2 = util_malloc((size_t)n_alleles * sizeof *al2);
    int a1 = fpd_spliced_allele(b->fpd, m, sample << 1);
    int a2 = fpd_spliced_allele(b->fpd, m, (sample << 1) | 1);
    if (a1 >= 0 && a2 >= 0) {
        unscaled_al_probs(b, m, 0, a1, a2, al1, n_alleles);
        unscaled_al_probs(b, m, 1, a1, a2, al2, n_alleles);
        float p1 = al1[a1] * al2[a2];
        float p2 = al1[a2] * al2[a1];
        b->decisions = trace_fold(b->decisions, (uint32_t)m);
        b->decisions = trace_fold(b->decisions, jnum_float_bits(p1));
        b->decisions = trace_fold(b->decisions, jnum_float_bits(p2));
        if (p1 < p2 || (p1 == p2 && jrandom_next_boolean(r))) {
            int tmp = a1;
            a1 = a2;
            a2 = tmp;
        }
    } else {
        a1 = impute_allele(b, m, 0, al1, n_alleles);
        a2 = impute_allele(b, m, 1, al2, n_alleles);
        b->decisions = trace_fold(b->decisions, (uint32_t)m);
        fold_probs(b, al1, n_alleles);
        fold_probs(b, al2, n_alleles);
    }
    set_phased_gt(b->s2, m, sample, a1, a2);
    free(al1);
    free(al2);
}

static uint64_t digest(const stage2_baum *b, int hap_bit) {
    uint64_t h = TRACE_FNV_BASIS;
    for (int m = 0; m < b->fpd->n_stage1; ++m) {
        const int *states = states_row(b, hap_bit, m);
        const float *probs = probs_row(b, hap_bit, m);
        for (int j = 0; j < b->n_states[hap_bit]; ++j) {
            h = trace_fold(h, (uint32_t)states[j]);
            h = trace_fold(h, jnum_float_bits(probs[j]));
        }
    }
    return h;
}

static void stage2_baum_phase(stage2_baum *b, int sample) {
    const fixed_phase_data *fpd = b->fpd;
    jrandom r;
    jrandom_init(&r, jrandom_seed_plus(phase_data_seed(b->pd), sample));
    b->todo.n = 0;
    if (trace_on()) {
        /* The T4b digest covers every stage-1 marker. */
        for (int m = 0; m < fpd->n_stage1; ++m) keep(b, m);
    }
    int start = 0;
    for (int j = 0; j < fpd->n_stage1; ++j) {
        int end = fpd->stage1_to2[j];
        scan_interval(b, sample, start, end);
        start = end + 1;
    }
    scan_interval(b, sample, start, fpd->n_markers);
    b->n_states[0] = b->n_states[1] = 0;
    b->digest[0] = b->digest[1] = 0;
    b->decisions = TRACE_FNV_BASIS;
    if (b->kept.n == 0) return;
    ensure_rows(b);
    hmm_state_probs_run(&b->hp, sample, b->slot, b->states, b->probs, b->n_states);
    if (trace_on()) {
        b->digest[0] = digest(b, 0);
        b->digest[1] = digest(b, 1);
    }
    for (int k = 0; k < b->todo.n; ++k) impute_marker(b, &r, sample, b->todo.v[k]);
    for (int k = 0; k < b->kept.n; ++k) b->slot[b->kept.v[k]] = -1;
    b->kept.n = 0;
}

static void stage2_haps_init(stage2_haps *s2, const phase_data *pd) {
    const fixed_phase_data *fpd = pd->fpd;
    s2->pd = pd;
    s2->carriers = util_malloc((size_t)fpd->n_markers * sizeof *s2->carriers);
    s2->locks = util_malloc((size_t)fpd->n_markers * sizeof *s2->locks);
    for (int m = 0; m < fpd->n_markers; ++m) {
        s2->carriers[m] = NULL;
        pthread_mutex_init(&s2->locks[m], NULL);
    }
    int start = 0;
    for (int j = 0; j <= fpd->n_stage1; ++j) {
        int end = j < fpd->n_stage1 ? fpd->stage1_to2[j] : fpd->n_markers;
        for (int m = start; m < end; ++m) {
            int n_alleles = fpd->carriers[m].n_alleles;
            s2->carriers[m] = util_malloc((size_t)n_alleles * sizeof **s2->carriers);
            for (int a = 0; a < n_alleles; ++a) s2->carriers[m][a] = (int_list){0};
        }
        start = end + 1;
    }
}

void stage2_haps_free(stage2_haps *s2) {
    if (s2->pd == NULL) return;
    for (int m = 0; m < s2->pd->fpd->n_markers; ++m) {
        if (s2->carriers[m] == NULL) continue;
        for (int a = 0; a < s2->pd->fpd->carriers[m].n_alleles; ++a) free(s2->carriers[m][a].v);
        free(s2->carriers[m]);
    }
    for (int m = 0; m < s2->pd->fpd->n_markers; ++m) pthread_mutex_destroy(&s2->locks[m]);
    free(s2->locks);
    free(s2->carriers);
    free(s2->major);
}

/* Stage2Haps.stage2Rec: the one high-frequency allele, or when missing
 * genotypes leave every allele rare, the first allele with the most carriers. */
static void set_major_alleles(stage2_haps *s2) {
    const fixed_phase_data *fpd = s2->pd->fpd;
    s2->major = util_malloc((size_t)fpd->n_markers * sizeof *s2->major);
    for (int m = 0; m < fpd->n_markers; ++m) {
        int major = -1;
        if (s2->carriers[m] != NULL) {
            int n_alleles = fpd->carriers[m].n_alleles;
            for (int a = 0; a < n_alleles; ++a) {
                if (!is_low_freq(fpd, m, a)) major = a;
            }
            if (major == -1) {
                major = 0;
                for (int a = 1; a < n_alleles; ++a) {
                    if (s2->carriers[m][a].n > s2->carriers[m][major].n) major = a;
                }
            }
        }
        s2->major[m] = major;
    }
}

static bool contains(const int_list *l, int x) {
    int lo = 0, hi = l->n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) >> 1;
        if (l->v[mid] < x) lo = mid + 1;
        else if (l->v[mid] > x) hi = mid - 1;
        else return true;
    }
    return false;
}

int stage2_haps_allele(const void *ctx, int m, int hap) {
    const stage2_haps *s2 = ctx;
    const fixed_phase_data *fpd = s2->pd->fpd;
    if (s2->carriers[m] == NULL) {
        const sample_phase *sp = &s2->pd->phase[hap >> 1];
        int j = fpd->prev_stage1_marker[m];
        return (hap & 1) == 0 ? sample_phase_allele1(sp, j) : sample_phase_allele2(sp, j);
    }
    for (int a = 0; a < fpd->carriers[m].n_alleles; ++a) {
        if (a != s2->major[m] && contains(&s2->carriers[m][a], hap)) return a;
    }
    return s2->major[m];
}

/* Trace seam T4c: at each stage-2 marker, the haplotypes carrying each
 * allele, "M" for the major allele that Stage2Haps.stage2Rec leaves implicit. */
static void trace_stage2_haps(const stage2_haps *s2) {
    const fixed_phase_data *fpd = s2->pd->fpd;
    kstring_t s = {0, 0, NULL};
    for (int m = 0; m < fpd->n_markers; ++m) {
        if (s2->carriers[m] == NULL) continue;
        s.l = 0;
        ksprintf(&s, "m\t%d", m);
        for (int a = 0; a < fpd->carriers[m].n_alleles; ++a) {
            kputc('\t', &s);
            if (a == s2->major[m]) {
                kputc('M', &s);
                continue;
            }
            for (int j = 0; j < s2->carriers[m][a].n; ++j) {
                if (j > 0) kputc(',', &s);
                kputw(s2->carriers[m][a].v[j], &s);
            }
        }
        trace_line("T4c", "%s", s.s);
    }
    free(s.s);
}

typedef struct {
    stage2_baum baum;
    int *n_states;      /* per target haplotype, for trace seam T4b */
    uint64_t *digest;   /* per target haplotype, T4b */
    uint64_t *decisions; /* per sample, T4d */
} stage2_worker;

static void stage2_task(void *worker, int sample) {
    stage2_worker *sw = worker;
    stage2_baum_phase(&sw->baum, sample);
    for (int i = 0; i < 2; ++i) {
        sw->n_states[(sample << 1) | i] = sw->baum.n_states[i];
        sw->digest[(sample << 1) | i] = sw->baum.digest[i];
    }
    sw->decisions[sample] = sw->baum.decisions;
}

void phase_ls_run_stage2(stage2_haps *s2, const phase_data *pd) {
    low_freq_phase_ibs ibs;
    low_freq_phase_ibs_init(&ibs, pd);
    stage2_haps_init(s2, pd);
    int n_haps = pd->n_samples << 1;
    int *n_states = util_malloc((size_t)n_haps * sizeof *n_states);
    uint64_t *digest = util_malloc((size_t)n_haps * sizeof *digest);
    uint64_t *decisions = util_malloc((size_t)pd->n_samples * sizeof *decisions);
    int n_threads = parallel_threads(pd->par->nthreads, pd->n_samples);
    stage2_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) {
        stage2_baum_init(&workers[t].baum, &ibs, s2);
        workers[t].n_states = n_states;
        workers[t].digest = digest;
        workers[t].decisions = decisions;
    }
    /* Each sample adds only its own haplotypes; the lists are sorted after,
     * as Java sorts them when read. */
    parallel_for(n_threads, pd->n_samples, workers, sizeof *workers, stage2_task);
    for (int t = 0; t < n_threads; ++t) stage2_baum_free(&workers[t].baum);
    free(workers);
    low_freq_phase_ibs_free(&ibs);
    for (int m = 0; m < pd->fpd->n_markers; ++m) {
        if (s2->carriers[m] == NULL) continue;
        for (int a = 0; a < pd->fpd->carriers[m].n_alleles; ++a) {
            int_list *l = &s2->carriers[m][a];
            if (l->n > 1) qsort(l->v, (size_t)l->n, sizeof *l->v, util_compare_ints);
        }
    }
    set_major_alleles(s2);
    if (trace_on()) {
        for (int s = 0; s < pd->n_samples; ++s) {
            for (int i = 0; i < 2; ++i) {
                int h = (s << 1) | i;
                trace_line("T4b", "hap\t%d\t%d\t%" PRIx64, h, n_states[h], digest[h]);
            }
            trace_line("T4d", "sample\t%d\t%" PRIx64, s, decisions[s]);
        }
        trace_stage2_haps(s2);
    }
    free(n_states);
    free(digest);
    free(decisions);
}
