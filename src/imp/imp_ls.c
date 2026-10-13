/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) imp/ImpLS.java, imp/ImpStates.java,
 * imp/ImpLSBaum.java, imp/StateProbsFactory.java and imp/StateProbs.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "imp/imp_ls.h"

#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>

#include "beagleutil/comp_hap_queue.h"
#include "blbutil/int_list.h"
#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "imp/imp_ibs.h"
#include "jcompat/jnum.h"
#include "jcompat/jrandom.h"
#include "phase/hmm_updater.h"

/* ImpStates: up to imp-states composite reference haplotypes built from the
 * target haplotype's IBS sets, as in BasicPhaseStates. */
typedef struct {
    const imp_data *id;
    const imp_ibs *ibs;
    int max_states;
    comp_hap_tracker t;
    int_list *comp_hap_hap;
    int_list *comp_hap_end;
    int *list_index;
    int *comp_hap_to_hap;
    int *ev_head;      /* per cluster: first state whose segment ends there, or -1 */
    int *ev_next;      /* per state: next state in the same cluster's list */
    int n_used;        /* lists to clear before the next haplotype */
} imp_states;

static void imp_states_init(imp_states *st, const imp_data *id, const imp_ibs *ibs, int max_states) {
    size_t n = (size_t)max_states;
    st->id = id;
    st->ibs = ibs;
    st->max_states = max_states;
    comp_hap_tracker_init(&st->t, max_states, id->n_haps);
    st->comp_hap_hap = util_malloc(n * sizeof *st->comp_hap_hap);
    st->comp_hap_end = util_malloc(n * sizeof *st->comp_hap_end);
    for (int j = 0; j < max_states; ++j) st->comp_hap_hap[j] = st->comp_hap_end[j] = (int_list){0};
    st->list_index = util_malloc(n * sizeof *st->list_index);
    st->comp_hap_to_hap = util_malloc(n * sizeof *st->comp_hap_to_hap);
    st->ev_head = util_malloc((size_t)id->n_clusters * sizeof *st->ev_head);
    for (int m = 0; m < id->n_clusters; ++m) st->ev_head[m] = -1;
    st->ev_next = util_malloc(n * sizeof *st->ev_next);
    st->n_used = 0;
}

static void imp_states_free(imp_states *st) {
    for (int j = 0; j < st->max_states; ++j) {
        free(st->comp_hap_hap[j].v);
        free(st->comp_hap_end[j].v);
    }
    free(st->comp_hap_hap);
    free(st->comp_hap_end);
    free(st->list_index);
    free(st->comp_hap_to_hap);
    free(st->ev_head);
    free(st->ev_next);
    comp_hap_tracker_free(&st->t);
}

/* ImpStates has no staleness rule: a segment is replaced only when the queue
 * is full, so min_steps is INT_MAX. */
static void update_fields(imp_states *st, int hap, int step) {
    comp_hap_change c = comp_hap_tracker_observe(&st->t, hap, step, INT_MAX);
    if (c.index >= 0) {
        int_list_add(&st->comp_hap_hap[c.index], hap);
        if (c.old_hap >= 0) int_list_add(&st->comp_hap_end[c.index], st->ibs->step_starts[c.end_step]);
    }
}

/* Used when no IBS set has a haplotype. Java seeds with the target haplotype
 * index and re-draws while the reference haplotype drawn equals it. */
static void fill_q_with_random_haps(imp_states *st, int hap) {
    int n_ref_haps = st->id->n_ref_haps;
    int n_states = n_ref_haps < st->max_states ? n_ref_haps : st->max_states;
    jrandom r;
    jrandom_init(&r, hap);
    for (int j = 0; j < n_states; ++j) {
        int h = jrandom_next_int_bound(&r, n_ref_haps);
        while (h == hap) h = jrandom_next_int_bound(&r, n_ref_haps);
        int_list_add(&st->comp_hap_hap[comp_hap_tracker_seed(&st->t, h)], h);
    }
}

/* Files state j under the cluster where its current segment ends. A state
 * whose segment ended before cluster m never moves again. */
static void segments_schedule(imp_states *st, int j, int end, int m) {
    if (end >= m && end < st->id->n_clusters) {
        st->ev_next[j] = st->ev_head[end];
        st->ev_head[end] = j;
    }
}

/* Each state's reference haplotype at cluster m. Callers visit every
 * cluster from 0 in increasing order; m == 0 restarts the walk. */
static const int *segments_at(imp_states *st, int m) {
    if (m == 0) {
        for (int c = 0; c < st->id->n_clusters; ++c) st->ev_head[c] = -1;
        for (int j = 0; j < st->n_used; ++j) {
            st->list_index[j] = 0;
            st->comp_hap_to_hap[j] = st->comp_hap_hap[j].v[0];
            segments_schedule(st, j, st->comp_hap_end[j].v[0], 0);
        }
    }
    int j = st->ev_head[m];
    st->ev_head[m] = -1;
    while (j >= 0) {
        int next = st->ev_next[j];
        int i = ++st->list_index[j];
        st->comp_hap_to_hap[j] = st->comp_hap_hap[j].v[i];
        segments_schedule(st, j, st->comp_hap_end[j].v[i], m + 1);
        j = next;
    }
    return st->comp_hap_to_hap;
}

/* ImpStates.ibsStates: builds each state's segment list, and sets
 * mismatch[c][j] to whether state j's sequence at cluster c differs from
 * the target's. When digest is non-NULL, folds each (hap, match) pair into
 * it for the T5c trace. */
static int imp_states_ibs_states(imp_states *st, int targ_hap, uint8_t **mismatch, uint64_t *digest) {
    const imp_data *id = st->id;
    for (int j = 0; j < st->n_used; ++j) st->comp_hap_hap[j].n = st->comp_hap_end[j].n = 0;
    comp_hap_tracker_clear(&st->t);
    for (int j = 0; j < st->ibs->n_steps; ++j) {
        const ibs_set *set = st->ibs->ibs[j][targ_hap];
        for (int k = 0; k < set->n; ++k) update_fields(st, set->haps[k], j);
    }
    if (comp_hap_tracker_size(&st->t) == 0) fill_q_with_random_haps(st, targ_hap);
    int n_comp_haps = comp_hap_tracker_size(&st->t);
    st->n_used = n_comp_haps;
    for (int j = 0; j < n_comp_haps; ++j) int_list_add(&st->comp_hap_end[j], id->n_clusters);
    uint64_t h = TRACE_FNV_BASIS;
    for (int m = 0; m < id->n_clusters; ++m) {
        const int *hap = segments_at(st, m);
        /* imp_data_seq with the cluster's map read once */
        const hap_seq_map *hs = &id->hap_to_seq[m];
        int targ_allele = hs->targ_seq[targ_hap];
        const int *ref_seq = hs->ref_seq;
        uint8_t *mis = mismatch[m];
        if (hs->group != NULL) {
            const int *group_seq = hs->group->hap_to_seq;
            for (int j = 0; j < n_comp_haps; ++j) mis[j] = ref_seq[group_seq[hap[j]]] != targ_allele;
        } else {
            for (int j = 0; j < n_comp_haps; ++j) mis[j] = ref_seq[hap[j]] != targ_allele;
        }
        if (digest != NULL) {
            for (int j = 0; j < n_comp_haps; ++j) h = trace_fold(trace_fold(h, (uint32_t)hap[j]), mis[j] ? 0u : 1u);
        }
    }
    if (digest != NULL) *digest = h;
    return n_comp_haps;
}

/* ImpLSBaum */
typedef struct {
    const imp_data *id;
    imp_states states;
    uint8_t **mismatch;
    float **fwd_val;
    float *bwd_val;
    int *kept;         /* keep_state_probs: the kept states of one cluster */
    kept_state *buf;   /* keep_state_probs: the kept states of one haplotype */
    size_t buf_cap;
} imp_ls_baum;

static void set_fwd_values(imp_ls_baum *b, int n_haps) {
    const imp_data *id = b->id;
    float em[2] = {1.0f - id->err_prob[0], id->err_prob[0]};
    float sum = 0.0f;
    for (int j = 0; j < n_haps; ++j) sum += b->fwd_val[0][j] = em[b->mismatch[0][j]];
    for (int m = 1; m < id->n_clusters; ++m) {
        em[0] = 1.0f - id->err_prob[m];
        em[1] = id->err_prob[m];
        hmm_lane lane = {b->fwd_val[m - 1], b->fwd_val[m], b->mismatch[m], em, &sum, id->p_recomb[m], n_haps};
        hmm_fwd_lanes(&lane, 1);
    }
}

static float set_bwd_value(imp_ls_baum *b, int m, int n_states, float last_sum) {
    const imp_data *id = b->id;
    float p_recomb = m + 1 < id->n_clusters ? id->p_recomb[m + 1] : 0.0f;
    float p_err = id->err_prob[m];
    float p_no_err = 1.0f - p_err;
    float scale = (1.0f - p_recomb) / last_sum;
    float shift = p_recomb / n_states;
    float bwd_val_sum = 0.0f;
    float state_sum = 0.0f;
    float *fwd = b->fwd_val[m];
    for (int j = 0; j < n_states; ++j) {
        b->bwd_val[j] = scale * b->bwd_val[j] + shift;
        fwd[j] *= b->bwd_val[j];
        state_sum += fwd[j];
        float em = b->mismatch[m][j] ? p_err : p_no_err;
        b->bwd_val[j] *= em;
        bwd_val_sum += b->bwd_val[j];
    }
    for (int j = 0; j < n_states; ++j) fwd[j] /= state_sum;
    return bwd_val_sum;
}

/* StateProbsFactory.stateProbs: the kept states of haplotype targ_hap, gathered
 * in b->buf and then copied once into the haplotype's own block. */
static void keep_state_probs(imp_ls_baum *b, int n_states, int targ_hap, state_probs *out) {
    int n = b->id->n_clusters;
    float threshold = 0.9999f / n_states;
    if (!(threshold < 0.005f)) threshold = 0.005f;
    cluster_probs *cps = &out->cluster[(size_t)targ_hap * (size_t)n];
    size_t total = 0;
    for (int m = 0; m < n; ++m) {
        int m_p1 = m < n - 1 ? m + 1 : m;
        const int *hap = segments_at(&b->states, m);
        const float *fwd = b->fwd_val[m], *fwd_p1 = b->fwd_val[m_p1];
        int kept = 0;
        for (int j = 0; j < n_states; ++j) {  /* branchless: few states pass */
            b->kept[kept] = j;
            kept += (fwd[j] > threshold) | (fwd_p1[j] > threshold);
        }
        if (total + (size_t)kept > b->buf_cap) {
            b->buf_cap = 2 * (total + (size_t)kept);
            b->buf = util_realloc(b->buf, b->buf_cap * sizeof *b->buf);
        }
        cps[m].n = kept;
        for (int k = 0; k < kept; ++k) {
            int j = b->kept[k];
            b->buf[total + (size_t)k] = (kept_state){hap[j], fwd[j], fwd_p1[j]};
        }
        total += (size_t)kept;
    }
    kept_state *block = util_malloc(total * sizeof *block);
    /* With no kept states buf may still be NULL, which memcpy may not receive. */
    if (total > 0) memcpy(block, b->buf, total * sizeof *block);
    for (int m = 0; m < n; ++m) {
        cps[m].s = block;
        block += cps[m].n;
    }
}

static int impute(imp_ls_baum *b, int targ_hap, state_probs *out, uint64_t *states_digest) {
    const imp_data *id = b->id;
    int n_states = imp_states_ibs_states(&b->states, targ_hap, b->mismatch, trace_on() ? states_digest : NULL);
    set_fwd_values(b, n_states);
    for (int j = 0; j < n_states; ++j) b->bwd_val[j] = 1.0f / n_states;
    float last_sum = 1.0f;
    for (int m = id->n_clusters - 1; m >= 0; --m) last_sum = set_bwd_value(b, m, n_states, last_sum);
    keep_state_probs(b, n_states, targ_hap, out);
    return n_states;
}

static uint64_t probs_digest(const state_probs *sp, int targ_hap) {
    uint64_t h = TRACE_FNV_BASIS;
    for (int m = 0; m < sp->n_clusters; ++m) {
        const cluster_probs *cp = state_probs_at(sp, targ_hap, m);
        h = trace_fold(h, (uint32_t)cp->n);
        for (int k = 0; k < cp->n; ++k) {
            h = trace_fold(h, (uint32_t)cp->s[k].ref_hap);
            h = trace_fold(h, jnum_float_bits(cp->s[k].prob));
            h = trace_fold(h, jnum_float_bits(cp->s[k].prob_p1));
        }
    }
    return h;
}

static void baum_init(imp_ls_baum *b, const imp_data *id, const imp_ibs *ibs, int max_states) {
    b->id = id;
    imp_states_init(&b->states, id, ibs, max_states);
    b->mismatch = util_malloc((size_t)id->n_clusters * sizeof *b->mismatch);
    b->fwd_val = util_malloc((size_t)id->n_clusters * sizeof *b->fwd_val);
    for (int m = 0; m < id->n_clusters; ++m) {
        b->mismatch[m] = util_malloc((size_t)max_states * sizeof **b->mismatch);
        b->fwd_val[m] = util_malloc((size_t)max_states * sizeof **b->fwd_val);
    }
    b->bwd_val = util_malloc((size_t)max_states * sizeof *b->bwd_val);
    b->kept = util_malloc((size_t)max_states * sizeof *b->kept);
    b->buf = NULL;
    b->buf_cap = 0;
}

static void baum_free(imp_ls_baum *b) {
    for (int m = 0; m < b->id->n_clusters; ++m) {
        free(b->mismatch[m]);
        free(b->fwd_val[m]);
    }
    free(b->mismatch);
    free(b->fwd_val);
    free(b->bwd_val);
    free(b->kept);
    free(b->buf);
    imp_states_free(&b->states);
}

typedef struct {
    imp_ls_baum baum;
    state_probs *sp;
    int *n_states;            /* per haplotype, for trace seam T5c */
    uint64_t *states_digest;  /* per haplotype, T5c */
} imp_worker;

static void imp_task(void *worker, int h) {
    imp_worker *iw = worker;
    iw->n_states[h] = impute(&iw->baum, h, iw->sp, &iw->states_digest[h]);
}

static state_probs *state_probs_create(int n_targ_haps, int n_clusters) {
    state_probs *sp = util_malloc(sizeof *sp);
    sp->n_targ_haps = n_targ_haps;
    sp->n_clusters = n_clusters;
    sp->cluster = util_malloc((size_t)n_targ_haps * (size_t)n_clusters * sizeof *sp->cluster);
    return sp;
}

state_probs *imp_ls_state_probs(const imp_data *id, const par *p) {
    imp_ibs ibs;
    imp_ibs_init(&ibs, id, p);
    int n = id->n_targ_haps;
    state_probs *sp = state_probs_create(n, id->n_clusters);
    int *n_states = util_malloc((size_t)n * sizeof *n_states);
    uint64_t *states_digest = util_malloc((size_t)n * sizeof *states_digest);
    int n_threads = parallel_threads(p->nthreads, n);
    imp_worker *workers = util_malloc((size_t)n_threads * sizeof *workers);
    for (int t = 0; t < n_threads; ++t) {
        baum_init(&workers[t].baum, id, &ibs, p->imp_states);
        workers[t].sp = sp;
        workers[t].n_states = n_states;
        workers[t].states_digest = states_digest;
    }
    /* Each haplotype's result depends only on the haplotype. */
    parallel_for(n_threads, n, workers, sizeof *workers, imp_task);
    for (int t = 0; t < n_threads; ++t) baum_free(&workers[t].baum);
    free(workers);
    imp_ibs_free(&ibs);
    if (trace_on()) {
        for (int h = 0; h < n; ++h) {
            trace_line("T5c", "hap\t%d\t%d\t%" PRIx64 "\t%" PRIx64, h, n_states[h], states_digest[h],
                    probs_digest(sp, h));
        }
    }
    free(n_states);
    free(states_digest);
    return sp;
}

void state_probs_free(state_probs *sp) {
    for (int h = 0; h < sp->n_targ_haps; ++h) free(state_probs_at(sp, h, 0)->s);
    free(sp->cluster);
    free(sp);
}
