/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) imp/ImputedVcfWriter.java,
 * imp/RefHapHash.java and WindowWriter.printImputed; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "imp/imputed_writer.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "blbutil/int_list.h"
#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jrandom.h"

/* RefHapHash: the reference haplotypes in a cluster's kept states, each
 * with a hash of its ALT alleles over all of the cluster's reference markers.
 * The hash decides which states set_al_probs sums together, so every piece of
 * a cluster must use the whole cluster's hash. */
typedef struct {
    int n;
    int *i2hap;
    int *hap2i;           /* per reference haplotype; its index in i2hap */
    int32_t *i2hash;
} ref_hap_hash;

/* Per haplotype in a ref_hap_hash, its (marker offset, allele) pairs for
 * the ALT alleles at reference markers [start, end), one piece of a cluster. */
typedef struct {
    int_list *alt_alleles;
    int start;
    int end;
} alt_lists;

/* The alleles at rec of the increasing haplotypes haps[0, n): one pass over
 * each non-major allele's increasing haplotype list instead of a search per
 * haplotype. */
static void hap_alleles(const ref_gt_rec *rec, const int *haps, int n, int *out) {
    if (rec->kind == REF_HAP) {
        for (int i = 0; i < n; ++i) out[i] = rec->seq_to_allele[rec->group->hap_to_seq[haps[i]]];
        return;
    }
    for (int i = 0; i < n; ++i) out[i] = rec->major_allele;
    int n_alleles = marker_n_alleles(&rec->marker);
    for (int a = 0; a < n_alleles; ++a) {
        if (a == rec->major_allele) continue;
        const int *v = rec->hap_lists[a];
        int len = rec->hap_list_len[a];
        for (int i = 0, k = 0; i < n && k < len;) {
            if (haps[i] < v[k]) ++i;
            else if (haps[i] > v[k]) ++k;
            else {
                out[i++] = a;
                ++k;
            }
        }
    }
}

/* Adds the ALT alleles at reference markers [start, end) of rh's haplotypes
 * to i2hash, seeded with start, if i2hash is not NULL, and to al's lists if al
 * is not NULL. */
static void scan_alleles(const ref_hap_hash *rh, const imp_data *id, int start, int end, int32_t *i2hash, alt_lists *al) {
    int n = rh->n;
    /* Java's two update paths draw the same numbers and give the same
     * hashes and lists; this is RefHapHash.standardUpdate. */
    jrandom r;
    jrandom_init(&r, start);
    int32_t *allele_hash = NULL;
    int *alleles = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *alleles);
    for (int m = start; m < end; ++m) {
        const ref_gt_rec *rec = id->w->ref[m];
        int n_alleles = marker_n_alleles(&rec->marker);
        if (i2hash != NULL) {
            allele_hash = util_realloc(allele_hash, (size_t)n_alleles * sizeof *allele_hash);
            allele_hash[0] = 0;
            for (int a = 1; a < n_alleles; ++a) allele_hash[a] = jrandom_next_int(&r);
        }
        hap_alleles(rec, rh->i2hap, n, alleles);
        for (int i = 0; i < n; ++i) {
            int allele = alleles[i];
            if (allele != 0) {
                if (i2hash != NULL) i2hash[i] = (int32_t)((uint32_t)i2hash[i] + (uint32_t)allele_hash[allele]);
                if (al != NULL) {
                    int_list_add(&al->alt_alleles[i], m - al->start);
                    int_list_add(&al->alt_alleles[i], allele);
                }
            }
        }
    }
    free(allele_hash);
    free(alleles);
}

static void alt_lists_init(alt_lists *al, const ref_hap_hash *rh, int start, int end) {
    al->alt_alleles = util_malloc((size_t)(rh->n > 0 ? rh->n : 1) * sizeof *al->alt_alleles);
    for (int i = 0; i < rh->n; ++i) al->alt_alleles[i] = (int_list){0};
    al->start = start;
    al->end = end;
}

static void alt_lists_free(alt_lists *al, int n) {
    for (int i = 0; i < n; ++i) free(al->alt_alleles[i].v);
    free(al->alt_alleles);
}

/* The hash of cluster, whose reference markers are [start, end); al, if not
 * NULL, receives the lists for all of them in the same pass. */
static void ref_hap_hash_init(ref_hap_hash *rh, alt_lists *al, const imp_data *id, const state_probs *sp, int cluster,
        int start, int end) {
    /* Flag each state's haplotype, then number the flagged haplotypes in
     * increasing order. */
    int *hap2i = util_malloc((size_t)(id->n_ref_haps > 0 ? id->n_ref_haps : 1) * sizeof *hap2i);
    memset(hap2i, 0, (size_t)id->n_ref_haps * sizeof *hap2i);
    for (int h = 0; h < id->n_targ_haps; ++h) {
        const cluster_probs *cp = state_probs_at(sp, h, cluster);
        for (int k = 0; k < cp->n; ++k) hap2i[cp->s[k].ref_hap] = 1;
    }
    int n = 0;
    for (int hap = 0; hap < id->n_ref_haps; ++hap) n += hap2i[hap];
    int *i2hap = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *i2hap);
    for (int hap = 0, i = 0; hap < id->n_ref_haps; ++hap) {
        if (hap2i[hap]) {
            hap2i[hap] = i;
            i2hap[i++] = hap;
        }
    }
    rh->n = n;
    rh->i2hap = i2hap;
    rh->hap2i = hap2i;
    rh->i2hash = util_malloc((size_t)(n > 0 ? n : 1) * sizeof *rh->i2hash);
    for (int i = 0; i < n; ++i) rh->i2hash[i] = 0;
    if (al != NULL) alt_lists_init(al, rh, start, end);
    scan_alleles(rh, id, start, end, rh->i2hash, al);
}

static void ref_hap_hash_free(ref_hap_hash *rh) {
    free(rh->i2hash);
    free(rh->i2hap);
    free(rh->hap2i);
}

static void set_alleles(const alt_lists *al, int index, int *alleles) {
    memset(alleles, 0, (size_t)(al->end - al->start) * sizeof *alleles);
    const int_list *il = &al->alt_alleles[index];
    for (int j = 0; j < il->n; j += 2) alleles[il->v[j]] = il->v[j + 1];
}

/* ImputedVcfWriter.setAlProbs: the haplotype's state probabilities, summed
 * over states whose haplotypes hash alike, spread over their alleles at the
 * reference markers of al. Markers before clust_end take the cluster's
 * probabilities; later ones weight them against the next cluster's. */
static void set_al_probs(const imp_data *id, const cluster_probs *cp, const ref_hap_hash *rh, const alt_lists *al,
        int clust_end, float **al_probs, int *alleles, int_list *indices, int_list *hashes, float *seq_probs,
        float *seq_probs_p1) {
    int ref_start = al->start, ref_end = al->end;
    int mid = clust_end < ref_start ? ref_start : clust_end;
    if (mid > ref_end) mid = ref_end;
    indices->n = hashes->n = 0;
    for (int j = 0; j < cp->n; ++j) {
        int index = rh->hap2i[cp->s[j].ref_hap];
        int hash = rh->i2hash[index];
        int i = 0;
        while (i < hashes->n && hashes->v[i] != hash) ++i;
        if (i == hashes->n) {
            int_list_add(indices, index);
            int_list_add(hashes, hash);
            seq_probs[i] = cp->s[j].prob;
            seq_probs_p1[i] = cp->s[j].prob_p1;
        } else {
            seq_probs[i] += cp->s[j].prob;
            seq_probs_p1[i] += cp->s[j].prob_p1;
        }
    }
    int n_seq = hashes->n;
    if (n_seq == 1) {
        set_alleles(al, indices->v[0], alleles);
        for (int m = ref_start; m < ref_end; ++m) al_probs[m - ref_start][alleles[m - ref_start]] = 1.0f;
        return;
    }
    for (int j = 0; j < n_seq; ++j) {
        set_alleles(al, indices->v[j], alleles);
        float prob = seq_probs[j];
        float prob_p1 = seq_probs_p1[j];
        for (int m = ref_start; m < mid; ++m) al_probs[m - ref_start][alleles[m - ref_start]] += prob;
        for (int m = mid; m < ref_end; ++m) {
            double wt = id->weight[m];
            al_probs[m - ref_start][alleles[m - ref_start]] += wt * prob + (1 - wt) * prob_p1;
        }
    }
}

/* Cluster c's reference markers [*ref_start, *ref_end) within the window's
 * [start, end); those from *clust_end on also use the next cluster's
 * probabilities. */
static void cluster_bounds(const imp_data *id, int c, int start, int end, int *ref_start, int *clust_end, int *ref_end) {
    *ref_start = c == 0 ? start : (start > id->ref_cluster_start[c] ? start : id->ref_cluster_start[c]);
    if (c < id->n_clusters - 1) {
        int tmp = start > id->ref_cluster_end[c] ? start : id->ref_cluster_end[c];
        *clust_end = tmp < end ? tmp : end;
        *ref_end = id->ref_cluster_start[c + 1] < end ? id->ref_cluster_start[c + 1] : end;
    } else {
        *clust_end = end;
        *ref_end = end;
    }
}

/* The most reference markers in one work item. A cluster with more is split
 * into near-equal pieces, so that one long cluster, where the target markers
 * are sparse, does not leave every other thread waiting for it. The output
 * does not depend on it: tests/check-piece-size.sh builds it at 1. */
#ifndef PIECE_RECORDS
#define PIECE_RECORDS 500
#endif

/* One work item: reference markers [start, end) of cluster. rh is the
 * cluster's hash when the cluster has more than one piece, and NULL when the
 * piece is the whole cluster. */
typedef struct {
    int cluster;
    int start, end;
    int clust_end;
    const ref_hap_hash *rh;
} piece;

/* The piece's records, into recs; returns how many. */
static int build_piece(out_rec *recs, const window_writer *ww, const imp_data *id, const state_probs *sp,
        const piece *pc) {
    const window *w = id->w;
    int cluster = pc->cluster, ref_start = pc->start, clust_end = pc->clust_end;
    ref_hap_hash own;
    alt_lists al;
    const ref_hap_hash *rh = pc->rh;
    if (rh == NULL) {
        ref_hap_hash_init(&own, &al, id, sp, cluster, ref_start, pc->end);
        rh = &own;
    } else {
        alt_lists_init(&al, rh, ref_start, pc->end);
        scan_alleles(rh, id, ref_start, pc->end, NULL, &al);
    }
    int n = pc->end - ref_start;
    float **a1 = util_malloc((size_t)n * sizeof *a1);
    float **a2 = util_malloc((size_t)n * sizeof *a2);
    bool *is_imputed = util_malloc((size_t)n * sizeof *is_imputed);
    /* One block for all the piece's records: a small block per record made
     * set_al_probs 2 to 3 times slower while other threads built pieces. */
    size_t n_probs = 0;
    for (int mm = 0; mm < n; ++mm) n_probs += (size_t)marker_n_alleles(&w->ref[ref_start + mm]->marker);
    float *probs1 = util_malloc(n_probs * sizeof *probs1);
    float *probs2 = util_malloc(n_probs * sizeof *probs2);
    for (int mm = 0, off = 0; mm < n; ++mm) {
        const marker *mk = &w->ref[ref_start + mm]->marker;
        int n_al = marker_n_alleles(mk);
        is_imputed[mm] = w->indices.marker_to_targ_marker[ref_start + mm] == -1;
        window_writer_rec_begin(ww, &recs[mm], mk, is_imputed[mm] ? OUT_IMPUTED : OUT_GENOTYPED);
        a1[mm] = probs1 + off;
        a2[mm] = probs2 + off;
        off += n_al;
        for (int a = 0; a < n_al; ++a) a1[mm][a] = a2[mm][a] = 0.0f;
    }
    int *alleles = util_malloc((size_t)n * sizeof *alleles);
    int_list indices = {0};
    int_list hashes = {0};
    int max_states = 0;
    for (int h = 0; h < id->n_targ_haps; ++h) {
        int n_kept = state_probs_at(sp, h, cluster)->n;
        if (n_kept > max_states) max_states = n_kept;
    }
    float *seq_probs = util_malloc((size_t)(max_states > 0 ? max_states : 1) * sizeof *seq_probs);
    float *seq_probs_p1 = util_malloc((size_t)(max_states > 0 ? max_states : 1) * sizeof *seq_probs_p1);
    for (int h = 0; h < id->n_targ_haps; h += 2) {
        set_al_probs(id, state_probs_at(sp, h, cluster), rh, &al, clust_end, a1, alleles, &indices, &hashes,
                seq_probs, seq_probs_p1);
        set_al_probs(id, state_probs_at(sp, h + 1, cluster), rh, &al, clust_end, a2, alleles, &indices, &hashes,
                seq_probs, seq_probs_p1);
        for (int mm = 0; mm < n; ++mm) {
            int n_al = marker_n_alleles(&w->ref[ref_start + mm]->marker);
            if (!is_imputed[mm]) {
                /* ImputedVcfWriter.setToObsAlleles: the phased alleles, using
                 * haplotype h + 1 even for a haploid sample. */
                for (int a = 0; a < n_al; ++a) a1[mm][a] = a2[mm][a] = 0.0f;
                int targ_m = w->indices.marker_to_targ_marker[ref_start + mm];
                a1[mm][id->targ_allele(id->targ_ctx, targ_m, h)] = 1.0f;
                a2[mm][id->targ_allele(id->targ_ctx, targ_m, h + 1)] = 1.0f;
            }
            window_writer_rec_probs(&recs[mm], a1[mm], a2[mm]);
            for (int a = 0; a < n_al; ++a) a1[mm][a] = a2[mm][a] = 0.0f;
        }
    }
    free(is_imputed);
    free(a1);
    free(a2);
    free(probs1);
    free(probs2);
    free(alleles);
    free(indices.v);
    free(hashes.v);
    free(seq_probs);
    free(seq_probs_p1);
    alt_lists_free(&al, rh->n);
    if (rh == &own) ref_hap_hash_free(&own);
    return n;
}

/* A cluster split into more than one piece, and its hash, which its pieces
 * share. */
typedef struct {
    int cluster;
    int start, end;
    ref_hap_hash rh;
} split_cluster;

typedef struct {
    const window_writer *ww;
    const imp_data *id;
    const state_probs *sp;
    const piece *pieces;
    split_cluster *splits;
} build_ctx;

static void hash_task(void *ctx, int split) {
    const build_ctx *c = ctx;
    split_cluster *sc = &c->splits[split];
    ref_hap_hash_init(&sc->rh, NULL, c->id, c->sp, sc->cluster, sc->start, sc->end);
}

static int build_item(const void *ctx, int item, out_rec *recs) {
    const build_ctx *c = ctx;
    return build_piece(recs, c->ww, c->id, c->sp, &c->pieces[item]);
}

static void trace_item(const void *ctx, int item, int j, const out_rec *r) {
    const build_ctx *c = ctx;
    window_writer_rec_trace(r, c->pieces[item].cluster, c->pieces[item].start + j);
}

/* Splits the window's clusters into pieces of at most PIECE_RECORDS
 * reference markers, in cluster order, and lists the clusters that split.
 * The pieces of a split cluster point to its hash, which is built later. */
static piece *make_pieces(const imp_data *id, int start, int end, int *n_pieces, split_cluster **splits, int *n_splits) {
    /* each split cluster has more than PIECE_RECORDS markers, and adds
     * fewer than n / PIECE_RECORDS pieces beyond its first */
    int max_splits = (end - start) / PIECE_RECORDS + 1;
    piece *pieces = util_malloc((size_t)(id->n_clusters + max_splits) * sizeof *pieces);
    split_cluster *sc = util_malloc((size_t)max_splits * sizeof *sc);
    int np = 0, ns = 0;
    for (int c = 0; c < id->n_clusters; ++c) {
        int ref_start, clust_end, ref_end;
        cluster_bounds(id, c, start, end, &ref_start, &clust_end, &ref_end);
        if (ref_start >= ref_end) continue;
        int n = ref_end - ref_start;
        int k = (n + PIECE_RECORDS - 1) / PIECE_RECORDS;
        const ref_hap_hash *rh = NULL;
        if (k > 1) {
            sc[ns] = (split_cluster){c, ref_start, ref_end, {0}};
            rh = &sc[ns++].rh;
        }
        for (int j = 0; j < k; ++j) {
            pieces[np++] = (piece){c, ref_start + (int)((int64_t)n * j / k), ref_start + (int)((int64_t)n * (j + 1) / k),
                    clust_end, rh};
        }
    }
    *n_pieces = np;
    *n_splits = ns;
    *splits = sc;
    return pieces;
}

void imputed_writer_print(window_writer *ww, const imp_data *id, const state_probs *sp, const par *p, int start, int end) {
    int n_pieces, n_splits;
    split_cluster *splits;
    piece *pieces = make_pieces(id, start, end, &n_pieces, &splits, &n_splits);
    build_ctx ctx = {ww, id, sp, pieces, splits};
    parallel_for(p->nthreads, n_splits, &ctx, 0, hash_task);
    window_writer_print_records(ww, n_pieces, PIECE_RECORDS, build_item, trace_on() ? trace_item : NULL, &ctx);
    for (int s = 0; s < n_splits; ++s) ref_hap_hash_free(&splits[s].rh);
    free(splits);
    free(pieces);
}
