/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) main/WindowWriter.java,
 * vcf/VcfWriter.java, vcf/VcfRecBuilder.java and imp/ImputedRecBuilder.java;
 * modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "main/window_writer.h"

#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "blbutil/parallel.h"
#include "blbutil/trace.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"

/* ImputedRecBuilder.DS_VALS ("#.##" of j/100) and R2_VALS ("0.00" where
 * DS_VALS is not 4 characters long). */
static char ds_vals[201][8];
static char r2_vals[101][8];

static void init_tables(void) {
    for (int j = 0; j < 201; ++j) jnum_format_hash2(ds_vals[j], sizeof ds_vals[j], j / 100.0);
    for (int i = 0; i < 101; ++i) {
        if (strlen(ds_vals[i]) != 4) jnum_format_fixed(r2_vals[i], sizeof r2_vals[i], i / 100.0, 2);
        else memcpy(r2_vals[i], ds_vals[i], sizeof r2_vals[i]);
    }
}

/* (int) Math.rint(100 * x) as a table index. */
static const char *ds(float x) {
    return ds_vals[jnum_d2i(rint((double)(100 * x)))];
}

typedef struct {
    int a1, a2;
    const char *ds;
    const float *p1, *p2;
} out_call;

static void scale(float *fa, int n) {
    float sum = 0.0f;
    for (int j = 0; j < n; ++j) sum += fa[j];
    for (int j = 0; j < n; ++j) fa[j] /= sum;
}

static int max_index(const float *fa, int n) {
    int max = 0;
    for (int j = 1; j < n; ++j) {
        if (fa[j] > fa[max]) max = j;
    }
    return max;
}

static uint64_t fold_probs(uint64_t h, const float *fa, int n) {
    for (int j = 0; j < n; ++j) h = trace_fold(h, jnum_float_bits(fa[j]));
    return h;
}

static void put_sep(kstring_t *s, int first) {
    kputc(first ? ':' : ',', s);
}

/* ImputedRecBuilder.addSampleData for a diploid sample: returns what its
 * fields printed, with a1 and a2 scaled to sum to 1 unless they were
 * certain. */
static out_call add_diploid(out_rec *rb, float *a1, float *a2, bool ap, bool gp, const kstring_t *hom_ref) {
    int n = rb->n_alleles;
    kstring_t *s = &rb->fields;
    rb->hap_cnt += 2;
    if (trace_on()) rb->digest = fold_probs(fold_probs(rb->digest, a1, n), a2, n);
    if (a1[0] == 1.0f && a2[0] == 1.0f && n < 5) {
        kputs(hom_ref[n].s, s);
        return (out_call){0, 0, n > 1 ? ds_vals[0] : NULL, NULL, NULL};
    }
    scale(a1, n);
    scale(a2, n);
    int allele1 = max_index(a1, n), allele2 = max_index(a2, n);
    kputc('\t', s);
    kputw(allele1, s);
    kputc('|', s);
    kputw(allele2, s);
    const char *ds1 = NULL;
    for (int a = 1; a < n; ++a) {
        float dose = a1[a] + a2[a];
        float dose2 = a1[a] * a1[a] + a2[a] * a2[a];
        rb->sum_al_probs[a] += dose;
        rb->sum_al_probs2[a] += dose2;
        put_sep(s, a == 1);
        const char *dose_str = ds(dose);
        if (a == 1) ds1 = dose_str;
        kputs(dose_str, s);
    }
    if (ap) {
        for (int a = 1; a < n; ++a) {
            put_sep(s, a == 1);
            kputs(ds(a1[a]), s);
        }
        for (int a = 1; a < n; ++a) {
            put_sep(s, a == 1);
            kputs(ds(a2[a]), s);
        }
    }
    if (gp) {
        for (int i2 = 0; i2 < n; ++i2) {
            for (int i1 = 0; i1 <= i2; ++i1) {
                float prob = a1[i1] * a2[i2];
                if (i1 != i2) prob += a1[i2] * a2[i1];
                put_sep(s, i2 == 0);
                kputs(ds(prob), s);
            }
        }
    }
    return (out_call){allele1, allele2, ds1, a1, a2};
}

static out_call add_haploid(out_rec *rb, float *a1, bool ap) {
    int n = rb->n_alleles;
    kstring_t *s = &rb->fields;
    ++rb->hap_cnt;
    if (trace_on()) rb->digest = fold_probs(rb->digest, a1, n);
    scale(a1, n);
    int allele1 = max_index(a1, n);
    kputc('\t', s);
    kputw(allele1, s);
    const char *ds1 = NULL;
    for (int a = 1; a < n; ++a) {
        float dose = a1[a];
        float dose2 = a1[a] * a1[a];
        rb->sum_al_probs[a] += dose;
        rb->sum_al_probs2[a] += dose2;
        put_sep(s, a == 1);
        const char *dose_str = ds(dose);
        if (a == 1) ds1 = dose_str;
        kputs(dose_str, s);
    }
    if (ap) {
        for (int a = 1; a < n; ++a) {
            put_sep(s, a == 1);
            kputs(ds(a1[a]), s);
        }
    }
    return (out_call){allele1, 0, ds1, a1, NULL};
}

static float r2(const out_rec *rb, int allele, int n_input_targ_haps) {
    float sum = rb->sum_al_probs[allele];
    if (sum == 0.0f) return 0.0f;
    float sum2 = rb->sum_al_probs2[allele];
    float mean_term = sum * sum / n_input_targ_haps;
    float num = sum2 - mean_term;
    float den = sum - mean_term;
    return num <= 0 ? 0.0f : num / den;
}

/* Trace seam T5d: a record's cluster and reference marker, the digest of
 * its samples' raw allele probabilities, and per ALT allele the dose sums
 * and DR2 as raw bits. */
void window_writer_rec_trace(const out_rec *rb, int cluster, int m) {
    int n_input_targ_haps = rb->writer->n_haps;
    kstring_t s = {0, 0, NULL};
    ksprintf(&s, "rec\t%d\t%d\t%" PRIx64, cluster, m, rb->digest);
    for (int a = 1; a < rb->n_alleles; ++a) {
        ksprintf(&s, "\t%" PRIx32 ":%" PRIx32 ":%" PRIx32, jnum_float_bits(rb->sum_al_probs[a]),
                jnum_float_bits(rb->sum_al_probs2[a]), jnum_float_bits(r2(rb, a, n_input_targ_haps)));
    }
    trace_line("T5d", "%s", s.s);
    free(s.s);
}

static void put_info(const out_rec *rb, bool is_imputed, int n_input_targ_haps, kstring_t *l) {
    if (rb->n_alleles == 1) {
        if (is_imputed) kputs("IMP", l);
        return;
    }
    char buf[32];
    for (int a = 1; a < rb->n_alleles; ++a) {
        kputs(a == 1 ? "DR2=" : ",", l);
        kputs(r2_vals[jnum_d2i(rint((double)(100 * r2(rb, a, n_input_targ_haps))))], l);
    }
    for (int a = 1; a < rb->n_alleles; ++a) {
        kputs(a == 1 ? ";AF=" : ",", l);
        jnum_format_fixed(buf, sizeof buf, (double)(rb->sum_al_probs[a] / n_input_targ_haps), 4);
        kputs(buf, l);
    }
    /* ImputedRecBuilder.extractEnd */
    span info = marker_info(rb->mk);
    for (int j = 0; j + 4 <= info.n; ++j) {
        if (memcmp(info.s + j, "END=", 4) == 0) {
            int k = j + 4;
            while (k < info.n && info.s[k] != ';') ++k;
            kputc(';', l);
            kputsn(info.s + j, (size_t)(k - j), l);
            break;
        }
    }
    if (is_imputed) kputs(";IMP", l);
}

/* ImputedRecBuilder.defaultHomRefFields and homRefFields */
static void hom_ref_fields(kstring_t *fields, bool ap, bool gp) {
    kputs("\t0|0", &fields[1]);
    kputs("\t0|0:0", &fields[2]);
    for (int j = 3; j < 5; ++j) {
        kputsn(fields[j - 1].s, fields[j - 1].l, &fields[j]);
        kputs(",0", &fields[j]);
    }
    for (int n_al = 1; n_al < 5; ++n_al) {
        kstring_t *s = &fields[n_al];
        if (ap) {
            for (int rep = 0; rep < 2; ++rep) {
                for (int a = 1; a < n_al; ++a) {
                    put_sep(s, a == 1);
                    kputs(ds_vals[0], s);
                }
            }
        }
        if (gp) {
            kputc(':', s);
            kputs(ds_vals[100], s);
            for (int i2 = 1; i2 < n_al; ++i2) {
                for (int i1 = 0; i1 <= i2; ++i1) {
                    kputc(',', s);
                    kputs(ds_vals[0], s);
                }
            }
        }
    }
}

struct out_worker {
    bgen_scratch *scratch;   /* NULL without bgen= */
};

/* Writes ww->line followed by a newline and empties it. */
static void write_line(window_writer *ww) {
    kputc('\n', &ww->line);
    if (bgzf_write(ww->out, ww->line.s, ww->line.l) < 0) util_exit("Error writing %s", ww->path);
    ww->n_written += ww->line.l;
    ww->line.l = 0;
}

void window_writer_open(window_writer *ww, const par *p, const run_outputs *out, const samples *s) {
    static pthread_once_t tables_once = PTHREAD_ONCE_INIT;
    pthread_once(&tables_once, init_tables);
    ww->samples = s;
    ww->ap = p->ap;
    ww->gp = p->gp;
    ww->nthreads = p->nthreads;
    ww->n_haps = 0;
    for (int j = 0; j < s->n; ++j) ww->n_haps += s->is_diploid[j] ? 2 : 1;
    for (int j = 0; j < 5; ++j) ww->hom_ref[j] = (kstring_t){0};
    hom_ref_fields(ww->hom_ref, p->ap, p->gp);
    ww->line = (kstring_t){0, 0, NULL};
    ww->bgen = p->bgen != BGEN_NONE ? bgen_writer_open(p, out, s) : NULL;
    ww->path = run_outputs_path(out, RUN_OUTPUT_VCF);
    ww->n_written = 0;
    ww->out = bgzf_open(ww->path, "w");
    if (ww->out == NULL) util_exit("Error opening %s", ww->path);
    if (p->nthreads > 1 && bgzf_mt(ww->out, p->nthreads, 256) != 0) util_exit("Error opening %s", ww->path);

    char date[16];
    time_t now = time(NULL);
    strftime(date, sizeof date, "%Y%m%d", localtime(&now));
    kstring_t *l = &ww->line;
    kputs("##fileformat=VCFv4.2", l);
    write_line(ww);
    ksprintf(l, "##filedate=%s", date);
    write_line(ww);
    kputs("##source=\"beagle.27Feb25.75f.jar\"", l);
    write_line(ww);
    kputs("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">", l);
    write_line(ww);
    kputs("##INFO=<ID=DR2,Number=A,Type=Float,Description=\"Dosage R-Squared: estimated squared correlation between "
            "estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">", l);
    write_line(ww);
    kputs("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">", l);
    write_line(ww);
    kputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", l);
    write_line(ww);
    kputs("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">", l);
    write_line(ww);
    if (p->ap) {
        kputs("##FORMAT=<ID=AP1,Number=A,Type=Float,Description=\"estimated ALT dose on first haplotype\">", l);
        write_line(ww);
        kputs("##FORMAT=<ID=AP2,Number=A,Type=Float,Description=\"estimated ALT dose on second haplotype\">", l);
        write_line(ww);
    }
    if (p->gp) {
        kputs("##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">", l);
        write_line(ww);
    }
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", l);
    for (int j = 0; j < s->n; ++j) {
        kputc('\t', l);
        kputs(s->ids[j], l);
    }
    write_line(ww);
    ww->index = p->tbi ? vcf_index_new(ww->path, run_outputs_path(out, RUN_OUTPUT_TBI), ww->n_written) : NULL;
}

void window_writer_begin_window(window_writer *ww, const window *w) {
    if (ww->bgen != NULL) bgen_writer_check_chrom(ww->bgen, marker_chrom(&w->targ[0]->marker));
}

static void put_span(span sp, kstring_t *l) {
    kputsn(sp.s, (size_t)sp.n, l);
}

void window_writer_rec_begin(const window_writer *ww, out_rec *r, const marker *mk, out_rec_kind kind) {
    r->writer = ww;
    r->mk = mk;
    r->kind = kind;
    r->n_alleles = marker_n_alleles(mk);
    r->n_samples = r->hap_cnt = 0;
    r->digest = TRACE_FNV_BASIS;
    r->info.l = 0;
    r->fields.l = 0;
    /* allocated even when empty, as the outputs memcpy from info.s */
    if (ks_resize(&r->info, 32) != 0) util_oom();
    if (ww->bgen != NULL) r->bgen = bgen_rec_begin(ww->bgen, r->bgen, mk);
    if (kind == OUT_PHASED) {
        put_span(marker_info(mk), &r->info);
        kputs("GT", &r->fields);
    } else {
        int n = r->n_alleles;
        r->sum_al_probs = util_realloc(r->sum_al_probs, (size_t)n * sizeof *r->sum_al_probs);
        r->sum_al_probs2 = util_realloc(r->sum_al_probs2, (size_t)n * sizeof *r->sum_al_probs2);
        for (int a = 0; a < n; ++a) r->sum_al_probs[a] = r->sum_al_probs2[a] = 0.0f;
        kputs("GT:DS", &r->fields);
        if (ww->ap) kputs(":AP1:AP2", &r->fields);
        if (ww->gp) kputs(":GP", &r->fields);
    }
}

void window_writer_rec_free(out_rec *r) {
    free(r->sum_al_probs);
    free(r->sum_al_probs2);
    free(r->info.s);
    free(r->fields.s);
    bgen_rec_free(r->bgen);
    *r = (out_rec){0};
}

/* Capture before the producer reuses its probability buffers. PLINK2 keeps
 * the printed DS string; phased BGEN quantises the scaled probabilities. */
static void capture_call(out_rec *r, out_call c) {
    if (r->bgen != NULL) bgen_rec_set_call(r->bgen, r->n_samples, c.a1, c.a2, c.ds, c.p1, c.p2);
    ++r->n_samples;
}

void window_writer_rec_gt(out_rec *r, int a1, int a2) {
    kputc('\t', &r->fields);
    kputw(a1, &r->fields);
    bool diploid = r->writer->samples->is_diploid[r->n_samples];
    r->hap_cnt += diploid ? 2 : 1;
    if (diploid) {
        kputc('|', &r->fields);
        kputw(a2, &r->fields);
    }
    capture_call(r, (out_call){a1, a2, NULL, NULL, NULL});
}

void window_writer_rec_probs(out_rec *r, float *p1, float *p2) {
    const window_writer *ww = r->writer;
    out_call c = ww->samples->is_diploid[r->n_samples]
        ? add_diploid(r, p1, p2, ww->ap, ww->gp, ww->hom_ref)
        : add_haploid(r, p1, ww->ap);
    capture_call(r, c);
}

out_worker *window_writer_worker_new(const window_writer *ww) {
    out_worker *wk = util_malloc(sizeof *wk);
    wk->scratch = ww->bgen != NULL ? bgen_scratch_new(ww->bgen) : NULL;
    return wk;
}

void window_writer_worker_free(out_worker *wk) {
    if (wk->scratch != NULL) bgen_scratch_free(wk->scratch);
    free(wk);
}

void window_writer_encode(out_worker *wk, out_rec *r) {
    if (r->n_samples != r->writer->samples->n || r->hap_cnt != r->writer->n_haps)
        util_exit("java.lang.IllegalStateException: inconsistent data");
    if (r->kind != OUT_PHASED) put_info(r, r->kind == OUT_IMPUTED, r->writer->n_haps, &r->info);
    if (r->bgen != NULL) bgen_writer_encode(r->writer->bgen, wk->scratch, (span){r->info.s, (int)r->info.l}, r->bgen);
}

/* Marker.toString, then the record's fields. */
void window_writer_put(window_writer *ww, out_rec *r) {
    kstring_t *l = &ww->line;
    const marker *mk = r->mk;
    kputs(marker_chrom(mk), l);
    kputc('\t', l);
    kputw(mk->pos, l);
    kputc('\t', l);
    put_span(marker_id(mk), l);
    kputc('\t', l);
    put_span(marker_alleles(mk), l);
    kputs("\t.\t", l);
    kputs(r->kind == OUT_PHASED ? "." : "PASS", l);
    kputc('\t', l);
    kputsn(r->info.s, r->info.l, l);
    size_t info_end = l->l;
    kputc('\t', l);
    kputsn(r->fields.s, r->fields.l, l);
    /* + 1 for the newline write_line adds */
    if (ww->index != NULL) vcf_index_add(ww->index, l->s, info_end, ww->n_written + l->l + 1);
    write_line(ww);
    if (ww->bgen != NULL) bgen_writer_put(ww->bgen, r->bgen);
}

/* WindowWriter.printPhased(Stage2Haps, ...) formats 50 records per task. At
 * most PHASED_SLOTS_PER_THREAD * nthreads tasks' records are held at once. */
#define PHASED_STEP 50
#define PHASED_SLOTS_PER_THREAD 2

typedef struct {
    window_writer *ww;
    const window *w;
    int start, end;
    phased_allele_fn allele;
    const void *ctx;
    int window;
    out_rec *recs;  /* task i's records in slot i % window, PHASED_STEP per slot */
} phased_ctx;

typedef struct {
    const phased_ctx *c;
    out_worker *wk;
} phased_worker;

/* Task item's first record, and its target markers in [*m_start, *m_end). */
static out_rec *phased_task(const phased_ctx *c, int item, int *m_start, int *m_end) {
    *m_start = c->start + item * PHASED_STEP;
    *m_end = c->end - *m_start < PHASED_STEP ? c->end : *m_start + PHASED_STEP;
    return &c->recs[(size_t)(item % c->window) * PHASED_STEP];
}

static void phased_build(void *worker, int item) {
    const phased_worker *pw = worker;
    const phased_ctx *c = pw->c;
    const window_writer *ww = c->ww;
    int m_start, m_end;
    out_rec *r = phased_task(c, item, &m_start, &m_end);
    for (int m = m_start; m < m_end; ++m, ++r) {
        window_writer_rec_begin(ww, r, &c->w->targ[m]->marker, OUT_PHASED);
        for (int s = 0; s < ww->samples->n; ++s) {
            int h1 = s << 1;
            int a1 = c->allele(c->ctx, m, h1), a2 = 0;
            if (ww->samples->is_diploid[s]) a2 = c->allele(c->ctx, m, h1 | 1);
            window_writer_rec_gt(r, a1, a2);
        }
        window_writer_encode(pw->wk, r);
    }
}

static void phased_put(void *ctx, int item) {
    const phased_ctx *c = ctx;
    int m_start, m_end;
    out_rec *r = phased_task(c, item, &m_start, &m_end);
    for (int m = m_start; m < m_end; ++m, ++r) window_writer_put(c->ww, r);
}

void window_writer_print_phased(window_writer *ww, const window *w, int start, int end, phased_allele_fn allele, const void *ctx) {
    int n_threads = ww->nthreads;
    int window = PHASED_SLOTS_PER_THREAD * n_threads;
    size_t n_recs = (size_t)window * PHASED_STEP;
    out_rec *recs = util_malloc(n_recs * sizeof *recs);
    for (size_t j = 0; j < n_recs; ++j) recs[j] = (out_rec){0};
    phased_ctx c = {ww, w, start, end, allele, ctx, window, recs};
    phased_worker *pws = util_malloc((size_t)n_threads * sizeof *pws);
    for (int t = 0; t < n_threads; ++t) pws[t] = (phased_worker){&c, window_writer_worker_new(ww)};
    int n_tasks = (end - start + PHASED_STEP - 1) / PHASED_STEP;
    parallel_ordered(n_threads, n_tasks, window, pws, sizeof *pws, phased_build, &c, phased_put);
    for (int t = 0; t < n_threads; ++t) window_writer_worker_free(pws[t].wk);
    for (size_t j = 0; j < n_recs; ++j) window_writer_rec_free(&recs[j]);
    free(pws);
    free(recs);
}

void window_writer_close(window_writer *ww) {
    if (bgzf_close(ww->out) != 0) util_exit("Error writing %s", ww->path);
    if (ww->index != NULL) vcf_index_write(ww->index, ww->n_written);
    if (ww->bgen != NULL) bgen_writer_close(ww->bgen);
    for (int j = 0; j < 5; ++j) free(ww->hom_ref[j].s);
    free(ww->line.s);
}
