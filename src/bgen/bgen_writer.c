/*
 * Copyright (C) 2005-2026 Shaun Purcell, Christopher Chang
 * Ported to C from PLINK 2.0 v2.0.0-a.7.8 plink2_import.cc (VCF dosage and
 * hardcall import), plink2_data.cc (ApplyHardCallThreshPhased),
 * plink2_pvar.cc (--extract-if-info), plink2_filter.cc (--maf) and
 * plink2_export.cc (--export bgen-1.2, .sample); modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "bgen/bgen_writer.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/kstring.h>
#include <libdeflate.h>

#include "bgen/bgen_files.h"
#include "bgen/plink2_num.h"
#include "blbutil/utilities.h"

/* Dosages in 1/16384 ALT alleles (plink2_common.h). */
enum { DOSAGE_MAX = 32768, DOSAGE_MID = 16384, DOSAGE_4TH = 8192 };
/* --hard-call-threshold 0.1 */
enum { HARD_CALL_HALFDIST = DOSAGE_4TH - DOSAGE_MID / 10 };
enum { BGEN_HEADER_LEN = 20, BGEN_FLAGS = 0x80000009 };

/* A sample after plink2's import: hardcall g (3 = missing), the stored
 * dosage, hardcall phase (pi: the GT was 1|0) and dosage phase. */
typedef struct {
    int g;
    bool has_dosage, hp, pi, has_dphase;
    int32_t d, dphase;
} plink2_geno;

struct bgen_writer {
    bgen_mode mode;
    FILE *bgen, *info_file;
    const char *bgen_path, *sample_path, *info_path;
    const samples *samples;
    bool has_min_dr2;
    double min_dr2, min_maf;
    int bits;
    uint32_t max;
    uint32_t n_variants;
    long n_multiallelic;
    uint32_t *n_missing;
    int chrom_index;
    int autosome_ct;
    kstring_t rec;
};

/* One sample's fields as the VCF record prints them: its GT alleles and its
 * DS string, NULL when not printed. For bgen=phased, hap_q holds each
 * haplotype's probabilities of alleles 0 to K-2 from bgen_quantise,
 * haplotype by haplotype; NULL means the GT alleles are certain. */
typedef struct {
    int a1, a2;
    const char *ds;
    const uint32_t *hap_q;
} bgen_call;

/* What bgen_writer_put does with an encoded record. */
typedef enum { BGEN_REC_WRITE, BGEN_REC_FILTERED, BGEN_REC_MULTIALLELIC } bgen_rec_kind;

/* The calls live from bgen_rec_begin to bgen_writer_encode, so records
 * waiting to be put hold only their output. For BGEN_REC_WRITE, bytes holds
 * the record header and compressed probability data, info the .info row,
 * and missing the n_missing samples plink2 counts as missing in the .sample
 * file. */
struct bgen_rec {
    const marker *mk;
    bgen_call *calls;  /* per sample */
    uint32_t *hap_q;   /* 2 * (n_alleles - 1) per sample for bgen=phased; NULL otherwise */
    int n_alleles;
    uint32_t max;      /* the value of probability 1 in hap_q */
    bgen_rec_kind kind;
    uint8_t *bytes;
    size_t len;
    kstring_t info;
    int *missing;
    int n_missing;
};

struct bgen_scratch {
    struct libdeflate_compressor *compressor;
    plink2_geno *geno;
    uint32_t *probs;
    size_t probs_len;
    kstring_t block, rec, info;
    uint8_t *zbuf;
    size_t zbuf_len;
};

static void put_u16(kstring_t *s, uint32_t v) {
    if (v > 0xffff) util_exit("bgen: field of %u bytes is too long for BGEN", v);
    kputc_((char)(v & 0xff), s);
    kputc_((char)(v >> 8), s);
}

static void put_u32(kstring_t *s, uint32_t v) {
    for (int k = 0; k < 4; ++k) kputc_((char)((v >> (8 * k)) & 0xff), s);
}

static void write_bytes(bgen_writer *bw, const void *p, size_t n) {
    if (fwrite(p, 1, n, bw->bgen) != n) util_exit("Error writing %s", bw->bgen_path);
}

/* GetChrCodeRaw for the autosomes, an optional "chr" and one or two digits:
 * the code of chrom if GetChrCode makes it an autosome under --chr-set
 * autosome_ct, or 0 for any other chromosome. */
static int autosome_code(const bgen_writer *bw, const char *chrom) {
    const char *c = chrom;
    if ((c[0] | 0x20) == 'c' && (c[1] | 0x20) == 'h' && (c[2] | 0x20) == 'r') c += 3;
    int code = 0;
    if (c[0] >= '0' && c[0] <= '9') {
        if (c[1] == '\0') code = c[0] - '0';
        else if (c[1] >= '0' && c[1] <= '9' && c[2] == '\0') code = 10 * (c[0] - '0') + (c[1] - '0');
    }
    return code <= bw->autosome_ct ? code : 0;
}

void bgen_writer_check_chrom(bgen_writer *bw, const char *chrom) {
    if (bw->mode == BGEN_PLINK2 && autosome_code(bw, chrom) == 0) {
        util_exit("bgen=plink2 supports autosomes 1-%d only, not chromosome %s: plink2 needs sex information for "
                "chrX, chrY and MT, and --allow-extra-chr for other contigs (bgen-chr-set=N sets the autosome count)",
                bw->autosome_ct, chrom);
    }
}

bgen_writer *bgen_writer_open(const par *p, const run_outputs *out, const samples *s) {
    bgen_writer *bw = util_malloc(sizeof *bw);
    *bw = (bgen_writer){0};
    bw->mode = p->bgen;
    bw->bgen_path = run_outputs_path(out, RUN_OUTPUT_BGEN);
    bw->sample_path = run_outputs_path(out, RUN_OUTPUT_SAMPLE);
    bw->info_path = run_outputs_path(out, RUN_OUTPUT_INFO);
    bw->samples = s;
    bw->has_min_dr2 = p->has_bgen_min_dr2;
    bw->min_dr2 = p->bgen_min_dr2;
    /* EnforceFreqConstraints: defend against floating point error */
    bw->min_maf = p->bgen_min_maf * (1.0 - 0x1p-44);
    bw->autosome_ct = p->bgen_chr_set;
    bw->bits = p->bgen_bits;
    bw->max = (1u << p->bgen_bits) - 1;
    bw->n_missing = util_malloc((size_t)(s->n > 0 ? s->n : 1) * sizeof *bw->n_missing);
    for (int j = 0; j < s->n; ++j) bw->n_missing[j] = 0;
    bw->chrom_index = -1;
    bw->bgen = bgen_files_open(bw->bgen_path);
    if (bw->bgen == NULL) util_exit("Error opening %s", bw->bgen_path);
    bw->info_file = bgen_files_open(bw->info_path);
    if (bw->info_file == NULL) util_exit("Error opening %s", bw->info_path);
    fputs("CHROM\tPOS\tID\tREF\tALT\tDR2\tAF\tIMP\n", bw->info_file);

    uint32_t sample_block_len = 8;
    for (int j = 0; j < s->n; ++j) sample_block_len += 2 + (uint32_t)strlen(s->ids[j]);
    kstring_t *h = &bw->rec;
    put_u32(h, BGEN_HEADER_LEN + sample_block_len);
    put_u32(h, BGEN_HEADER_LEN);
    put_u32(h, 0);  /* variant count, written at close */
    put_u32(h, (uint32_t)s->n);
    kputsn("bgen", 4, h);
    put_u32(h, BGEN_FLAGS);
    put_u32(h, sample_block_len);
    put_u32(h, (uint32_t)s->n);
    for (int j = 0; j < s->n; ++j) {
        size_t len = strlen(s->ids[j]);
        put_u16(h, (uint32_t)len);
        kputsn(s->ids[j], len, h);
    }
    write_bytes(bw, h->s, h->l);
    h->l = 0;
    return bw;
}

bgen_scratch *bgen_scratch_new(const bgen_writer *bw) {
    bgen_scratch *sc = util_malloc(sizeof *sc);
    *sc = (bgen_scratch){0};
    sc->compressor = libdeflate_alloc_compressor(6);
    if (sc->compressor == NULL) util_exit("libdeflate_alloc_compressor failed");
    sc->geno = util_malloc((size_t)(bw->samples->n > 0 ? bw->samples->n : 1) * sizeof *sc->geno);
    return sc;
}

void bgen_scratch_free(bgen_scratch *sc) {
    libdeflate_free_compressor(sc->compressor);
    free(sc->geno);
    free(sc->probs);
    free(sc->block.s);
    free(sc->rec.s);
    free(sc->info.s);
    free(sc->zbuf);
    free(sc);
}

static void release_calls(bgen_rec *rec) {
    free(rec->calls);
    free(rec->hap_q);
    rec->calls = NULL;
    rec->hap_q = NULL;
}

bgen_rec *bgen_rec_begin(const bgen_writer *bw, bgen_rec *rec, const marker *mk) {
    if (rec == NULL) {
        rec = util_malloc(sizeof *rec);
        *rec = (bgen_rec){0};
    }
    release_calls(rec);
    size_t n = (size_t)bw->samples->n;
    int n_alleles = marker_n_alleles(mk);
    rec->mk = mk;
    rec->calls = util_malloc(n * sizeof *rec->calls);
    rec->hap_q = bw->mode == BGEN_PHASED && n_alleles > 1 ? util_malloc(n * 2 * (size_t)(n_alleles - 1) * sizeof *rec->hap_q) : NULL;
    rec->n_alleles = n_alleles;
    rec->max = bw->max;
    return rec;
}

void bgen_rec_set_call(bgen_rec *rec, int s, int a1, int a2, const char *ds, const float *p1, const float *p2) {
    int n = rec->n_alleles;
    uint32_t *q = rec->hap_q != NULL && p1 != NULL ? rec->hap_q + (size_t)s * 2 * (size_t)(n - 1) : NULL;
    rec->calls[s] = (bgen_call){a1, a2, ds, q};
    if (q == NULL) return;
    bgen_quantise(p1, n, rec->max, q);
    if (p2 != NULL) bgen_quantise(p2, n, rec->max, q + (n - 1));
}

void bgen_rec_free(bgen_rec *rec) {
    if (rec == NULL) return;
    release_calls(rec);
    free(rec->bytes);
    free(rec->info.s);
    free(rec->missing);
    free(rec);
}

/* BiallelicDosageHalfdist */
static int32_t halfdist(int32_t d) {
    int32_t h = (d & (DOSAGE_MID - 1)) - DOSAGE_4TH;
    return h < 0 ? -h : h;
}

/* A phased biallelic sample through VCF import with dosage=DS and then
 * ApplyHardCallThreshPhased. */
static plink2_geno import_call(const bgen_call *c, bool diploid) {
    plink2_geno x = {0};
    if (diploid) {
        x.g = c->a1 + c->a2;
        if (x.g == 1) {
            x.hp = true;
            x.pi = c->a1 == 1;
        }
    } else {
        x.g = 2 * c->a1;
    }
    if (c->ds == NULL) return x;
    double v;
    if (plink2_scanadv_double(c->ds, &v) == NULL || v < 0.0) util_exit("bgen=plink2: invalid dosage %s", c->ds);
    if (!diploid) v *= 2;
    if (v > 2.0) util_exit("bgen=plink2: invalid dosage %s", c->ds);
    int32_t d = (int32_t)(v * DOSAGE_MID + 0.5);
    int32_t hd = halfdist(d);
    if (hd == DOSAGE_4TH) {
        /* an integer dosage is not stored, but overrides the hardcall */
        x.g = (d + DOSAGE_4TH) / DOSAGE_MID;
        if (x.hp && x.g != 1) x.hp = false;
        return x;
    }
    x.has_dosage = true;
    x.d = d;
    int new_g = hd < HARD_CALL_HALFDIST ? 3 : (d + DOSAGE_4TH) / DOSAGE_MID;
    if (new_g != x.g) {
        if (x.hp) {
            x.hp = false;
            x.has_dphase = true;
            x.dphase = d <= DOSAGE_MID ? d : DOSAGE_MAX - d;
            if (!x.pi) x.dphase = -x.dphase;
        }
        x.g = new_g;
    }
    return x;
}

/* ComputeAlleleFreqs and GetTypedFreq on an autosome: the minor allele
 * frequency from dosages, or hardcalls where no dosage is stored. */
static double maf(const plink2_geno *geno, int n) {
    uint64_t alt_sum = 0, ref_sum = 0;
    for (int j = 0; j < n; ++j) {
        int32_t alt;
        if (geno[j].has_dosage) alt = geno[j].d;
        else if (geno[j].g != 3) alt = geno[j].g * DOSAGE_MID;
        else continue;
        alt_sum += (uint64_t)alt;
        ref_sum += (uint64_t)(DOSAGE_MAX - alt);
    }
    uint64_t total = alt_sum + ref_sum;
    double ref_freq = total == 0 ? 0.5 : (double)ref_sum * (1.0 / (double)total);
    double alt_freq = 1.0 - ref_freq;
    return alt_freq < ref_freq ? alt_freq : ref_freq;
}

static uint32_t scaled(uint32_t num, uint32_t max) {
    return (num * max + DOSAGE_4TH) / DOSAGE_MID;
}

void bgen_pack_bits(const uint32_t *v, size_t n, int bits, uint8_t *out) {
    uint32_t acc = 0;
    int acc_bits = 0;
    for (size_t k = 0; k < n; ++k) {
        acc |= v[k] << acc_bits;
        acc_bits += bits;
        for (; acc_bits >= 8; acc_bits -= 8) {
            *out++ = (uint8_t)acc;
            acc >>= 8;
        }
    }
    if (acc_bits > 0) *out = (uint8_t)acc;
}

/* Appends v[0, n) packed at bits each. */
static void put_bits(kstring_t *b, const uint32_t *v, size_t n, int bits) {
    size_t len = (n * (size_t)bits + 7) / 8;
    if (ks_resize(b, b->l + len) != 0) util_oom();
    bgen_pack_bits(v, n, bits, (uint8_t *)b->s + b->l);
    b->l += len;
}

/* sc->probs with room for n values. */
static uint32_t *probs_buf(bgen_scratch *sc, size_t n) {
    if (n > sc->probs_len) {
        sc->probs = util_realloc(sc->probs, n * sizeof *sc->probs);
        sc->probs_len = n;
    }
    return sc->probs;
}

/* The uncompressed probability data of a diploid autosome record
 * (plink2_export.cc, the !is_haploid path). */
static void probability_block(const bgen_writer *bw, bgen_scratch *sc, const plink2_geno *geno, int n) {
    kstring_t *b = &sc->block;
    uint32_t max = bw->max;
    /* ConstructBgen13LookupTables: half_val_roundeven */
    uint32_t half = ((max + 1) / 2) & ~1u;
    bool phased = false;
    bool dphase_on = false;
    for (int j = 0; j < n; ++j) {
        const plink2_geno *x = &geno[j];
        if (x->hp) phased = true;
        if (x->has_dphase) {
            uint32_t l = (uint32_t)((x->d + x->dphase) >> 1) * max + DOSAGE_4TH;
            uint32_t r = (uint32_t)((x->d - x->dphase) >> 1) * max + DOSAGE_4TH;
            /* plink2 skips this test at 15 and 16 bits, where any dosage
             * phase already passes it */
            if ((l ^ r) / DOSAGE_MID) dphase_on = true;
        }
    }
    phased = phased || dphase_on;
    put_u32(b, (uint32_t)n);
    put_u16(b, 2);
    kputc_(2, b);
    kputc_(2, b);
    for (int j = 0; j < n; ++j) kputc_((char)(2 | (geno[j].g == 3 && !geno[j].has_dosage ? 0x80 : 0)), b);
    kputc_(phased ? 1 : 0, b);
    kputc_((char)bw->bits, b);
    uint32_t *probs = probs_buf(sc, 2 * (size_t)n);
    for (int j = 0; j < n; ++j) {
        const plink2_geno *x = &geno[j];
        uint32_t d = (uint32_t)x->d;
        uint32_t v0, v1;
        if (!phased) {
            /* P(ALT/ALT), P(ALT/REF) */
            if (!x->has_dosage) {
                v0 = x->g == 2 ? max : 0;
                v1 = x->g == 1 ? max : 0;
            } else if (d > DOSAGE_MID) {
                v0 = scaled(d - DOSAGE_MID, max);
                v1 = max - v0;
            } else {
                v0 = 0;
                v1 = scaled(d, max);
            }
        } else if (!x->has_dosage) {
            /* P(hap1 = ALT), P(hap2 = ALT) */
            if (x->g == 1) {
                if (x->hp) {
                    v0 = x->pi ? max : 0;
                    v1 = x->pi ? 0 : max;
                } else {
                    v0 = v1 = half;
                }
            } else {
                v0 = v1 = x->g == 2 ? max : 0;
            }
        } else if (dphase_on && x->has_dphase) {
            v0 = scaled((uint32_t)((x->d + x->dphase) >> 1), max);
            v1 = scaled((uint32_t)((x->d - x->dphase) >> 1), max);
        } else if (x->hp) {
            if (d > DOSAGE_MID) {
                v0 = scaled(d - DOSAGE_MID, max);
                v1 = max;
            } else {
                v0 = 0;
                v1 = scaled(d, max);
            }
            if (x->pi) {
                uint32_t t = v0;
                v0 = v1;
                v1 = t;
            }
        } else {
            v0 = v1 = (d * max + DOSAGE_MID) / DOSAGE_MAX;
        }
        probs[2 * j] = v0;
        probs[2 * j + 1] = v1;
    }
    put_bits(b, probs, 2 * (size_t)n, bw->bits);
}

void bgen_quantise(const float *p, int n, uint32_t max, uint32_t *out) {
    uint32_t q[n];
    double rem[n];
    uint32_t sum = 0;
    for (int a = 0; a < n; ++a) {
        if (!isfinite(p[a])) util_exit(PROGRAM ": bgen=phased: cannot encode a non-finite allele probability");
        double v = (double)p[a] * max;
        q[a] = (uint32_t)v;
        rem[a] = v - q[a];
        sum += q[a];
    }
    while (sum < max) {
        int best = 0;
        for (int a = 1; a < n; ++a) {
            if (rem[a] > rem[best]) best = a;
        }
        ++q[best];
        rem[best] = -1.0;
        ++sum;
    }
    for (int a = 0; a + 1 < n; ++a) out[a] = q[a];
}

/* Deflates sc->block into sc->zbuf and returns the compressed length. */
static size_t compress_block(bgen_scratch *sc) {
    size_t bound = libdeflate_zlib_compress_bound(sc->compressor, sc->block.l);
    if (bound > sc->zbuf_len) {
        sc->zbuf = util_realloc(sc->zbuf, bound);
        sc->zbuf_len = bound;
    }
    size_t zlen = libdeflate_zlib_compress(sc->compressor, sc->block.s, sc->block.l, sc->zbuf, bound);
    if (zlen == 0) util_exit("bgen: libdeflate_zlib_compress failed");
    return zlen;
}

/* The uncompressed probability data of a bgen=phased record: layout 2,
 * phased, one value per haplotype for alleles 0 to K-2. */
static void phased_block(const bgen_writer *bw, bgen_scratch *sc, const bgen_rec *rec, int n_alleles) {
    kstring_t *b = &sc->block;
    const samples *s = bw->samples;
    int min_ploidy = 2, max_ploidy = 1;
    for (int j = 0; j < s->n; ++j) {
        int ploidy = s->is_diploid[j] ? 2 : 1;
        if (ploidy < min_ploidy) min_ploidy = ploidy;
        if (ploidy > max_ploidy) max_ploidy = ploidy;
    }
    put_u32(b, (uint32_t)s->n);
    put_u16(b, (uint32_t)n_alleles);
    kputc_((char)min_ploidy, b);
    kputc_((char)max_ploidy, b);
    for (int j = 0; j < s->n; ++j) kputc_(s->is_diploid[j] ? 2 : 1, b);
    kputc_(1, b);
    kputc_((char)bw->bits, b);
    uint32_t *probs = probs_buf(sc, 2 * (size_t)s->n * (size_t)(n_alleles - 1));
    size_t n_probs = 0;
    for (int j = 0; j < s->n; ++j) {
        const bgen_call *c = &rec->calls[j];
        for (int h = 0; h < (s->is_diploid[j] ? 2 : 1); ++h) {
            int allele = h == 0 ? c->a1 : c->a2;
            for (int a = 0; a + 1 < n_alleles; ++a) {
                probs[n_probs++] = c->hap_q != NULL ? c->hap_q[h * (n_alleles - 1) + a] : (allele == a ? bw->max : 0);
            }
        }
    }
    put_bits(b, probs, n_probs, bw->bits);
}

span bgen_info_value(span info, const char *key) {
    size_t n = strlen(key);
    const char *p = info.s, *end = info.s + info.n;
    while (p < end) {
        const char *stop = memchr(p, ';', (size_t)(end - p));
        if (stop == NULL) stop = end;
        if ((size_t)(stop - p) > n && memcmp(p, key, n) == 0 && p[n] == '=') return (span){p + n + 1, (int)(stop - p - n - 1)};
        p = stop + 1;
    }
    return (span){".", 1};
}

bool bgen_info_flag(span info, const char *key) {
    size_t n = strlen(key);
    const char *p = info.s, *end = info.s + info.n;
    while (p < end) {
        const char *stop = memchr(p, ';', (size_t)(end - p));
        if (stop == NULL) stop = end;
        if ((size_t)(stop - p) == n && memcmp(p, key, n) == 0) return true;
        p = stop + 1;
    }
    return false;
}

/* The record header sc->rec, then the zlen compressed bytes in sc->zbuf,
 * as a BGEN_REC_WRITE record, and its .info row: the fields BGEN has no
 * place for, as the VCF prints them. */
static void finish_rec(bgen_scratch *sc, size_t zlen, span info, bgen_rec *rec) {
    kstring_t *r = &sc->rec;
    if (sc->block.l > UINT32_MAX || zlen > UINT32_MAX - 4) util_exit("bgen: record too large for BGEN");
    put_u32(r, (uint32_t)(zlen + 4));
    put_u32(r, (uint32_t)sc->block.l);
    rec->kind = BGEN_REC_WRITE;
    rec->len = r->l + zlen;
    rec->bytes = util_realloc(rec->bytes, rec->len);
    memcpy(rec->bytes, r->s, r->l);
    memcpy(rec->bytes + r->l, sc->zbuf, zlen);

    const marker *mk = rec->mk;
    span id = marker_id(mk), alleles = marker_alleles(mk);
    const char *tab = memchr(alleles.s, '\t', (size_t)alleles.n);
    int ref_len = (int)(tab - alleles.s);
    span dr2 = bgen_info_value(info, "DR2"), af = bgen_info_value(info, "AF");
    ksprintf(&rec->info, "%s\t%d\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%d\n", marker_chrom(mk), mk->pos, id.n, id.s, ref_len,
            alleles.s, alleles.n - ref_len - 1, tab + 1, dr2.n, dr2.s, af.n, af.s, bgen_info_flag(info, "IMP"));
}

/* A bgen=phased record: every allele in VCF order and the chromosome as
 * the VCF names it. */
static void phased_encode(const bgen_writer *bw, bgen_scratch *sc, span info, bgen_rec *rec) {
    const marker *mk = rec->mk;
    int n_alleles = marker_n_alleles(mk);
    sc->block.l = 0;
    phased_block(bw, sc, rec, n_alleles);
    size_t zlen = compress_block(sc);

    kstring_t *r = &sc->rec;
    r->l = 0;
    span id = marker_id(mk);
    span alleles = marker_alleles(mk);
    const char *chrom = marker_chrom(mk);
    put_u16(r, 0);
    put_u16(r, (uint32_t)id.n);
    kputsn(id.s, (size_t)id.n, r);
    put_u16(r, (uint32_t)strlen(chrom));
    kputs(chrom, r);
    put_u32(r, (uint32_t)mk->pos);
    put_u16(r, (uint32_t)n_alleles);
    /* REF, a tab, then the ALT alleles separated by commas */
    const char *a = alleles.s, *end = alleles.s + alleles.n;
    for (int k = 0; k < n_alleles; ++k) {
        const char *stop = a;
        while (stop < end && *stop != '\t' && *stop != ',') ++stop;
        put_u32(r, (uint32_t)(stop - a));
        kputsn(a, (size_t)(stop - a), r);
        a = stop + 1;
    }
    finish_rec(sc, zlen, info, rec);
}

/* The bgen=plink2 encoding of rec, or its kind when plink2 would not
 * export it. */
static void plink2_encode(const bgen_writer *bw, bgen_scratch *sc, span info, bgen_rec *rec) {
    const marker *mk = rec->mk;
    if (marker_n_alleles(mk) > 2) {
        rec->kind = BGEN_REC_MULTIALLELIC;
        return;
    }
    if (bw->has_min_dr2) {
        /* InfoConditionSatisfiedInternal for "DR2 >= x": the whole value
         * must parse */
        span v = bgen_info_value(info, "DR2");
        sc->info.l = 0;
        kputsn(v.s, (size_t)v.n, &sc->info);
        double x;
        if (plink2_scanadv_double(sc->info.s, &x) != sc->info.s + v.n || !(x >= bw->min_dr2)) return;
    }
    int n = bw->samples->n;
    plink2_geno *geno = sc->geno;
    for (int j = 0; j < n; ++j) geno[j] = import_call(&rec->calls[j], bw->samples->is_diploid[j]);
    if (bw->min_maf != 0.0 && maf(geno, n) < bw->min_maf) return;

    sc->block.l = 0;
    probability_block(bw, sc, geno, n);
    size_t zlen = compress_block(sc);

    kstring_t *r = &sc->rec;
    r->l = 0;
    span id = marker_id(mk);
    span alleles = marker_alleles(mk);
    const char *tab = memchr(alleles.s, '\t', (size_t)alleles.n);
    int ref_len = (int)(tab - alleles.s);
    char chrom[12];
    snprintf(chrom, sizeof chrom, "%d", autosome_code(bw, marker_chrom(mk)));
    put_u16(r, 0);
    put_u16(r, (uint32_t)id.n);
    kputsn(id.s, (size_t)id.n, r);
    put_u16(r, (uint32_t)strlen(chrom));
    kputs(chrom, r);
    put_u32(r, (uint32_t)mk->pos);
    put_u16(r, 2);
    put_u32(r, (uint32_t)(alleles.n - ref_len - 1));  /* ALT first */
    kputsn(tab + 1, (size_t)(alleles.n - ref_len - 1), r);
    put_u32(r, (uint32_t)ref_len);
    kputsn(alleles.s, (size_t)ref_len, r);
    finish_rec(sc, zlen, info, rec);

    for (int j = 0; j < n; ++j) rec->n_missing += geno[j].g == 3 && !geno[j].has_dosage;
    if (rec->n_missing > 0) {
        rec->missing = util_realloc(rec->missing, (size_t)rec->n_missing * sizeof *rec->missing);
        for (int j = 0, k = 0; j < n; ++j) {
            if (geno[j].g == 3 && !geno[j].has_dosage) rec->missing[k++] = j;
        }
    }
}

void bgen_writer_encode(const bgen_writer *bw, bgen_scratch *sc, span info, bgen_rec *rec) {
    rec->kind = BGEN_REC_FILTERED;
    rec->len = 0;
    rec->info.l = 0;
    rec->n_missing = 0;
    if (bw->mode == BGEN_PHASED) phased_encode(bw, sc, info, rec);
    else plink2_encode(bw, sc, info, rec);
    release_calls(rec);
}

void bgen_writer_put(bgen_writer *bw, bgen_rec *rec) {
    const marker *mk = rec->mk;
    if (bw->mode == BGEN_PLINK2 && mk->chrom_index != bw->chrom_index) {
        bgen_writer_check_chrom(bw, marker_chrom(mk));
        bw->chrom_index = mk->chrom_index;
    }
    if (rec->kind == BGEN_REC_MULTIALLELIC) ++bw->n_multiallelic;
    if (rec->kind == BGEN_REC_WRITE) {
        write_bytes(bw, rec->bytes, rec->len);
        if (fwrite(rec->info.s, 1, rec->info.l, bw->info_file) != rec->info.l) util_exit("Error writing %s", bw->info_path);
        ++bw->n_variants;
        for (int k = 0; k < rec->n_missing; ++k) ++bw->n_missing[rec->missing[k]];
    }
}

/* ExportOxSample: no phenotypes, FID 0, sex unknown. */
static void write_sample_file(bgen_writer *bw) {
    FILE *f = bgen_files_open(bw->sample_path);
    if (f == NULL) util_exit("Error opening %s", bw->sample_path);
    fputs("ID_1 ID_2 missing sex\n0 0 0 D\n", f);
    double recip = bw->n_variants == 0 ? 0.0 : 1.0 / (double)bw->n_variants;
    char buf[16];
    for (int j = 0; j < bw->samples->n; ++j) {
        plink2_dtoa_g_unit((double)bw->n_missing[j] * recip, buf);
        fprintf(f, "0 %s %s NA\n", bw->samples->ids[j], buf);
    }
    bool write_failed = ferror(f);
    if (fclose(f) != 0 || write_failed) util_exit("Error writing %s", bw->sample_path);
}

void bgen_writer_close(bgen_writer *bw) {
    if (bw->n_multiallelic > 0) {
        fprintf(stderr, "bgen=plink2: skipped %ld multiallelic record%s, as plink2 --import-max-alleles 2 would\n",
                bw->n_multiallelic, bw->n_multiallelic == 1 ? "" : "s");
    }
    if (bw->mode == BGEN_PLINK2 && bw->n_variants == 0) util_exit("bgen=plink2: no variants remaining after the bgen filters");
    bw->rec.l = 0;
    put_u32(&bw->rec, bw->n_variants);
    if (fseek(bw->bgen, 8, SEEK_SET) != 0) util_exit("Error writing %s", bw->bgen_path);
    write_bytes(bw, bw->rec.s, bw->rec.l);
    if (fclose(bw->bgen) != 0) util_exit("Error writing %s", bw->bgen_path);
    if (fclose(bw->info_file) != 0) util_exit("Error writing %s", bw->info_path);
    write_sample_file(bw);
    bgen_files_complete();
    free(bw->rec.s);
    free(bw->n_missing);
    free(bw);
}
