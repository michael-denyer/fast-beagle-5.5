/*
 * Copyright (C) 2014-2021 Brian L. Browning
 * Ported to C from Beagle 5.5 (27Feb25) vcf/Marker.java, vcf/MarkerParser.java
 * and vcf/MarkerUtils.java; modified 2026.
 *
 * This file is part of fast-beagle, a C port of Beagle. It is free software:
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version. See LICENSE.
 */
#include "vcf/marker.h"

#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "beagleutil/chrom_ids.h"
#include "blbutil/utilities.h"
#include "jcompat/jnum.h"

#include <htslib/kstring.h>

#define ID_STORED ((uint16_t)(1 << 15))
#define ALLELES_STORED ((uint16_t)(1 << 14))
#define END_STORED ((uint16_t)(1 << 10))
#define STORED_N_ALLELES_MASK 0xff
#define INDEXED_N_ALLELES_MASK 0x7
#define SNV_INDEX_MASK 0x7f

#define N_SNV_PERMS 100
#define MAX_SNV_ALLELES 5

/* MarkerUtils.snvPerms(): REF '\t' ALT for every ordering of {*, A, C, G, T}
 * whose REF is not '*', plus the four REF-only records, sorted. */
static char snv_perms[N_SNV_PERMS][2 * MAX_SNV_ALLELES];
static int n_snv_perms;

static void add_snv_perm(const char a[MAX_SNV_ALLELES]) {
    if (a[0] == '*') return;
    char *s = snv_perms[n_snv_perms++];
    int k = 0;
    s[k++] = a[0];
    for (int j = 1; j < MAX_SNV_ALLELES; ++j) {
        s[k++] = j == 1 ? '\t' : ',';
        s[k++] = a[j];
    }
    s[k] = '\0';
}

static int compare_str(const void *a, const void *b) {
    return strcmp((const char *)a, (const char *)b);
}

static void init_snv_perms(void) {
    /* Heap's algorithm visits every ordering; the qsort below fixes their order. */
    char a[MAX_SNV_ALLELES] = {'*', 'A', 'C', 'G', 'T'};
    int c[MAX_SNV_ALLELES] = {0};
    add_snv_perm(a);
    for (int i = 1; i < MAX_SNV_ALLELES;) {
        if (c[i] < i) {
            int j = i % 2 == 0 ? 0 : c[i];
            char t = a[j];
            a[j] = a[i];
            a[i] = t;
            add_snv_perm(a);
            ++c[i];
            i = 1;
        } else {
            c[i++] = 0;
        }
    }
    for (const char *b = "ACGT"; *b != '\0'; ++b) {
        char *s = snv_perms[n_snv_perms++];
        s[0] = *b;
        s[1] = '\t';
        s[2] = '.';
        s[3] = '\0';
    }
    qsort(snv_perms, N_SNV_PERMS, sizeof snv_perms[0], compare_str);
}

/* String.compareTo of a table row with the len bytes at s. s is a Java String,
 * so it can hold NUL bytes and has no terminator. */
static int compare_perm(const char *perm, const char *s, size_t len) {
    size_t perm_len = strlen(perm);
    int c = memcmp(perm, s, perm_len < len ? perm_len : len);
    return c != 0 ? c : (perm_len > len) - (perm_len < len);
}

/* MarkerParser.snvIndex: Arrays.binarySearch gives the match or the insertion
 * point, which for distinct sorted entries is the first entry >= s. */
static int snv_index(const char *s, size_t len) {
    int lo = 0, hi = N_SNV_PERMS;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (compare_perm(snv_perms[mid], s, len) < 0) lo = mid + 1;
        else hi = mid;
    }
    if (lo == N_SNV_PERMS) return -1;
    return strlen(snv_perms[lo]) >= len && memcmp(snv_perms[lo], s, len) == 0 ? lo : -1;
}

static int truncate80(size_t len) {
    return len < 80 ? (int)len : 80;
}

static const char *index_of(const char *s, const char *end, char c) {
    const void *p = memchr(s, c, (size_t)(end - s));
    return p;
}

/* MarkerParser.nAlleles for alleles stored as text. */
static int n_text_alleles(const char *alt, const char *end) {
    if (end - alt == 1 && alt[0] == '.') return 1;
    int n = 2;
    for (const char *c = alt; c < end; ++c) n += *c == ',';
    return n;
}

static bool is_java_whitespace(char c) {
    return c == ' ' || (c >= '\t' && c <= '\r') || (c >= 0x1c && c <= 0x1f);
}

static int32_t parse_pos(const char *s, const char *end) {
    int32_t pos;
    if (!jnum_parse_int(s, (size_t)(end - s), &pos)) {
        util_exit("java.lang.NumberFormatException: For input string: \"%.*s\"", (int)(end - s), s);
    }
    return pos;
}

void marker_parse(marker *m, const char *rec, size_t len) {
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    pthread_once(&once, init_snv_perms);
    const char *rec_end = rec + len;
    const char *tabs[8];
    const char *p = rec;
    int n_tabs = 0;
    while (n_tabs < 8 && (p = index_of(p, rec_end, '\t')) != NULL) tabs[n_tabs++] = p++;
    if (n_tabs < 8) {
        int n = len < 800 ? (int)len : 800;
        util_exit("VCF record does not contain 8 tabs:%.*s", n, rec);
    }

    /* MarkerUtils.chromIndex */
    size_t chrom_len = (size_t)(tabs[0] - rec);
    if (chrom_len == 0 || (chrom_len == 1 && rec[0] == '.')) {
        util_exit("ERROR: missing chromosome: %.*s", truncate80(len), rec);
    }
    for (size_t j = 0; j < chrom_len; ++j) {
        if (is_java_whitespace(rec[j])) {
            util_exit("ERROR: CHROM field contains whitespace: %.*s", truncate80(len), rec);
        }
    }
    m->chrom_index = chrom_ids_index(rec, chrom_len);
    m->pos = parse_pos(tabs[0] + 1, tabs[1]);

    if (m->chrom_index >= INT16_MAX) util_exit("java.lang.IndexOutOfBoundsException: %d", m->chrom_index);

    /* MarkerParser.storeMarkerFields: ID, then REF and ALT, then INFO/END, each
     * appended to `fields` with a tab when `fields` is not empty. */
    kstring_t sb = {0, 0, NULL};
    uint16_t info = 0;
    const char *id = tabs[1] + 1;
    size_t id_len = (size_t)(tabs[2] - id);
    if (!(id_len == 1 && id[0] == '.')) {
        kputsn(id, id_len, &sb);
        info |= ID_STORED;
    }

    const char *alleles = tabs[2] + 1;
    size_t alleles_len = (size_t)(tabs[4] - alleles);
    int snv = snv_index(alleles, alleles_len);
    if (snv >= 0) {
        bool ref_only = alleles_len >= 2 && alleles[alleles_len - 2] == '\t' && alleles[alleles_len - 1] == '.';
        info |= (uint16_t)(snv << 3);
        info |= (uint16_t)(ref_only ? 1 : (int)((alleles_len + 1) >> 1));
    } else {
        if (tabs[3] + 1 == tabs[4]) {
            util_exit("ERROR: missing ALT field: %.*s", truncate80(len), rec);
        }
        int n_alleles = n_text_alleles(tabs[3] + 1, tabs[4]);
        if (n_alleles > STORED_N_ALLELES_MASK) {
            util_exit("java.lang.IndexOutOfBoundsException: %d alleles: %.*s", n_alleles, truncate80(len), rec);
        }
        if (sb.l > 0) kputc('\t', &sb);
        kputsn(alleles, alleles_len, &sb);
        info |= (uint16_t)n_alleles;
        info |= ALLELES_STORED;
    }

    /* MarkerParser.storeInfo with storeInfo false: the first ';'-separated
     * subfield starting "END=", if the INFO field is not ".". */
    const char *inf = tabs[6] + 1;
    const char *inf_end = tabs[7];
    if (!(inf_end - inf == 1 && inf[0] == '.')) {
        const char *f = inf;
        while (f <= inf_end) {
            const char *semi = index_of(f, inf_end, ';');
            const char *f_end = semi == NULL ? inf_end : semi;
            if (f_end - f >= 4 && memcmp(f, "END=", 4) == 0) {
                if (sb.l > 0) kputc('\t', &sb);
                kputsn(f, (size_t)(f_end - f), &sb);
                info |= END_STORED;
                break;
            }
            f = f_end + 1;
        }
    }
    m->field_info = info;
    m->fields_len = sb.l;
    if (sb.l == 0) {
        free(sb.s);
        m->fields = NULL;
    } else {
        m->fields = sb.s;
    }
}

void marker_free(marker *m) {
    free(m->fields);
    m->fields = NULL;
}

const char *marker_chrom(const marker *m) {
    return chrom_ids_id(m->chrom_index);
}

/* String.indexOf(ch, from) on fields. */
static long find(const marker *m, char c, long from) {
    if (from < 0) from = 0;
    if ((size_t)from >= m->fields_len) return -1;
    const char *p = memchr(m->fields + from, c, m->fields_len - (size_t)from);
    return p == NULL ? -1 : p - m->fields;
}

static const char *fields_or_npe(const marker *m) {
    if (m->fields == NULL) util_exit("java.lang.NullPointerException: Marker.fields is null");
    return m->fields;
}

/* String.substring(start, end) on fields. */
static span substring(const marker *m, long start, long end) {
    if (start < 0 || end > (long)m->fields_len || start > end) {
        util_exit("java.lang.StringIndexOutOfBoundsException: begin %ld, end %ld, length %zu", start, end, m->fields_len);
    }
    return (span){m->fields + start, (int)(end - start)};
}

int marker_n_alleles(const marker *m) {
    if (m->field_info & ALLELES_STORED) return m->field_info & STORED_N_ALLELES_MASK;
    return m->field_info & INDEXED_N_ALLELES_MASK;
}

bool marker_has_id(const marker *m) {
    return (m->field_info & ID_STORED) != 0;
}

span marker_id(const marker *m) {
    if (!(m->field_info & ID_STORED)) return (span){".", 1};
    fields_or_npe(m);
    long end = find(m, '\t', 0);
    return substring(m, 0, end < 0 ? (long)m->fields_len : end);
}

span marker_alleles(const marker *m) {
    if (m->field_info & ALLELES_STORED) {
        long start = 0;
        fields_or_npe(m);
        if (m->field_info & ID_STORED) start = find(m, '\t', 0) + 1;
        long end = find(m, '\t', start);
        end = find(m, '\t', end + 1);
        return substring(m, start, end < 0 ? (long)m->fields_len : end);
    }
    int snv = (m->field_info >> 3) & SNV_INDEX_MASK;
    int n = m->field_info & INDEXED_N_ALLELES_MASK;
    const char *perm = snv_perms[snv];
    return (span){perm, n == 1 ? (int)strlen(perm) : 2 * n - 1};
}

/* Marker.qualStartIndex: skip the stored ID, REF and ALT. */
static long qual_start(const marker *m) {
    long start = 0;
    if (m->field_info & ID_STORED) start = find(m, '\t', 0) + 1;
    if (m->field_info & ALLELES_STORED) {
        start = find(m, '\t', start) + 1;
        start = find(m, '\t', start) + 1;
    }
    return start;
}

span marker_end_value(const marker *m) {
    if (!(m->field_info & END_STORED)) return (span){"", 0};
    fields_or_npe(m);
    long from = qual_start(m);
    const char *hit = NULL;
    for (long j = from < 0 ? 0 : from; j + 4 <= (long)m->fields_len; ++j) {
        if (memcmp(m->fields + j, "END=", 4) == 0) {
            hit = m->fields + j;
            break;
        }
    }
    long start = (hit == NULL ? -1 : hit - m->fields) + 4;  /* Java: indexOf("END=", from) + 4 */
    long end = find(m, ';', start);
    return substring(m, start, end < 0 ? (long)m->fields_len : end);
}

int marker_bits_per_allele(const marker *m) {
    int n_alleles = marker_n_alleles(m);
    int bits = 0;
    while ((1 << bits) < n_alleles) ++bits;
    return bits;
}

/* Marker.info: QUAL and FILTER are not stored, so INFO starts where QUAL would. */
span marker_info(const marker *m) {
    if (!(m->field_info & END_STORED)) return (span){".", 1};
    fields_or_npe(m);
    long start = qual_start(m);
    long end = find(m, '\t', start);
    return substring(m, start, end < 0 ? (long)m->fields_len : end);
}
