/*
 * vcf_index against htslib's tbx_index_build3 (tabix -p vcf): the two .tbi
 * files must be byte-identical, for VCFs written as window_writer writes them
 * by one thread and by several.
 *
 * Usage: vcf_index_test <dir>                  the synthetic cases, in <dir>
 *        vcf_index_test tabix <vcf.gz> <tbi>   writes tabix's index to <tbi>
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#include "main/vcf_index.h"

static int failures;

typedef struct {
    BGZF *out;
    vcf_index *index;
    uint64_t n_written;
    kstring_t line;
    uint32_t rng;
} writer;

static void write_line(writer *w) {
    kputc('\n', &w->line);
    if (bgzf_write(w->out, w->line.s, w->line.l) < 0) exit(2);
    w->n_written += w->line.l;
    w->line.l = 0;
}

static void header(writer *w, size_t pad_to) {
    kputs("##fileformat=VCFv4.2", &w->line);
    write_line(w);
    const char *chrom = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3";
    if (pad_to > 0) {
        size_t n = pad_to - w->n_written - strlen(chrom) - 1;
        kputs("##x=", &w->line);
        for (size_t j = 5; j < n; ++j) kputc('x', &w->line);
        write_line(w);
    }
    kputs(chrom, &w->line);
    write_line(w);
    if (pad_to > 0 && w->n_written != pad_to) exit(2);
}

static uint32_t next_rand(writer *w) {
    w->rng = w->rng * 1103515245u + 12345u;
    return w->rng >> 16;
}

/* A record with n_samples sample fields. line_len > 0 pads ID so that the
 * line, with its newline, is line_len bytes long. */
static void record(writer *w, const char *chrom, int pos, const char *ref, const char *info, int n_samples, size_t line_len) {
    kstring_t *l = &w->line;
    ksprintf(l, "%s\t%d\t", chrom, pos);
    size_t id_at = l->l;
    ksprintf(l, "id%d\t%s\tA\t.\tPASS\t%s", pos, ref, info);
    size_t info_end = l->l;
    kputs("\tGT:DS", l);
    for (int s = 0; s < n_samples; ++s) {
        uint32_t r = next_rand(w);
        ksprintf(l, "\t%u|%u:%u.%02u", r & 1, r >> 1 & 1, r >> 2 & 1, r % 100);
    }
    if (line_len > 0) {
        size_t pad = line_len - l->l - 1;
        if (line_len < l->l + 1) exit(2);
        ks_resize(l, l->l + pad + 1);
        memmove(l->s + id_at + pad, l->s + id_at, l->l - id_at + 1);
        memset(l->s + id_at, 'p', pad);
        l->l += pad;
        info_end += pad;
    }
    vcf_index_add(w->index, l->s, info_end, w->n_written + l->l + 1);
    write_line(w);
}

/* A record that ends the line exactly at the next multiple of the BGZF block
 * size. */
static void record_to_boundary(writer *w, const char *chrom, int pos) {
    size_t to = (w->n_written / BGZF_BLOCK_SIZE + 1) * BGZF_BLOCK_SIZE;
    if (to - w->n_written < 200) record(w, chrom, pos++, "C", ".", 3, 0);
    to = (w->n_written / BGZF_BLOCK_SIZE + 1) * BGZF_BLOCK_SIZE;
    record(w, chrom, pos, "C", ".", 3, to - w->n_written);
}

static void open_writer(writer *w, const char *path, int nthreads) {
    *w = (writer){.out = bgzf_open(path, "w"), .rng = 7};
    if (w->out == NULL || (nthreads > 1 && bgzf_mt(w->out, nthreads, 256) != 0)) exit(2);
}

static void start_index(writer *w, const char *path) {
    kstring_t tbi = {0, 0, NULL};
    ksprintf(&tbi, "%s.tbi", path);
    w->index = vcf_index_new(path, tbi.s, w->n_written);
    free(tbi.s);
}

static void close_writer(writer *w) {
    if (bgzf_close(w->out) != 0) exit(2);
    vcf_index_write(w->index, w->n_written);
    free(w->line.s);
}

static char *slurp(const char *path, size_t *n) {
    FILE *f = fopen(path, "rb");
    if (f == NULL) return NULL;
    kstring_t s = {0, 0, NULL};
    char buf[65536];
    size_t k;
    while ((k = fread(buf, 1, sizeof buf, f)) > 0) kputsn(buf, k, &s);
    fclose(f);
    *n = s.l;
    return s.s != NULL ? s.s : calloc(1, 1);
}

static void check(const char *name, const char *path, int nthreads) {
    kstring_t ours = {0, 0, NULL}, tabix = {0, 0, NULL};
    ksprintf(&ours, "%s.tbi", path);
    ksprintf(&tabix, "%s.tabix.tbi", path);
    size_t n1 = 0, n2 = 0;
    char *a = NULL, *b = NULL;
    if (tbx_index_build3(path, tabix.s, 0, 0, &tbx_conf_vcf) != 0) {
        printf("FAIL %s nthreads=%d: tabix cannot index %s\n", name, nthreads, path);
        ++failures;
    } else if ((a = slurp(ours.s, &n1)) == NULL || (b = slurp(tabix.s, &n2)) == NULL) {
        printf("FAIL %s nthreads=%d: missing index\n", name, nthreads);
        ++failures;
    } else if (n1 != n2 || memcmp(a, b, n1) != 0) {
        printf("FAIL %s nthreads=%d: %s differs from tabix's %s\n", name, nthreads, ours.s, tabix.s);
        ++failures;
    } else {
        printf("pass %s nthreads=%d\n", name, nthreads);
    }
    free(a);
    free(b);
    free(ours.s);
    free(tabix.s);
}

/* Chromosomes 1, 2 and 10; lines longer than a block; records ending on block
 * boundaries; a multi-base REF crossing a 16 kb linear-index window; INFO END=
 * at the start, after another entry, at or before POS (ignored), and as a
 * suffix of another key. */
static void chroms(const char *path, int nthreads) {
    writer w;
    open_writer(&w, path, nthreads);
    header(&w, 0);
    start_index(&w, path);
    int pos = 100;
    for (int j = 0; j < 400; ++j) record(&w, "1", pos += 37, "A", "DR2=0.5;AF=0.1", 40, 0);
    record_to_boundary(&w, "1", pos += 11);
    record(&w, "1", pos = 32767, "ACGTACGT", "DR2=0.9", 3, 0);
    record(&w, "1", pos += 5, "A", "END=900000", 3, 0);
    record(&w, "1", pos += 5, "A", "DR2=0.9;END=950000", 3, 0);
    record(&w, "1", pos += 5, "A", "END=5", 3, 0);
    record(&w, "1", pos += 5, "A", "XEND=990000", 3, 0);
    record(&w, "1", pos += 5, "AT", "DR2=0.1;END=.", 3, 0);
    record(&w, "1", pos += 5, "A", ".", 30000, 0);
    record_to_boundary(&w, "1", pos += 5);
    pos = 1;
    for (int j = 0; j < 3000; ++j) record(&w, "2", pos += 1000 + j, j % 7 == 0 ? "ACG" : "T", "IMP", 3, 0);
    record_to_boundary(&w, "2", pos += 50);
    record_to_boundary(&w, "10", 1);
    for (int j = 0; j < 50; ++j) record(&w, "10", 100000000 + j * 20000, "G", ".", 3, 0);
    close_writer(&w);
    check("chroms", path, nthreads);
}

/* A two-base deletion straddles a 16 kb bin boundary. Its short END must
 * not remove the bin that contains the second REF base from the index. */
static void short_end(const char *path, int nthreads) {
    writer w;
    open_writer(&w, path, nthreads);
    header(&w, 0);
    start_index(&w, path);
    record(&w, "1", 16384, "AA", "END=16384", 3, 0);
    close_writer(&w);
    tbx_t *tbx = tbx_index_load(path);
    if (tbx == NULL) exit(2);
    hts_itr_t *it = tbx_itr_queryi(tbx, 0, 16384, 16385);
    if (it == NULL || it->n_off == 0) {
        fprintf(stderr, "FAIL short END loses the second REF base at 1:16385\n");
        ++failures;
    }
    hts_itr_destroy(it);
    tbx_destroy(tbx);
}

static void header_on_boundary(const char *path, int nthreads, int n_recs) {
    writer w;
    open_writer(&w, path, nthreads);
    header(&w, BGZF_BLOCK_SIZE);
    start_index(&w, path);
    for (int j = 0; j < n_recs; ++j) record(&w, "chr5", 1000 + j, "AC", ".", 3, 0);
    close_writer(&w);
    check(n_recs == 0 ? "header-on-boundary-only" : "header-on-boundary", path, nthreads);
}

static void header_only(const char *path, int nthreads) {
    writer w;
    open_writer(&w, path, nthreads);
    header(&w, 0);
    start_index(&w, path);
    close_writer(&w);
    check("header-only", path, nthreads);
}

int main(int argc, char **argv) {
    if (argc == 4 && strcmp(argv[1], "tabix") == 0)
        return tbx_index_build3(argv[2], argv[3], 0, 0, &tbx_conf_vcf) != 0;
    if (argc != 2) {
        fprintf(stderr, "usage: vcf_index_test <dir> | vcf_index_test tabix <vcf.gz> <tbi>\n");
        return 2;
    }
    int threads[] = {1, 4};
    for (int t = 0; t < 2; ++t) {
        kstring_t path = {0, 0, NULL};
        ksprintf(&path, "%s/chroms.t%d.vcf.gz", argv[1], threads[t]);
        chroms(path.s, threads[t]);
        path.l = 0;
        ksprintf(&path, "%s/short-end.t%d.vcf.gz", argv[1], threads[t]);
        short_end(path.s, threads[t]);
        path.l = 0;
        ksprintf(&path, "%s/boundary.t%d.vcf.gz", argv[1], threads[t]);
        header_on_boundary(path.s, threads[t], 5);
        path.l = 0;
        ksprintf(&path, "%s/boundary-only.t%d.vcf.gz", argv[1], threads[t]);
        header_on_boundary(path.s, threads[t], 0);
        path.l = 0;
        ksprintf(&path, "%s/header.t%d.vcf.gz", argv[1], threads[t]);
        header_only(path.s, threads[t]);
        free(path.s);
    }
    return failures != 0;
}
