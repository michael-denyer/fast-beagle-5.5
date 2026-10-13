/* The window cutter over an in-memory source of target records. */
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "blbutil/utilities.h"
#include "vcf/gt_rec.h"
#include "vcf/sliding_window.h"
#include "vcf/vcf_header.h"

typedef struct {
    gt_rec **recs;
    int n, n_read;
    const samples *smp;
} array_it;

static const samples *array_samples(const void *self) {
    return ((const array_it *)self)->smp;
}

static void *array_next(void *self) {
    array_it *it = self;
    return it->n_read < it->n ? it->recs[it->n_read++] : NULL;
}

static void array_close(void *self) {
    array_it *it = self;
    for (int j = it->n_read; j < it->n; ++j) gt_rec_release(it->recs[j]);
}

static const marker *rec_marker(const void *r) {
    return &((const gt_rec *)r)->marker;
}

static void rec_release(void *r) {
    gt_rec_release(r);
}

static const sample_file_it_ops array_ops = {.samples = array_samples, .next = array_next, .close = array_close,
        .marker = rec_marker, .release = rec_release};

static vcf_header header;

static gt_rec *rec_at(const char *chrom, int32_t pos) {
    char line[128];
    int n = snprintf(line, sizeof line, "%s\t%d\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|1", chrom, pos);
    gt_rec *rec = util_malloc(sizeof *rec);
    gt_rec_parse(rec, line, (size_t)n, &header);
    rec->refs = 1;
    return rec;
}

/* window=2 cM and overlap=0.5 cM on the default map of 1 cM per Mb. */
static sliding_window *open_windows(gt_rec **recs, int n, array_it *src) {
    static const par p = {.window = 2.0f, .overlap = 0.5f, .window_markers = 1000000, .impute = true};
    *src = (array_it){.recs = recs, .n = n, .smp = &header.samples};
    return sliding_window_open_its(&p, (sample_file_it){&array_ops, src}, (sample_file_it){0}, genetic_map_open(NULL, NULL));
}

static void check_window(const window *w, int index, bool last, int chrom, int n_targ, int prev_splice, int overlap_end,
        int overlap_start, int next_splice, int32_t first_pos, int32_t last_pos) {
    assert(w->index == index);
    assert(w->last_window == last);
    assert(w->chrom_index == chrom);
    assert(w->n_targ == n_targ && w->n_ref == 0);
    assert(w->indices.n_targ_markers == n_targ && w->indices.n_markers == n_targ);
    assert(w->indices.targ_prev_splice == prev_splice && w->indices.prev_splice == prev_splice);
    assert(w->indices.targ_overlap_end == overlap_end && w->indices.overlap_end == overlap_end);
    assert(w->indices.targ_overlap_start == overlap_start && w->indices.overlap_start == overlap_start);
    assert(w->indices.targ_next_splice == next_splice && w->indices.next_splice == next_splice);
    assert(w->targ[0]->marker.pos == first_pos);
    assert(w->targ[n_targ - 1]->marker.pos == last_pos);
}

int main(void) {
    char *lines[] = {"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"};
    const char *first = "1\t100000\t.\tA\tC\t.\t.\t.\tGT\t0|1\t1|1";
    vcf_header_init(&header, "test.vcf", lines, 1, first, strlen(first), NULL);
    assert(header.samples.n == 2);

    /* Window 1 ends 2 cM past the first marker, at 2.1 Mb. Its overlap starts
     * at the last marker before 0.5 cM back from its last marker (1.4 Mb). Window 2
     * ends 1.5 cM past the first new marker, at 4.0 Mb, and carries the markers
     * from 2.5 Mb into window 3, which takes the rest and is the last window.
     * The splice points halve each overlap, so window 1 prints its markers 0 to
     * 2 and windows 2 and 3 their markers 1 to 3: each of the nine once. */
    gt_rec *recs[9];
    const int32_t pos[9] = {100000, 500000, 1000000, 1500000, 1900000, 2500000, 3100000, 3500000, 4000000};
    for (int j = 0; j < 9; ++j) recs[j] = rec_at("1", pos[j]);
    array_it src;
    sliding_window *sw = open_windows(recs, 9, &src);
    assert(sliding_window_targ_samples(sw)->n == 2);
    assert(sliding_window_n_ref_samples(sw) == 0);
    window *w = sliding_window_next(sw);
    check_window(w, 1, false, 0, 5, 0, 0, 2, 3, 100000, 1900000);
    window_free(w);
    w = sliding_window_next(sw);
    check_window(w, 2, false, 0, 6, 1, 3, 3, 4, 1000000, 3500000);
    assert(w->targ[3]->marker.pos == 2500000);
    window_free(w);
    w = sliding_window_next(sw);
    check_window(w, 3, true, 0, 4, 1, 3, 4, 4, 2500000, 4000000);
    window_free(w);
    assert(sliding_window_next(sw) == NULL);
    assert(sliding_window_cum_targ_markers(sw) == 9 && sliding_window_cum_markers(sw) == 9);
    sliding_window_close(sw);
    puts("PASS windows end window= cM past their first new marker and overlap by overlap= cM");

    /* A chromosome break ends the window with no overlap, and the next
     * chromosome's first window is a full window= cM again. */
    gt_rec *two[4] = {rec_at("1", 100000), rec_at("1", 500000), rec_at("2", 100000), rec_at("2", 300000)};
    sw = open_windows(two, 4, &src);
    w = sliding_window_next(sw);
    check_window(w, 1, false, 0, 2, 0, 0, 2, 2, 100000, 500000);
    window_free(w);
    w = sliding_window_next(sw);
    check_window(w, 2, true, 1, 2, 0, 0, 2, 2, 100000, 300000);
    window_free(w);
    assert(sliding_window_next(sw) == NULL);
    assert(sliding_window_cum_targ_markers(sw) == 4);
    sliding_window_close(sw);
    puts("PASS a chromosome break ends the window without an overlap");

    vcf_header_free(&header);
    return 0;
}
