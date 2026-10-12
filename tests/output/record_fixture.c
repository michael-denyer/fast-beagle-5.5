/* Exercise the same record construction interface as both output producers,
 * then the phased driver over enough markers to wrap its ring of records and
 * leave a partial last task. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main/window_writer.h"

#define PHASED_MARKERS 1013

/* Haplotype hap's allele at marker m: one bit of m per haplotype, so rows and
 * samples differ. tests/check_records.py recomputes it. */
static int allele(const void *ctx, int m, int hap) {
    (void)ctx;
    return (m >> (hap % 4)) & 1;
}

int main(int argc, char **argv) {
    if (argc != 4) return 2;
    char *ids[] = {"D1", "H", "D2"};
    bool diploid[] = {true, false, true};
    samples s = {3, ids, diploid};
    par p = {.out = argv[1], .nthreads = atoi(argv[3]), .ap = true, .gp = true,
        .bgen = (bgen_mode)atoi(argv[2]), .bgen_bits = 16, .bgen_chr_set = 22};
    run_outputs *out = run_outputs_new(&p);
    window_writer ww;
    window_writer_open(&ww, &p, out, &s);
    out_worker *wk = window_writer_worker_new(&ww);
    out_rec r = {0};
    const char *lines[] = {
        "22\t100\tphased\tA\tC\t.\t.\t.\tGT",
        "22\t200\timputed\tA\tC\t.\t.\t.\tGT",
        "22\t300\tobserved\tA\tC,G\t.\t.\tEND=301\tGT",
        "22\t400\tmono\tA\t.\t.\t.\t.\tGT",
        "22\t500\tfive\tA\tC,G,T,AC\t.\t.\t.\tGT"
    };
    for (int m = 0; m < 5; ++m) {
        marker mk;
        marker_parse(&mk, lines[m], strlen(lines[m]));
        out_rec_kind kind = m == 0 ? OUT_PHASED : m == 2 ? OUT_GENOTYPED : OUT_IMPUTED;
        window_writer_rec_begin(&ww, &r, &mk, kind);
        for (int j = 0; j < s.n; ++j) {
            if (m == 0) {
                window_writer_rec_gt(&r, j == 1 ? 1 : 0, 1);
            } else {
                float a1[5] = {1}, a2[5] = {1};
                if (m == 1 && j == 0) {
                    a1[0] = 0.6662f; a1[1] = 0.3338f;
                    a2[0] = 0.1234f; a2[1] = 0.8766f;
                } else if (m == 1 && j == 1) {
                    a1[0] = 3; a1[1] = 1;
                } else if (m == 2 && j < 2) {
                    a1[0] = 0; a1[j == 0 ? 2 : 1] = 1;
                }
                window_writer_rec_probs(&r, a1, a2);
                /* Both outputs must already own what they need. */
                memset(a1, 0, sizeof a1);
                memset(a2, 0, sizeof a2);
            }
        }
        window_writer_encode(wk, &r);
        window_writer_put(&ww, &r);
        marker_free(&mk);
    }
    window_writer_rec_free(&r);
    window_writer_worker_free(wk);

    gt_rec *recs = calloc(PHASED_MARKERS, sizeof *recs);
    gt_rec **targ = malloc(PHASED_MARKERS * sizeof *targ);
    for (int m = 0; m < PHASED_MARKERS; ++m) {
        char line[64];
        int len = snprintf(line, sizeof line, "22\t%d\tm%d\tA\tC\t.\t.\t.\tGT", 1000 + m, m);
        marker_parse(&recs[m].marker, line, (size_t)len);
        targ[m] = &recs[m];
    }
    window w = {.n_targ = PHASED_MARKERS, .targ = targ};
    window_writer_print_phased(&ww, &w, 0, PHASED_MARKERS, allele, NULL);
    for (int m = 0; m < PHASED_MARKERS; ++m) marker_free(&recs[m].marker);
    free(targ);
    free(recs);

    window_writer_close(&ww);
    run_outputs_free(out);
    return 0;
}
