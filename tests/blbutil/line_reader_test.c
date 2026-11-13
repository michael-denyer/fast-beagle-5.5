/* line_reader reads past a run of empty BGZF blocks. Beagle 5.5's BGZipIt
 * stops at a batch of 4 x nthreads blocks that all inflate to nothing and
 * drops the rest of the file; fast-beagle reads it all
 * (docs/beagle-divergences.md#input-files). */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include "blbutil/line_reader.h"

/* The empty BGZF block that ends a BGZF file. */
static const unsigned char EMPTY_BLOCK[28] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 0x42, 0x43, 2, 0,
        0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/* One BGZF member holding text, ending with its own empty block. */
static void append_member(const char *path, const char *text) {
    const char *tmp = "build/blbutil/line_reader_test_member.gz";
    BGZF *z = bgzf_open(tmp, "w");
    assert(z != NULL);
    assert(bgzf_write(z, text, strlen(text)) == (ssize_t)strlen(text));
    assert(bgzf_close(z) == 0);
    FILE *in = fopen(tmp, "rb");
    FILE *out = fopen(path, "ab");
    assert(in != NULL && out != NULL);
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, in)) > 0) assert(fwrite(buf, 1, n, out) == n);
    assert(fclose(in) == 0 && fclose(out) == 0);
    remove(tmp);
}

int main(void) {
    const char *path = "build/blbutil/line_reader_test.vcf.gz";
    remove(path);
    append_member(path, "1\n");
    FILE *out = fopen(path, "ab");
    assert(out != NULL);
    for (int j = 0; j < 16; ++j) assert(fwrite(EMPTY_BLOCK, 1, sizeof EMPTY_BLOCK, out) == sizeof EMPTY_BLOCK);
    assert(fclose(out) == 0);
    append_member(path, "2\n");

    for (int n_threads = 1; n_threads <= 4; n_threads += 3) {
        line_reader *r = line_reader_open(path, n_threads);
        kstring_t line = {0, 0, NULL};
        assert(line_reader_next(r, &line) && strcmp(line.s, "1") == 0);
        assert(line_reader_next(r, &line) && strcmp(line.s, "2") == 0);
        assert(!line_reader_next(r, &line));
        assert(!line_reader_next(r, &line));
        free(line.s);
        line_reader_close(r);
    }
    remove(path);
    puts("PASS line_reader reads past 16 empty BGZF blocks at 1 and 4 threads");
    return 0;
}
