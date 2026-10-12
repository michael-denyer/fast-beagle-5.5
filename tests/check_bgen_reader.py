# /// script
# requires-python = ">=3.10,<3.13"
# dependencies = ["bgen-reader==4.0.9"]
# [tool.uv]
# exclude-newer = "2026-09-26T00:00:00Z"
# ///
"""Check that bgen-reader decodes each BGEN as tests/check_bgen_phased.py does.

Usage: uv run --python 3.12 --script check_bgen_reader.py <file>.bgen...

bgen-reader (with its cbgen backend) must find the same sample IDs, the same
variants in order (ID, rsid, chromosome, position, alleles), and per variant
and sample the same phased flag, ploidy, missing flag and every probability
within 1e-6 of our stored units / (2^bits - 1). It writes its metadata cache
to <file>.bgen.metadata2.mmm. Exits 1 and prints the first difference. One
process checks every file, since importing bgen-reader costs more than a
check.
"""

import sys

import numpy as np
from bgen_reader import open_bgen
from check_bgen_phased import fail, read_bgen, unpack


def check(path):
    ids, variants = read_bgen(path)
    with open_bgen(path, metadata_filepath=path + ".metadata2.mmm", allow_complex=True, verbose=False) as bgen:
        if list(bgen.samples) != ids:
            fail(f"{path}: bgen-reader samples {list(bgen.samples)}, want {ids}")
        if bgen.nvariants != len(variants):
            fail(f"{path}: bgen-reader has {bgen.nvariants} variants, want {len(variants)}")
        probs, missings, ploidies = bgen.read(return_missings=True, return_ploidies=True)
        for v, (varid, rsid, chrom, vpos, alleles, block) in enumerate(variants):
            where = f"{path} {chrom}:{vpos}"
            got = (bgen.ids[v], bgen.rsids[v], bgen.chromosomes[v], int(bgen.positions[v]), bgen.allele_ids[v])
            if got != (varid, rsid, chrom, vpos, ",".join(alleles)):
                fail(f"{where}: bgen-reader variant {got}")
            phased, top, samples = unpack(block, len(ids), len(alleles))
            if bool(bgen.phased[v]) != bool(phased):
                fail(f"{where}: bgen-reader phased={bgen.phased[v]}, want {phased}")
            want_ploidies = [ploidy for ploidy, _, _ in samples]
            want_missings = [missing for _, missing, _ in samples]
            bad = np.flatnonzero((ploidies[:, v] != want_ploidies) | (missings[:, v] != want_missings))
            if bad.size:
                s = bad[0]
                fail(f"{where} sample {s}: bgen-reader ploidy {ploidies[s, v]} missing {missings[s, v]}")
            # A sample's probabilities, then NaN to the row's end; a missing sample is all NaN.
            want = np.full((len(ids), probs.shape[2]), np.nan)
            for s, (_, missing, groups) in enumerate(samples):
                if not missing:
                    row = [u / top for g in groups for u in g]
                    want[s, : len(row)] = row
            bad = np.flatnonzero(~np.isclose(probs[:, v], want, rtol=0, atol=1e-6, equal_nan=True).all(axis=1))
            if bad.size:
                s = bad[0]
                fail(f"{where} sample {s}: bgen-reader {probs[s, v]}, want {'missing' if samples[s][1] else want[s]}")
    print(f"{path}: {len(variants)} variants")


def main():
    for path in sys.argv[1:]:
        check(path)


if __name__ == "__main__":
    main()
