#!/usr/bin/env python3
"""Check completed records and captured probabilities through the writer."""

import gzip
import subprocess
import tempfile
from pathlib import Path

from check_bgen_phased import decode, read_bgen

ROOT = Path(__file__).resolve().parent.parent
# Wraps the 4-thread ring of 8 tasks of 50 records twice, with a partial last task.
PHASED_MARKERS = 1013


def main():
    with tempfile.TemporaryDirectory() as tmp:
        outputs = []
        for mode in range(3):
            for nthreads in (1, 4):
                out = Path(tmp) / f"{mode}.t{nthreads}"
                subprocess.run([ROOT / "build/output/record_fixture", out, str(mode), str(nthreads)], check=True)
                with gzip.open(f"{out}.vcf.gz", "rt") as f:
                    outputs.append([s.rstrip("\n").split("\t") for s in f if s[0] != "#"])
        assert all(o == outputs[0] for o in outputs[1:])
        phased, imputed, observed, mono, five = outputs[0][:5]
        assert phased[6:] == [".", ".", "GT", "0|1", "1", "0|1"]
        assert imputed[6] == "PASS" and imputed[7].endswith(";IMP")
        assert imputed[8] == "GT:DS:AP1:AP2:GP"
        assert imputed[9].split(":")[:4] == ["0|1", "1.21", "0.33", "0.88"]
        assert imputed[10:] == ["0:0.25:0.25", "0|0:0:0:0:1,0,0"]
        assert observed[7].endswith(";END=301")
        assert observed[9:] == [
            "2|0:0,1:0,1:0,0:0,0,0,1,0,0",
            "1:1,0:1,0",
            "0|0:0,0:0,0:0,0:1,0,0,0,0,0",
        ]
        assert mono[7] == "IMP" and mono[9:] == ["0|0:1", "0", "0|0:1"]
        assert five[9].startswith("0|0:0,0,0,0:0,0,0,0:0,0,0,0:1,")
        print("PASS phased, imputed, observed, haploid and multiallelic records")

        # The phased driver's rows, in marker order, each from the fixture's
        # allele function: haplotype hap's allele is bit hap % 4 of the marker.
        rows = outputs[0][5:]
        assert len(rows) == PHASED_MARKERS
        for m, row in enumerate(rows):
            a = [(m >> (hap % 4)) & 1 for hap in range(6)]
            assert row[1] == str(1000 + m), (m, row[:2])
            assert row[8:] == ["GT", f"{a[0]}|{a[1]}", f"{a[2]}", f"{a[4]}|{a[5]}"], (m, row[8:])
        print(f"PASS {PHASED_MARKERS} phased rows in marker order at 1 and 4 threads")

        ids, variants = read_bgen(Path(tmp) / "2.t1.bgen")
        assert ids == ["D1", "H", "D2"] and len(variants) == 5 + PHASED_MARKERS
        top, probs = decode(variants[1][-1], 3, 2)
        assert top == 65535
        # At 16 bits, rounding to the VCF's two decimals before capture fails.
        for got, want in zip(probs[0], [0.3338, 0.8766], strict=True):
            assert abs(got[1] / top - want) < 1 / top
        assert abs(probs[1][0][1] / top - 0.25) < 1 / top
        assert probs[2] == [[top, 0], [top, 0]]
        _, probs = decode(variants[2][-1], 3, 3)
        assert probs == [
            [[0, 0, top], [top, 0, 0]],
            [[0, top, 0]],
            [[top, 0, 0], [top, 0, 0]],
        ]
        print("PASS BGEN captures probabilities before buffer reuse and VCF rounding")


if __name__ == "__main__":
    main()
