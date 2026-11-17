#!/usr/bin/env python3
"""Check output/input collisions and numeric output boundaries."""

import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BEAGLE = Path(sys.argv.pop(1)).resolve() if len(sys.argv) > 1 else ROOT / "build/beagle"
PROGRAM = "fast-beagle"
OPTIONAL_SUFFIXES = (".bgen", ".info", ".sample", ".vcf.gz.tbi")
SUFFIXES = (".vcf.gz", ".log", *OPTIONAL_SUFFIXES)
VCF = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
    "1\t100\tv1\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\t0|0\n"
    "1\t200\tv2\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\t0|0\n"
)
INPUTS = {
    "gt": VCF,
    "ref": VCF,
    "map": "1 v1 0 100\n1 v2 1 200\n",
    "excludesamples": "S3\n",
    "excludemarkers": "absent-marker\n",
    "truth": VCF,
    "ped": "",
}


def run(args):
    return subprocess.run(
        [str(BEAGLE), "nthreads=2", "seed=-99999", *[f"{k}={v}" for k, v in args.items()]],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


class OutputFailures(unittest.TestCase):
    def collision(self, field, suffix, mode, alias="direct"):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            prefix = root / "result"
            dest = Path(str(prefix) + suffix)
            contents = INPUTS[field].encode()
            if suffix == ".vcf.gz":
                contents = gzip.compress(contents)
            source = dest
            if alias in ("symlink", "hardlink"):
                source = root / "source"
            source.write_bytes(contents)
            if alias == "symlink":
                dest.symlink_to(source)
            elif alias == "hardlink":
                os.link(source, dest)
            elif alias == "relative":
                (root / "sub").mkdir()
                prefix = root / "sub" / ".." / "result"
            args = {"gt": target, "out": prefix, field: source, "tbi": "true"}
            if mode is not None:
                args["bgen"] = mode
            proc = run(args)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            # Beagle refuses only a VCF path that equals the gt= or ref= path.
            if suffix == ".vcf.gz" and field in ("gt", "ref") and alias == "direct":
                self.assertIn(f"ERROR: VCF output file equals input file: {source}", proc.stderr)
            else:
                self.assertIn(f"{PROGRAM}: output file {prefix}{suffix} equals input file ", proc.stderr)
            self.assertEqual(source.read_bytes(), contents, "input was overwritten")
            if alias == "symlink":
                self.assertTrue(dest.is_symlink(), "input link was removed")
            for ext in SUFFIXES:
                if ext != suffix:
                    self.assertFalse(Path(str(prefix) + ext).exists(), f"unexpected {ext}")

    def test_input_collisions(self):
        for mode in ("phased", "plink2"):
            for field in INPUTS:
                for suffix in SUFFIXES:
                    with self.subTest(mode=mode, field=field, suffix=suffix):
                        self.collision(field, suffix, mode)

    def test_alias_collisions(self):
        for alias in ("relative", "symlink", "hardlink"):
            for suffix in SUFFIXES:
                with self.subTest(alias=alias, suffix=suffix):
                    self.collision("gt", suffix, "phased", alias)

    def test_output_alias_collisions(self):
        pairs = [(a, b) for j, a in enumerate(SUFFIXES) for b in SUFFIXES[j + 1 :]]
        for alias in ("symlink", "hardlink"):
            for first, second in pairs:
                with self.subTest(alias=alias, first=first, second=second), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    target = root / "target.vcf"
                    target.write_text(VCF)
                    prefix = root / "result"
                    kept, linked = Path(str(prefix) + first), Path(str(prefix) + second)
                    kept.write_text("untouched\n")
                    if alias == "symlink":
                        linked.symlink_to(kept)
                    else:
                        os.link(kept, linked)
                    proc = run({"gt": target, "out": prefix, "bgen": "phased", "tbi": "true"})
                    self.assertEqual(proc.returncode, 1, proc.stderr)
                    self.assertEqual(
                        proc.stderr.strip(),
                        f"{PROGRAM}: output file {prefix}{first} equals output file {prefix}{second}",
                    )
                    self.assertEqual(kept.read_text(), "untouched\n", "output was overwritten")
                    for ext in SUFFIXES:
                        if ext not in (first, second):
                            self.assertFalse(Path(str(prefix) + ext).exists(), f"unexpected {ext}")

    def test_dangling_output_aliases(self):
        for suffix in (".vcf.gz", ".info", ".sample"):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                target = root / "input.vcf"
                target.write_text(VCF)
                prefix = root / "result"
                shared = root / "shared"
                (root / "result.log").symlink_to("shared")
                Path(str(prefix) + suffix).symlink_to("result.log")
                args = {"gt": target, "out": prefix, "bgen": "phased"}
                proc = run(args)
                self.assertEqual(proc.returncode, 1, proc.stderr)
                self.assertIn("result.log", proc.stderr)
                self.assertIn("result" + suffix, proc.stderr)
                self.assertTrue((root / "result.log").is_symlink())
                self.assertFalse(shared.exists() and shared.stat().st_size > 0)
                Path(str(prefix) + suffix).unlink()
                proc = run(args)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertIn("finished", shared.read_text())

    def test_trace_preserves_inputs(self):
        for kind in ("direct", "hardlink", "symlink"):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                seam = root / "T2.txt"
                target = seam if kind == "direct" else root / "target.vcf"
                target.write_text(VCF)
                if kind == "hardlink":
                    os.link(target, seam)
                elif kind == "symlink":
                    seam.symlink_to(target)
                proc = run({"gt": target, "out": root / "result", "trace": root})
                self.assertEqual(proc.returncode, 1, proc.stderr)
                self.assertIn(str(seam), proc.stderr)
                self.assertEqual(target.read_text(), VCF)

    def test_vcf_replacement_invalidates_index(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target, index = root / "target.vcf", root / "result.vcf.gz.tbi"
            target.write_text(VCF)
            args = {"gt": target, "out": root / "result"}
            proc = run({**args, "tbi": "true"})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            old_index = index.read_bytes()
            target.write_text(VCF.replace("\n1\t", "\n2\t"))
            # Even when index generation is disabled, a declared input is protected.
            proc = run({**args, "truth": index})
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertEqual(index.read_bytes(), old_index)
            proc = run(args)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertFalse(index.exists())
            with gzip.open(root / "result.vcf.gz", "rt") as f:
                self.assertEqual([row.split("\t")[0] for row in f if not row.startswith("#")], ["2", "2"])

    def test_reserved_plink_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF.replace("S1", "0"))
            args = {"gt": target, "out": root / "result"}
            proc = run({**args, "bgen": "plink2"})
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn("sample ID cannot be '0'", proc.stderr)
            for suffix in (".bgen", ".sample", ".info"):
                self.assertFalse((root / ("result" + suffix)).exists())
            exclude = root / "exclude.txt"
            exclude.write_text("0\n")
            proc = run({**args, "bgen": "plink2", "excludesamples": exclude})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            proc = run({**args, "bgen": "phased"})
            self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_chromosome_whitespace(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            for char in ("\u1680", "\u2003", "\u2028", "\u3000", "\u00a0", "\u2007", "\u202f", "\U0001f9ec"):
                with self.subTest(char=hex(ord(char))):
                    chrom = f"chr{char}1"
                    target.write_text(VCF.replace("\n1\t", f"\n{chrom}\t"))
                    proc = run({"gt": target, "out": root / "result"})
                    if char in "\u1680\u2003\u2028\u3000":
                        self.assertEqual(proc.returncode, 1, proc.stderr)
                        self.assertIn("CHROM field contains whitespace", proc.stderr)
                    else:
                        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_outputs_may_share_a_device(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            prefix = root / "result"
            for suffix in (".info", ".sample"):
                Path(str(prefix) + suffix).symlink_to("/dev/null")
            proc = run({"gt": target, "out": prefix, "bgen": "phased"})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(gzip.decompress(Path(str(prefix) + ".vcf.gz").read_bytes()).count(b"\n1\t"), 2)

    def test_log_collision_without_bgen(self):
        self.collision("gt", ".log", None)

    def test_disabled_outputs_do_not_collide(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            for suffix in OPTIONAL_SUFFIXES:
                with self.subTest(suffix=suffix):
                    prefix = root / suffix[1:]
                    exclude = Path(str(prefix) + suffix)
                    exclude.write_text("S3\n")
                    proc = run({"gt": target, "out": prefix, "excludesamples": exclude})
                    self.assertEqual(proc.returncode, 0, proc.stderr)
                    self.assertEqual(exclude.read_text(), "S3\n")

    def test_large_population_size(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\n"
                "1\t20\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\n"
            )
            proc = run({"gt": target, "out": root / "result", "ne": "1e20", "burnin": 1, "iterations": 1})
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertIn("9223372036854775807", (root / "result.log").read_text())
            self.assertNotIn("runtime error:", proc.stderr)
            self.assertNotIn("Sanitizer", proc.stderr)

    def test_preflight_error_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "result"
            target = Path(str(prefix) + ".vcf.gz")
            target.write_bytes(gzip.compress(VCF.encode()))
            args = {"gt": target, "out": prefix, "window": 1, "overlap": 2}
            proc = run(args)
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertEqual(proc.stderr.strip(), f"ERROR: VCF output file equals input file: {target}")
            proc = run({**args, "unknown": 1})
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertEqual(proc.stderr.strip(), "Error: unrecognized parameter: unknown=1")
            self.assertEqual(gzip.decompress(target.read_bytes()).decode(), VCF)
            self.assertFalse(Path(str(prefix) + ".log").exists())

    def test_partial_output_lifecycle(self):
        for mode in ("phased", "plink2"):
            for blocked in (".info", ".vcf.gz", ".sample"):
                with self.subTest(mode=mode, blocked=blocked), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    target = root / "target.vcf"
                    target.write_text(VCF)
                    prefix = root / "result"
                    Path(str(prefix) + blocked).mkdir()
                    sample = Path(str(prefix) + ".sample")
                    index = Path(str(prefix) + ".vcf.gz.tbi")
                    if blocked != ".sample":
                        sample.write_text("untouched sample\n")
                    index.write_text("untouched index\n")
                    proc = run({"gt": target, "out": prefix, "bgen": mode, "tbi": "true"})
                    self.assertEqual(proc.returncode, 1, proc.stderr)
                    self.assertEqual(proc.stderr.strip(), f"Error opening {prefix}{blocked}")
                    self.assertEqual(Path(str(prefix) + ".log").read_text().splitlines()[-1], proc.stderr.strip())
                    self.assertFalse(Path(str(prefix) + ".bgen").exists())
                    if blocked != ".info":
                        self.assertFalse(Path(str(prefix) + ".info").exists())
                    if blocked == ".sample":
                        self.assertTrue(sample.is_dir())
                        self.assertTrue(Path(str(prefix) + ".vcf.gz").is_file())
                        self.assertNotEqual(index.read_bytes(), b"untouched index\n")
                    else:
                        self.assertEqual(sample.read_text(), "untouched sample\n")
                        self.assertEqual(index.read_text(), "untouched index\n")

    def test_read_ahead_bgen_error_order(self):
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"

        def rec(chrom, pos, last="1|0"):
            return f"{chrom}\t{pos}\t.\tA\tC\t.\t.\t.\tGT\t0|1\t{last}\n"

        for source in ("target", "reference"):
            for mode in ("phased", "plink2"):
                with self.subTest(source=source, mode=mode), tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    target = root / "target.vcf"
                    target.write_text(
                        header + rec(1, 100) + rec(2, 100) + rec(2, 200, "1|x" if source == "target" else "1|0")
                    )
                    prefix = root / "result"
                    sample = Path(str(prefix) + ".sample")
                    sample.write_text("untouched sample\n")
                    trace = root / "trace"
                    trace.mkdir()
                    args = {"gt": target, "out": prefix, "bgen": mode, "tbi": "true", "trace": trace}
                    if source == "reference":
                        reference = root / "reference.vcf"
                        reference.write_text(
                            header
                            + rec(1, 100)
                            + "".join(rec(2, pos) for pos in range(100, 100000, 100))
                            + "".join(rec(3, pos, "1|x" if pos == 5000 else "1|0") for pos in range(100, 10100, 100))
                        )
                        args["ref"] = reference
                    proc = run(args)
                    error = "java.lang.IllegalArgumentException: Window has only one position: CHROM=1 POS=100"
                    self.assertEqual(proc.returncode, 1, proc.stderr)
                    self.assertEqual(proc.stderr.strip(), error)
                    self.assertEqual(Path(str(prefix) + ".log").read_text().splitlines()[-1], error)
                    self.assertFalse(Path(str(prefix) + ".bgen").exists())
                    self.assertFalse(Path(str(prefix) + ".info").exists())
                    self.assertEqual(sample.read_text(), "untouched sample\n")

    def test_nonfinite_phased_probabilities(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp) / "result"
            proc = run(
                {
                    "gt": ROOT / "data/target.thin.vcf.gz",
                    "ref": ROOT / "data/ref.vcf.gz",
                    "out": prefix,
                    "err": 0,
                    "bgen": "phased",
                }
            )
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertIn(f"{PROGRAM}: bgen=phased: cannot encode a non-finite allele probability", proc.stderr)
            self.assertNotIn("runtime error:", proc.stderr)
            self.assertNotIn("Sanitizer", proc.stderr)
            for suffix in (".bgen", ".info", ".sample"):
                self.assertFalse(Path(str(prefix) + suffix).exists(), f"partial {suffix}")

    @unittest.skipUnless(Path("/dev/full").exists(), "requires Linux /dev/full")
    def test_completed_bgen_survives_log_close_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "target.vcf"
            target.write_text(VCF)
            prefix = root / "result"
            Path(str(prefix) + ".log").symlink_to("/dev/full")
            proc = run({"gt": target, "out": prefix, "bgen": "phased", "tbi": "true"})
            self.assertEqual(proc.returncode, 1, proc.stderr)
            self.assertEqual(proc.stderr.strip(), f"Error writing {prefix}.log")
            for suffix in (".bgen", ".info", ".sample", ".vcf.gz", ".vcf.gz.tbi"):
                self.assertGreater(Path(str(prefix) + suffix).stat().st_size, 0, suffix)


if __name__ == "__main__":
    unittest.main()
