# Byte identity with Beagle 5.5

fast-beagle writes the same VCF text as the Java release of Beagle 5.5, `beagle.27Feb25.75f.jar`, when both run with the same arguments, seed and `nthreads=`. The gate proves this on 35 recorded cases at 1, 2 and 18 threads, and compares 20 intermediate trace seams with an instrumented copy of the Java source. It also runs 200 generated inputs and 84 invalid-parameter inputs through both tools. It checks the Java library functions that Beagle calls against a real JVM, and runs every check on macOS arm64 and on Linux x86_64. The proof covers the tested inputs, parameters, thread counts and platforms. The [limits](#limits) section lists what it does not cover.

## What byte-identical means

The compared output is the decompressed VCF text. `tests/check-oracle.sh` decompresses `<out>.vcf.gz` and removes the `##source` and `##filedate` header lines. It then takes the first 16 hex digits of the SHA-256 of the rest. fast-beagle writes the same `##source="beagle.27Feb25.75f.jar"` line as the jar. `##filedate` holds the date of the run, so it differs between runs on different days.

The compressed `.vcf.gz` bytes are not compared, and they differ. On case `gt` at 2 threads, fast-beagle wrote a 37,177-byte file and the jar wrote a 39,810-byte file, with the same decompressed text.

Every comparison runs both tools at the same `nthreads=`. Beagle's own output depends on the thread count when a window is longer than 4 cM, because its PBWT phaser splits such windows by thread count. fast-beagle uses the same partitions, so it matches the jar run with the same thread count. On such windows, neither tool's output at one thread count matches its output at another.

## Evidence

Each check below runs in the pre-merge gate, `tests/gate-steps.sh`, except the chr20 benchmark.

| Evidence | What it checks | Scale | Rerun with |
|---|---|---|---|
| [Oracle hashes](#oracle-hashes) | The VCF hash of each case against the hash that the jar writes | 35 cases, each at 1, 2 and 18 threads, for the jar, a build of the Java source, the Java trace build and fast-beagle | `tests/check-oracle.sh build/beagle` |
| [Trace seams](#trace-seams) | Intermediate values at 20 points in the pipeline, against an instrumented Java build | 42 cases at 2 threads, and 9 seams again at 1 and 18 threads on the 2 thread-dependent cases | `tests/check-trace.sh T1a T1b ... T5d` |
| [Differential fuzzing](#differential-fuzzing) | Generated inputs and options run through the jar and fast-beagle | 200 examples at 1, 2, 3 or 5 threads, plus 2 examples for each of 42 invalid-parameter changes, plus 1 saved regression | `uv run --python 3.12 --script tests/check_fuzz.py --examples 200` |
| [Refused inputs](#refused-inputs) | Inputs that Beagle rejects, with the exit code and Java's message | 19 cases for fast-beagle, 18 of them also on the jar, and 49 fast-beagle-only output collision checks | `tests/check-failures.sh build/beagle` |
| [Java library fixtures](#java-library-fixtures) | C reproductions of the Java library behaviour that Beagle depends on, against a real JVM | 8 fixture sets | `make check-jcompat` |
| [Thread safety and ordering](#thread-safety-and-ordering) | Memory errors, undefined behaviour and data races on the oracle cases, the ordered writer protocol, and the imputation work-item size | Sanitizers on 35 cases, TLC on 8 model sizes, piece size on 20 cases at 2 and 18 threads | `tests/check-sanitizers.sh`, `tests/check-tsan.sh`, `tests/check-tla.sh`, `make check-piece-size` |
| [Platforms](#platforms) | The whole gate on two operating systems and two CPU architectures | macOS arm64 and Linux x86_64 | `tests/check-local.sh` |
| [chr20 benchmark](#chr20-benchmark) | A realistic imputation run, outside the gate | 1000 Genomes chr20, 6 runs per tool at 18 threads | `tests/bench/bench.sh <dir> 3` |

### Oracle hashes

`tests/oracle-cases.txt` has 35 cases. All but one run on Beagle's public test data, `test.beagle.vcf.gz`: 1,356 markers and 191 samples in 100 kb of chromosome 22. `tests/fetch-fixtures.sh` downloads it with the jar, checks both against SHA-256 checksums, and derives the fixtures. The cases cover phasing only, imputation, a VCF and a bref3 reference, chrX with haploid samples, multiallelic markers, missing genotypes, exclusion lists, a `chrom=` interval, genetic maps, `ap=` and `gp=`, and the state and segment limits. The other case, `imp-nan`, imputes with `err=0` from a 3-marker, 2-sample reference and writes `AF=NaN`. [Fixtures](testing.md#fixtures) lists them.

`tests/check-oracle.sh` runs each case at 1, 2 and 18 threads, 105 runs per implementation, with `seed=-99999`. Two cases, `gt-ibs2` and `gt-ibs2-miss`, use a steep genetic map with windows longer than 4 cM. They record one hash per thread count. The other 33 cases record one hash for all three thread counts.

The gate runs `tests/check-oracle.sh` four times, once for each of these implementations:

1. the release jar, `java -ea -jar data/beagle.27Feb25.75f.jar`
2. the unmodified Java source in `java/src/`, compiled by the gate
3. the Java source with `java/trace.patch` applied, with tracing on
4. fast-beagle, `build/beagle`

The first run proves on every gate run that the recorded hashes are the jar's output. The third run proves that the trace hooks do not change Java's output.

### Trace seams

A matching final hash can hide two errors that cancel. The trace seams compare the pipeline's intermediate state. `java/trace.patch` adds hooks to a copy of the Java source, and `make java-trace` builds it in `build/java-trace/`. With `-Dbeagle.trace=<dir>`, the Java build writes each seam to `<dir>/<seam>.txt`. fast-beagle writes the same files with `trace=<dir>`.

The 20 seams are `T1a T1b T1c T1d T2 T2b T3a T3b0 T3b1 T3b T3c T3d T4a T4b T4c T4d T5a T5b T5c T5d`. They cover input records and samples (T1), windows and fixed phasing data (T2), stage-1 phasing (T3), stage-2 phasing (T4) and imputation (T5). Seams that hold floating-point values print them as raw bits. Seams over large state print an FNV-1a digest of it.

`tests/check-trace.sh` runs both builds at 2 threads on the 35 oracle cases and the 7 input edge cases in `tests/trace-cases.txt`, 42 cases in all. The 7 extra cases cover CRLF line ends, Latin-1 bytes, a missing final newline, 600 markers 10 bp apart, a bref3 edge case, and two runs that both tools refuse. The script compares each seam file with `cmp`. It also requires the same exit code per case and fails if fast-beagle writes a seam file that Java does not. The gate then runs the 9 thread-dependent seams, `T3b0 T3b1 T3b T3c T3d T4a T4b T4c T4d`, again at 1 and 18 threads on `gt-ibs2` and `gt-ibs2-miss`.

### Differential fuzzing

`tests/check_fuzz.py` uses Hypothesis to generate small inputs. Each example has a target VCF of 1 to 400 markers and 1 to 30 samples, on chromosome `20`, `chr7` or `X`. Some examples add a reference panel of 2 to 40 samples and a genetic map. The generator varies missing genotypes, multiallelic markers, phased input, haploid chrX samples and the density of target markers. It also varies the options `seed`, `nthreads` (1, 2, 3 or 5), `burnin`, `iterations`, `phase-states`, `window`, `overlap`, `imp-states`, `impute`, `cluster`, `gp` and `ap`. Some fast-beagle runs add `bgen=phased` or `bgen=plink2`, and the VCF must still match.

For each example, both tools must succeed with the same VCF hash, or both must fail. A failure counts as agreement when Java prints an exception or `ERROR` line and fast-beagle exits 1 with a message. The two messages can differ, and the summary counts such examples. A fast-beagle crash never counts as agreement. When Java repeats a window that cannot advance without end, fast-beagle must exit 1 with its "does not advance" error.

The invalid-parameter examples change one parameter of a generated input so that Beagle rejects it. There are 42 changes: `out` naming an input or a directory, 15 `window` and `overlap` pairs, 21 values out of bounds or not numbers, and an unknown parameter. Both tools must exit 1 with the same message, and the input files must stay unchanged. Java may prefix the message with the exception class.

The full gate runs a fixed set of 200 examples and 2 examples per invalid-parameter change, and the nightly CI run fuzzes 1000 new examples on each runner. Each directory in `tests/fuzz-regressions/` runs first. It holds one input today, `window-stall`.

### Refused inputs

`tests/check-failures.sh` runs inputs that Beagle rejects at 2 threads. Each run must exit 1 and print the same message that the jar prints, and must write no output. The 19 cases are:

- 6 parameter errors: `out=` equal to the `gt=` file (three times: as given, with a doubled slash in `out=`, and with a doubled slash in `gt=`), `out=` equal to the `ref=` file, `out=` naming a directory, and `window` less than 1.1 times `overlap`. The `out=` messages must name the input by its normalized path, as `java.io.File` prints it.
- a `window` shorter than the marker spacing, which leaves a window with no marker
- 8 genetic map errors
- a one-character non-ASCII GT allele
- a bref3 SNV allele code with a negative permutation index
- `imp-segment` below half of `imp-step`, which divides by zero in Java
- a bref3 header whose sample count overflows when doubled

The gate runs the script on the jar and on fast-beagle. The bref3 sample-count case runs for fast-beagle only, because the jar's result depends on its heap size.

On fast-beagle only, the script also runs 49 checks of the [output collision refusal](beagle-divergences.md#errors-and-exit-status), which Beagle does not have. In 48 runs an output is an existing input file: through `./`, `..`, a relative path, a symlink, a hard link and a directory symlink, and each output file against each input file parameter. Each run must print the `fast-beagle:` message, leave the input unchanged and write no file. The 49th run names an input like an output that is not enabled, and must succeed.

### Java library fixtures

Beagle's output depends on the exact behaviour of Java library code. `src/jcompat/` reproduces that code in C. `make check-jcompat` compiles `tests/jcompat/JcompatFixtures.java`, runs it on the JVM, and requires each C fixture to print the same text. The 8 fixture sets are:

- `random`: `java.util.Random` sequences for 10 seeds and 14 bounds
- `math`: `StrictMath.log`, `log10`, `pow` and `expm1` from the fdlibm 5.3 port, on special values and 20,000 random inputs, and each Beagle call site of `Math.log`, `Math.log10`, `Math.pow` and `Math.expm1` evaluated as Beagle evaluates it
- `numbers`: `Math.round` edges, `double` and `float` casts to `int` and `long`, `DecimalFormat` patterns including NaN and the infinities, `Double.toString` and `Float.toString`
- `parse`: `Double.parseDouble`
- `parseint`: `Integer.parseInt`, including digits in other scripts
- `utf8`: Java's UTF-8 decoding of byte strings
- `pqueue`: `java.util.PriorityQueue` order with tied keys
- `search`: `Arrays.binarySearch` on arrays with runs of equal elements

### Thread safety and ordering

`tests/check-sanitizers.sh` builds fast-beagle and the unit tests with AddressSanitizer and UndefinedBehaviorSanitizer at `-O1`. It runs every oracle case at 1 and 2 threads, as VCF only, with `bgen=plink2` and with `bgen=phased`. Each run must exit 0 with the oracle hash and no sanitizer report. The `bgen=plink2` runs skip the 3 cases on non-autosomes. On Linux, LeakSanitizer also runs. On macOS, `tests/check-tsan.sh` builds fast-beagle with ThreadSanitizer and runs every oracle case at 18 threads, and the 2 thread-dependent cases again with tracing on.

`tests/check-tla.sh` model-checks [tla/ParallelOrdered.tla](../tla/ParallelOrdered.tla), the protocol of `parallel_ordered` in `src/blbutil/parallel.c`, which the phased and the imputed writer print records through. Workers build records in parallel, and the calling thread writes them in order. TLC checks 8 model sizes, from 1 worker and 3 items to 4 workers and 8 items, with windows of 1 to 5. It checks that items are consumed in order, that no slot is overwritten before it is consumed, that there is no deadlock, and that every run consumes every item. The model includes spurious wakeups, so a lost wakeup fails the check.

`make check-piece-size` builds a second binary that splits long imputation clusters into work items of one marker, not 500. It runs both binaries on the 20 oracle cases with `ref=` at 2 and 18 threads. Both must write the same VCF and the same trace seam T5d. The work-item size is a tuning value, and the check proves it does not change the output.

### Platforms

`tests/check-local.sh` runs the lint hooks, then every check in `tests/gate-steps.sh` natively on macOS arm64 and again on Linux x86_64 in docker. The Linux leg uses the `eclipse-temurin:21-jdk` image with the distribution's gcc. The GitHub Actions workflow `.github/workflows/gate.yml` runs the same list on `ubuntu-latest` for pull requests, and on `ubuntu-latest` and `macos-latest` for pushes to `main`. Every oracle hash is the same on both platforms, because both run the same case table.

The release workflow runs `tests/check-oracle.sh` on the binaries that users install: the conda package and the static binary, each on Linux x86_64, Linux arm64 and macOS arm64. [Releases](release.md) describes it.

### chr20 benchmark

`tests/bench/bench.sh` runs the jar and fast-beagle on a 1000 Genomes chr20 imputation. The target has 321 samples at 15,879 markers, and the reference has 2,881 samples at 1,642,181 markers. The script records the output hash as `tests/check-oracle.sh` computes it. In the recorded run, all 6 Java runs and all 6 fast-beagle runs at 18 threads wrote hash `38d4ed96cec67a44`. [perf-baseline.md](perf-baseline.md) has the method, the machine and the timings.

## Why identical output is possible

Floating-point results depend on the order of operations. fast-beagle keeps Java's order and Java's library behaviour.

- The Makefile appends `-ffp-contract=off -fno-fast-math -fwrapv` after any user `CFLAGS`, so every build uses them. `-ffp-contract=off` stops the compiler from fusing a multiply and an add, which Java never does. `-fno-fast-math` stops it from reordering floating-point operations. `-fwrapv` makes signed integer overflow wrap, as Java's `int` and `long` do.
- The C source has no `fma` calls and no SIMD intrinsics. Its only pragma is `#pragma omp simd`, on loops that have no sum, so vectorising them changes no value.
- Beagle calls `Math.log`, `Math.log10`, `Math.pow` and `Math.expm1`. fast-beagle calls the fdlibm 5.3 functions that define `StrictMath`, from `src/jcompat/fdlibm/`. At each Beagle call site the result is cast to `float` or floored, and the `math` fixture checks that `Math` and the fdlibm port then give the same value.
- `src/jcompat/jrandom.c` reproduces `java.util.Random`. `src/jcompat/jnum.c` reproduces Java's narrowing casts, `Math.round`, number parsing and number formatting.
- `nthreads=` sets both the thread count and the partitions that Beagle's output depends on, as in Beagle.
- Threaded loops write results in a fixed order. The phased and the imputed writer build records on worker threads and write them in item order.

## Limits

- The proof covers the inputs, parameters and thread counts that the checks run. The oracle cases use 100 kb of chromosome 22 with at most 5,068 reference samples. Only the chr20 benchmark runs at realistic scale, and it is outside the gate.
- Thread counts 1, 2 and 18 are checked on the oracle cases, and 1, 2, 3 and 5 by the fuzzer. Only 2 oracle cases have windows longer than 4 cM, where the thread count changes the output.
- The benchmark proves identity for one input, at 18 threads, on one macOS arm64 machine.
- The fuzzer accepts a failure in both tools with different messages. Only the invalid-parameter examples and `tests/check-failures.sh` require the same message.
- The `math` fixture checks `Math` against the fdlibm port at sampled inputs for each call site, not at every input.
- The gate runs on macOS arm64 with Apple clang and on Linux x86_64 with gcc. Linux arm64, macOS x86_64, Windows and other compilers are untested.
- The gate runs the jar on Java 21 only.
- The proof covers `beagle.27Feb25.75f.jar` only. Any other Beagle version is out of scope.
- The hash is the first 64 bits of a SHA-256. The compressed `.vcf.gz` bytes differ from the jar's.
- Console output is compared only for error messages.
- The TLA+ model is a hand-written description of `src/blbutil/parallel.c`. TLC checks the model at small sizes, not the C code.

BGEN and tabix outputs have no Java counterpart, so they have their own oracles. `bgen=plink2` files must match, byte for byte, what the pinned plink2 build writes from fast-beagle's VCF. `bgen=phased` files are checked against the VCF and decoded with bgen-reader ([BGEN checks](testing.md#bgen-checks)). `tbi=true` indexes must match htslib's `tbx_index_build3` byte for byte ([tabix index checks](testing.md#tabix-index-checks)).

## Reproduce the evidence

You need a C11 compiler, `make`, htslib with `bgzip` and `tabix`, Java 21 (`java` and `javac`), `curl`, `perl`, Python 3 and [uv](https://docs.astral.sh/uv/). The BGEN check also needs the [pinned plink2 build](testing.md#pinned-plink2-build).

Run these commands from the repository root. To fetch the fixtures and the jar and verify their checksums:

```bash
tests/fetch-fixtures.sh
```

To build fast-beagle and the Java trace build:

```bash
make build/beagle java-trace
```

To run the checks one at a time:

```bash
tests/check-oracle.sh java -ea -jar data/beagle.27Feb25.75f.jar
tests/check-oracle.sh build/beagle
tests/check-trace.sh T1a T1b T1c T1d T2 T2b T3a T3b0 T3b1 T3b T3c T3d T4a T4b T4c T4d T5a T5b T5c T5d
tests/check-failures.sh java -ea -jar data/beagle.27Feb25.75f.jar
tests/check-failures.sh build/beagle
uv run --python 3.12 --script tests/check_fuzz.py --examples 200
make check-jcompat check-piece-size
tests/check-sanitizers.sh
tests/check-tsan.sh        # macOS only
tests/check-tla.sh
```

To run the whole gate on one platform, set `PLINK2` to the pinned plink2 binary:

```bash
PLINK2=/path/to/plink2 tests/gate-steps.sh "$PWD"
```

To run the gate on both platforms, you also need [prek](https://github.com/j178/prek) and a docker that runs `linux/amd64` images. Set `PLINK2_ARM64` and `PLINK2_LINUX` to the macOS arm64 and Linux x86_64 plink2 binaries:

```bash
PLINK2_ARM64=/path/to/mac/plink2 PLINK2_LINUX=/path/to/linux/plink2 tests/check-local.sh
```

`tests/check-local.sh` prints one pass or fail line per check. Each native check writes its log to `build/check-<name>.log`. To try new fuzz examples, run `uv run --python 3.12 --script tests/check_fuzz.py --examples 1000 --random`. To rerun the chr20 benchmark, see [perf-baseline.md](perf-baseline.md).
