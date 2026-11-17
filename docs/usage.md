# Parameters and output files

fast-beagle takes Beagle 5.5's `key=value` arguments and writes Beagle's output. This page lists what differs from Beagle and the parameters that fast-beagle adds. For every other parameter, see the [Beagle 5.5 documentation](https://faculty.washington.edu/browning/beagle/beagle.html). [How fast-beagle differs from Beagle 5.5](beagle-divergences.md) lists every other difference.

## Run a phasing or imputation job

```bash
build/beagle gt=target.vcf.gz ref=ref.vcf.gz map=plink.map out=result
```

The run writes `result.vcf.gz`. Without `ref=` it phases the target only. The VCF has Beagle's DS, DR2, AF and optional AP1/AP2 and GP fields.

## Thread count

`nthreads=` sets both the thread count and the partitions that Beagle's output depends on. The output therefore matches Beagle run with the same thread count. Beagle's output depends on the thread count when a window is longer than 4 cM.

## Reference panels in bref3 format

`ref=` also reads bref3 files (`.bref3`). fast-beagle writes the output that Beagle writes from the same bref3 file. In Beagle and in fast-beagle alike, that output can differ from the output with a VCF reference of the same panel. It differs in two ways:

- First, bref3 stores the sequence coding computed over every sample and marker in the file. Beagle codes a VCF reference after `excludesamples=` and `excludemarkers=` remove records and samples. With either parameter the coding groups differ. Imputation clusters break at group boundaries, so imputed values change (cases `imp-bref3-excl` and `imp-bref3-excls`).
- Second, bref3 keeps INFO/END only for markers that are not SNVs, so imputed SNV records lose `END=` (case `imp-bref3-end`).

## BGEN output

fast-beagle has two BGEN modes. `bgen=plink2` writes the file plink2 would write from the VCF. `bgen=phased` writes each haplotype's probabilities before the VCF rounds them. Both modes write `<out>.bgen` and `<out>.sample` next to the VCF, and also [`<out>.info`](#the-info-file).

### Parameters

| Parameter | Values | Needs | Effect |
|---|---|---|---|
| `bgen=` | `plink2` or `phased` | | Writes the BGEN files. |
| `bgen-bits=N` | 1 to 16, default 8 | `bgen=` | Sets the bits per stored probability, in both BGEN modes. |
| `bgen-min-dr2=x` | a number | `bgen=plink2` | Adds `--extract-if-info "DR2 >= x"`. |
| `bgen-min-maf=x` | a number | `bgen=plink2` | Adds `--maf x`. |
| `bgen-chr-set=N` | 1 to 95 | `bgen=plink2` | Adds `--chr-set N`, for genomes with N autosomes. |

A `bgen-bits` value outside 1 to 16, or `bgen-bits` without `bgen=`, fails with a message.

### `bgen=plink2`

`bgen=plink2` writes `<out>.bgen` and `<out>.sample`. These files are byte-identical to the output of `plink2 --vcf <out>.vcf.gz dosage=DS --export bgen-1.2 bits=N`.

The writer starts from the values the VCF prints (rounded DS, printed DR2, GT phase). It ports plink2's import, filter and export rules, its number parser and its `.sample` number format. As a result, `bgen-min-dr2=0.3` drops a printed `DR2=0.30`, as plink2 does. The VCF is unchanged. Three cases differ from a plain VCF run:

- The BGEN leaves out multiallelic records, as `--import-max-alleles 2` does. plink2 without that flag aborts on them. The run reports how many records it skipped.
- `bgen=plink2` supports only autosomes, 1-22 or 1-N with `bgen-chr-set=N`. plink2 needs sex information for chrX (chromosome `23` by default, `N+1` with `bgen-chr-set=N`), chrY and MT. It needs `--allow-extra-chr` for other contigs. `bgen=plink2` exits at the first such window.
- If the filters remove every variant, the run exits with an error and writes no `.bgen`, as plink2 does.

PLINK 2 reserves sample ID `0`. A `bgen=plink2` run refuses that ID among the samples selected for export.

plink2 compresses BGEN with libdeflate at level 6, and zlib's own deflate gives different bytes. fast-beagle therefore compresses with the libdeflate 1.25 files in `third_party/libdeflate/`.

### `bgen=phased`

`bgen=phased` writes `<out>.bgen` and `<out>.sample` with each haplotype's allele probabilities before the VCF rounds them. The file is BGEN v1.2 layout 2, phased, with `bgen-bits` bits per probability, and zlib-compressed. It holds every record the VCF has, including multiallelic records and every chromosome, under the VCF's chromosome name and ID. Alleles are in VCF order with REF first, so plink2 reads the file with `--bgen <out>.bgen ref-first`. Haploid samples have ploidy 1. Genotyped markers and phasing-only runs store the GT alleles as certain. At a reference-only site (`ALT=.`), the BGEN allele labels are REF and `.`, with probability one on REF for each haplotype.

The writer scales each haplotype's probabilities to 2^N - 1 (255 at the default 8 bits). It rounds every value down, then adds one to the values with the largest remainders. On ties, lower alleles come first. The BGEN v1.2 specification suggests this method. The filters apply to `bgen=plink2` only.

### The info file

Both `bgen=` modes also write `<out>.info`, the fields BGEN has no place for. The file is tab-separated, with the header `CHROM POS ID REF ALT DR2 AF IMP`. It has one row per variant in the `.bgen`, in the same order.

- The first five columns are the VCF record's.
- `DR2` and `AF` are the values the VCF prints, or `.` where it has none (phasing-only runs).
- `IMP` is 1 for imputed markers and 0 otherwise.
- `CHROM` is the VCF's name even where `bgen=plink2` renames the chromosome in the `.bgen`.

### Failed runs

If a `bgen=` run fails before the BGEN files are complete, it removes the `.bgen`, `.info` and `.sample` it has created. `bgen=phased` refuses nonfinite allele probabilities, which can arise with `err=0`, because BGEN cannot represent them. VCF-only output retains Beagle's `NaN` values.

Before writing, fast-beagle refuses any log, VCF, BGEN, sidecar or tabix destination that names an input file, including aliases through symbolic or hard links. It also refuses two outputs that resolve to the same regular file, even if a symlink initially has a missing target. Each file opened by `trace=` is checked against the inputs and the enabled output destinations before truncation.

## Log file and console output

fast-beagle prints Beagle's progress report to standard output and writes the same text to `<out>.log`. The report has Beagle's lines: the start time, the command line with `nthreads=` added when you did not set it, the sample counts, the markers in each window, the estimated `ne` and `err`, the time of each phasing iteration and imputation step, and the totals. `tests/check-log.sh` checks that these lines match Beagle's.

The report differs from Beagle's in these lines:

| Line | fast-beagle |
|---|---|
| Banner | `fast-beagle: a C port of beagle.27Feb25.75f.jar (version 5.5)` and the copyright line |
| `Command line:` | The program path as you ran it, in place of `java -Xmx<heap>m -jar beagle.27Feb25.75f.jar` |
| `Total time:` | Includes closing the output files |
| `CPU time:` | User plus system CPU time. Beagle does not print it. |
| `Max memory:` | The max resident memory in MB. Beagle does not print it. |
| Last line | `fast-beagle finished` |

An error message goes to standard error without Java's exception class or stack trace. When the run fails after it has created `<out>.log`, fast-beagle also writes the message at the end of the log, so the log shows why the run stopped.

## Tabix index

`tbi=true` (default `false`) also writes `<out>.vcf.gz.tbi`, byte-identical to the index `tabix -p vcf <out>.vcf.gz` builds, without a second pass over the VCF. Replacing the VCF invalidates any existing tabix index at that prefix, including on a run with `tbi=false`. [Architecture](architecture.md#tabix-index-in-the-same-pass) describes how `src/main/vcf_index.c` builds the index.
