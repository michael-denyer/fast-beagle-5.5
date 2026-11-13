# How fast-beagle differs from Beagle 5.5

fast-beagle is a C port of Java Beagle 5.5 (`beagle.27Feb25.75f.jar`). Run with the same arguments and the same `nthreads=`, it writes a VCF whose text is byte-identical to Beagle's, header lines included. Everything else on this page is a difference you can see when you switch. It refuses `ped=` and runs without a JVM or a heap limit. It also adds BGEN output, a tabix index and a trace parameter.

## Differences at a glance

| Area | Java Beagle 5.5 | fast-beagle |
|---|---|---|
| VCF text | The reference output | Byte-identical at the same `nthreads=` |
| `.vcf.gz` file bytes | BGZF from Beagle's own writer | BGZF from htslib, so the compressed bytes differ |
| `<out>.log` and standard output | Beagle's progress report | The same report, naming fast-beagle, with `CPU time:` and `Max memory:` added and a failed run's error at the end of the log. See [Log file and console output](usage.md#log-file-and-console-output). |
| No arguments | Prints the usage text and exits 0 | Prints `missing gt argument` and exits 1 |
| Error output | The message, often with the exception class, a stack trace, the usage text or the banner | The message only |
| Exit status on error | 1 | 1 |
| `ped=` | Accepted and ignored | Refused |
| Added parameters | None | `bgen=`, `bgen-bits=`, `bgen-min-dr2=`, `bgen-min-maf=`, `bgen-chr-set=`, `tbi=`, `trace=` |
| Default `nthreads=` | The processor count the JVM reports | The processor count `sysconf` reports |
| Memory limit | The JVM heap size (`-Xmx`) | None |
| Requirements | A Java runtime | A native build linked to htslib |
| bref3 files | Read by `ref=`. A separate bref3 jar writes them. | Read by `ref=`. fast-beagle has no bref3 writer. |

## Output files

### The VCF

The decompressed VCF is byte-identical to Beagle's when both run with the same arguments and `nthreads=`. The header is the same too. fast-beagle writes `##source="beagle.27Feb25.75f.jar"`, and `##filedate=` holds the run date in both tools. A downstream tool therefore cannot tell from the header which program wrote the file.

The compressed `.vcf.gz` files differ byte for byte. Beagle compresses with its own BGZF writer and fast-beagle compresses with htslib. Both files are valid BGZF, and `bgzip -t` accepts both. Compare the decompressed text, not the `.vcf.gz` bytes.

### BGEN files and the tabix index

fast-beagle can also write `<out>.bgen`, `<out>.sample`, `<out>.info` and `<out>.vcf.gz.tbi`. Beagle writes none of these. [Parameters and output files](usage.md) describes them.

## No arguments

With no arguments, Beagle prints its usage text and exits with status 0. fast-beagle treats an empty command line like any other command line without `gt=`. It prints `missing gt argument` to standard error and exits with status 1.

## Errors and exit status

Both tools exit with status 1 when they refuse a run. fast-beagle prints the same message text as Beagle, on standard error. It leaves out the lines Beagle adds around that message:

- For a parameter error, Beagle prints the message after `Exception in thread "main" java.lang.IllegalArgumentException:` and follows it with a stack trace.
- For an input file that does not exist, Beagle adds `java.lang.Throwable: File does not exist` and a stack trace.
- For an unrecognized parameter, a missing input file and several input errors, Beagle ends with a blank line and `Terminating program.`
- For an `out=` that is a directory or names an input file, Beagle prints the full usage text before the message.
- For a `window=` less than 1.1 times `overlap=`, Beagle prints its banner before the message and `Exiting program.` after it.

Beagle refuses an `out=` whose VCF path names the `gt=` or `ref=` file after `java.io.File` normalizes the path text, and fast-beagle prints the same message. fast-beagle then also refuses any output that is an existing input file, even through `./`, `..`, a relative path, a symlink or a hard link, and prints `fast-beagle: output file <output> equals input file <input>`. Beagle overwrites the input in these cases. The check compares the VCF, log, and enabled BGEN and tabix outputs with every input file parameter before any output is opened.

fast-beagle also refuses two enabled outputs that are the same existing file through a symbolic or hard link, such as `<out>.log` linked to `<out>.vcf.gz`, and prints `fast-beagle: output file <output> equals output file <output>`. Beagle writes both outputs into the one file and exits 0.

An `imp-step=` so small that adding it to a marker cluster's genetic position leaves the position unchanged, such as `1e-20`, stops fast-beagle with `fast-beagle: imp-step=<value> is too small to advance from genetic position <position>`. Beagle's step loop never ends in this case and allocates until it runs out of memory.

When an input file has several malformed records, fast-beagle reports the first. At `nthreads=` above 1, Beagle's parse threads race, so the record Beagle names can change from run to run.

Both tools read the next window while they phase the current one. When the current window fails, for example because it has only one position, and the next window cannot be read, Beagle reports either error, depending on thread timing. fast-beagle always reports the current window's error. Both tools parse the reference in blocks of 1024 lines and report a malformed line when its block is read, which can be before an earlier window fails. Beagle also parses the first 1024 target lines at startup. fast-beagle reports a malformed target line only when reading the windows reaches it.

When a command line has more than one unrecognized parameter, both tools list them in one message. Beagle lists them in hash-map order. fast-beagle lists them in the order you gave them.

A few errors exist only in fast-beagle, such as a refused `ped=` or a failure to start a thread. Their messages start with `fast-beagle:`. So do its refusals of nonfinite probabilities in phased BGEN output and of log, VCF, BGEN and tabix output destinations that are input files, including symbolic and hard-link aliases. The one exception is a VCF output path equal to the `gt=` or `ref=` path, which both tools refuse with Beagle's message. [Failed runs](usage.md#failed-runs) describes these checks.

## Parameters

fast-beagle accepts every parameter that Beagle's parser accepts, with the same default, the same valid range and the same meaning. That includes the parameters Beagle's usage text leaves out (`initial-lr`, `step-scale`, `rare`, `imp-segment`, `imp-step`, `imp-nsteps`, `buffer` and `truth`). The one exception is `ped=`. Beagle checks that the file exists and then ignores it. fast-beagle stops with `fast-beagle: the ped= parameter is not supported`. Remove `ped=` from a Beagle command line before you run it with fast-beagle.

fast-beagle adds these parameters. Beagle refuses each of them as an unrecognized parameter.

| Parameter | Effect |
|---|---|
| `bgen=`, `bgen-bits=`, `bgen-min-dr2=`, `bgen-min-maf=`, `bgen-chr-set=` | Write BGEN v1.2 files next to the VCF. See [BGEN output](usage.md#bgen-output). |
| `tbi=true` | Writes `<out>.vcf.gz.tbi`. See [Tabix index](usage.md#tabix-index). |
| `trace=<dir>` | Writes internal trace files to an existing directory, for comparison with a traced Java build. See [Compare trace seams](testing.md#compare-trace-seams). |

## Threads

Beagle's output depends on `nthreads=` when a window is longer than 4 cM. fast-beagle reproduces that dependence, so its VCF matches Beagle's run with the same `nthreads=`, and a run at another thread count can differ from both. See [Thread count](usage.md#thread-count).

Without `nthreads=`, both tools use the processor count. Beagle reads it from the JVM and fast-beagle reads it from `sysconf`. To reproduce a Beagle run exactly, set `nthreads=` to the value in the Beagle log's command line.

## Runtime and memory

fast-beagle is a native program. It needs no Java runtime and no `-Xmx` setting. The build needs a C11 compiler, `make` and htslib, and the binary links to htslib at run time. See [Install](../README.md#install). The [pre-merge gate](testing.md#run-the-pre-merge-gate) builds and tests it on macOS arm64 and on Linux x86_64.

The JVM's maximum heap caps Beagle's memory. Without `-Xmx`, the JVM sets the cap itself. On a 64 GB machine it chose 16 GB, and Beagle's log showed the command line as `java -Xmx16384m`. fast-beagle has no such cap. [Performance](perf-baseline.md) compares the two tools' wall time, CPU time and max memory, from a 321-sample benchmark to an imputation of about 100,000 samples.

## Input files

fast-beagle opens `gt=` and `ref=` files with Beagle's rules. The rules cover gzip and BGZF detection, bref3 references and the reference file-name extensions. A bref3 reference can give different imputed values than a VCF of the same panel in both tools. See [Reference panels in bref3 format](usage.md#reference-panels-in-bref3-format).

fast-beagle has no bref3 writer. To convert a VCF reference to bref3, use Beagle's `bref3.27Feb25.75f.jar`.

Beagle stops reading a BGZF file at a run of empty BGZF blocks and drops the rest of the file without an error. Its `BGZipIt` inflates blocks in batches of 4 × `nthreads`, and a batch in which every block is empty ends the input. fast-beagle reads the whole file, as htslib does. A file made by `bgzip` or by concatenating `.vcf.gz` files has at most one empty block between members, which neither tool stops at. A file whose middle holds 8 empty blocks gives this difference at `nthreads=1`: Beagle writes the records before the run, and fast-beagle writes all of them. `make check-line-reader` tests that fast-beagle reads past such a run.
