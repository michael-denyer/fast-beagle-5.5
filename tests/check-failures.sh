#!/bin/bash
# Runs a Beagle implementation on arguments that Beagle refuses. Each run must
# exit 1 and print Java's message, and a refused run must leave its input
# files unchanged and write no output.
#
# Usage: tests/check-failures.sh <command...>
#   tests/check-failures.sh java -ea -jar data/beagle.27Feb25.75f.jar
#   tests/check-failures.sh build/beagle
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

fail=0
failed() { echo "FAIL $*"; fail=1; }

# ABSENT names an output the refused run must not write.
check() {  # name message args...
  local name=$1 message=$2; shift 2
  "${BEAGLE[@]}" "$@" seed=$SEED nthreads=2 > "$OUT/$name.log" 2>&1
  refused "$OUT/$name.log" $? "$message" ${ABSENT:+"$ABSENT"} || { failed "$name $VERDICT"; return; }
  echo "PASS $name: $message"
}

unchanged() {  # name file hash
  [ "$($SHA "$2" | cut -c1-16)" = "$3" ] || failed "$1 changed $(basename "$2")"
}

BEAGLE=("$@")
cp "$DATA/target.vcf.gz" "$OUT/in.vcf.gz"
cp "$DATA/ref.vcf.gz" "$OUT/ref.vcf.gz"
in_hash=$($SHA "$OUT/in.vcf.gz" | cut -c1-16)
ref_hash=$($SHA "$OUT/ref.vcf.gz" | cut -c1-16)
mkdir "$OUT/dir"

# Main.parameters and Main.checkOutputPrefix. File.equals compares normalized
# paths, so a doubled slash still names the input, and the message prints the
# input's normalized path.
norm=$(tr -s / <<< "$OUT")
check out-equals-gt "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT/in.vcf.gz" out="$OUT/in"
unchanged out-equals-gt "$OUT/in.vcf.gz" "$in_hash"
check out-equals-gt-slashes "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT/in.vcf.gz" out="$OUT//in"
unchanged out-equals-gt-slashes "$OUT/in.vcf.gz" "$in_hash"
check gt-slashes "ERROR: VCF output file equals input file: $norm/in.vcf.gz" gt="$OUT//in.vcf.gz" out="$OUT/in"
unchanged gt-slashes "$OUT/in.vcf.gz" "$in_hash"
check out-equals-ref "ERROR: VCF output file equals input file: $norm/ref.vcf.gz" \
  ref="$OUT/ref.vcf.gz" gt="$DATA/target.thin.vcf.gz" out="$OUT/ref"
unchanged out-equals-ref "$OUT/ref.vcf.gz" "$ref_hash"
ABSENT="$OUT/dir.vcf.gz" check out-directory "ERROR: \"out\" parameter cannot be a directory" \
  gt="$DATA/target.vcf.gz" out="$OUT/dir"
ABSENT="$OUT/window.vcf.gz" check window-overlap \
  "ERROR: The \"window\" parameter must be at least 1.1 times the \"overlap\" parameter" \
  gt="$DATA/target.vcf.gz" out="$OUT/window" window=1 overlap=1
# A window shorter than the marker spacing holds no marker, and BasicGT
# indexes the empty marker array.
check empty-window "java.lang.ArrayIndexOutOfBoundsException: Index 0 out of bounds for length 0" \
  gt="$DATA/target.thin.vcf.gz" window=0.0000001 overlap=0.00000001 out="$OUT/empty-window"

# VcfRecGTParser reads a one-character allele as c - '0' and longer ones with
# Integer.parseInt, so a lone Arabic-Indic one (U+0661) is allele 1585.
gzip -dc "$DATA/target.vcf.gz" | LC_ALL=C awk 'BEGIN {OFS = "\t"} /^#/ {print; next}
  !done {$10 = sprintf("%c%c", 217, 161) substr($10, 2); done = 1} {print}' | gzip > "$OUT/one-char.vcf.gz"
check one-char-allele "ERROR: Invalid allele [$(printf '\331\241')]" gt="$OUT/one-char.vcf.gz" out="$OUT/one-char-out"

# VcfIt reads a block of 1024 lines before it parses them, so when the stream
# ends early inside the block that holds a bad line, the read error comes first.
awk 'BEGIN {
  srand(7)
  print "##fileformat=VCFv4.2"
  h = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
  for (s = 1; s <= 100; s++) h = h "\tS" s
  print h
  for (r = 1; r <= 900; r++) {
    line = "1\t" (r * 100) "\t.\tA\tC\t.\t.\t.\tGT"
    for (s = 1; s <= 100; s++) line = line "\t" (r == 2 && s == 1 ? "1|x" : int(rand() * 2) "|" int(rand() * 2))
    print line
  }
}' | gzip > "$OUT/whole.vcf.gz"
head -c $(($(wc -c < "$OUT/whole.vcf.gz") * 3 / 4)) "$OUT/whole.vcf.gz" > "$OUT/in-bad-cut.vcf.gz"
check bad-line-then-cut "Error reading $OUT/in-bad-cut.vcf.gz" gt="$OUT/in-bad-cut.vcf.gz" out="$OUT/bad-cut"

# Window 1 is a lone marker, which MarkerMap rejects, and window 2 holds a bad
# allele, in the target or in the reference's second block of 1024 lines. The
# port reads window 2 while window 1 is phased but reports window 1's error.
# Beagle names the bad allele: it parses the target's first 1024 lines at
# startup, and the reference block when its reader thread reaches it, which
# can be before window 1 fails.
if [ "$1" != java ]; then
  vcf_head() { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n'; }
  rec() { printf '%s\t%s\t.\tA\tC\t.\t.\t.\tGT\t0|1\t%s\n' "$1" "$2" "${3:-1|0}"; }
  lone="java.lang.IllegalArgumentException: Window has only one position: CHROM=1 POS=100"
  { vcf_head; rec 1 100; rec 2 100; rec 2 200 '1|x'; } | gzip > "$OUT/in-next-bad.vcf.gz"
  check next-window-target "$lone" gt="$OUT/in-next-bad.vcf.gz" out="$OUT/next-target"
  { vcf_head; rec 1 100; rec 2 100; rec 2 200; } | gzip > "$OUT/in-next-targ.vcf.gz"
  { vcf_head; rec 1 100
    for p in $(seq 100 100 99900); do rec 2 "$p"; done
    for p in $(seq 100 100 10000); do if [ "$p" = 5000 ]; then rec 3 "$p" '1|x'; else rec 3 "$p"; fi; done
  } | gzip > "$OUT/in-next-ref.vcf.gz"
  check next-window-ref "$lone" ref="$OUT/in-next-ref.vcf.gz" gt="$OUT/in-next-targ.vcf.gz" out="$OUT/next-ref"

  # Marker validation used to call exit directly, bypassing the reader's
  # util_try. Reject any extra next-window diagnostic as well as a wrong exit.
  first_window_error() {  # name args...
    local name=$1; shift
    check "$name" "$lone" "$@"
    if grep -Eq '^Duplicate marker:|^markers not in chromosomal order:' "$OUT/$name.log"; then
      failed "$name reported the next window's marker error"
    fi
  }
  { vcf_head; rec 1 100; rec 2 100; rec 2 200; rec 2 200; } | gzip > "$OUT/in-next-duplicate.vcf.gz"
  { vcf_head; rec 1 100; rec 2 100; rec 2 300; rec 2 200; } | gzip > "$OUT/in-next-order.vcf.gz"
  for bad in duplicate order; do
    first_window_error "next-window-$bad-target" \
      gt="$OUT/in-next-$bad.vcf.gz" out="$OUT/next-$bad-target"
    first_window_error "next-window-$bad-ref" \
      ref="$OUT/in-next-$bad.vcf.gz" gt="$OUT/in-next-targ.vcf.gz" out="$OUT/next-$bad-ref"
  done
fi

# PlinkGenMap throws IllegalArgumentException, and prints a genetic position
# with Double.toString.
map_check() {  # name message map-lines...
  local name=$1 message=$2; shift 2
  printf '%s\n' "$@" > "$OUT/$name.map"
  check "$name" "java.lang.IllegalArgumentException: $message" gt="$DATA/target.vcf.gz" map="$OUT/$name.map" out="$OUT/$name"
}
map_check map-same-cm "All loci in genetic map have the same genetic position [0.480887]: 22 . 0.480887 20000000" \
  "22 . 0.480887 20000000" "22 . 0.480887 20100000"
map_check map-same-cm-sci "All loci in genetic map have the same genetic position [1.0E-4]: 22 . 0.0001 20000000" \
  "22 . 0.0001 20000000" "22 . 0.00010 20100000"
map_check map-infinite "invalid map position: Infinity" "22 . 0.1 20000000" "22 . +Infinity 20100000"
map_check map-duplicate "duplication position: 22 . 0.2 20000000" "22 . 0.1 20000000" "22 . 0.2 20000000"
map_check map-descending "map positions not in ascending order: 22 . 0.05 20100000" "22 . 0.1 20000000" "22 . 0.05 20100000"
map_check map-short-line "Map file format error: 22 . 0.1" "22 . 0.1"
map_check map-long-line "Map file format error: 22 . 0.1 20000000 x" "22 . 0.1 20000000 x" "22 . 0.2 20100000"
map_check map-missing-chrom "missing genetic map for chromosome 22" "21 . 0.1 20000000" "21 . 0.2 20100000"

# Bref3Reader.readMarker indexes SNV_PERMS with alleleCode >> 2, so allele code
# byte 0x80 on the first marker (20000086, rs138720731) is permutation -32.
LC_ALL=C perl -0777 -pe 's/(\x01\x31\x2d\x56\x01\x00\x0brs138720731)./$1\x80/s or die "no first marker\n"' \
  "$DATA/ref.bref3" > "$OUT/snv-code.bref3" || failed "bref3-snv-code cannot edit $DATA/ref.bref3"
check bref3-snv-code "java.lang.ArrayIndexOutOfBoundsException" \
  ref="$OUT/snv-code.bref3" gt="$DATA/target.thin.vcf.gz" out="$OUT/snv-code"

# A bref3 header with 2^30 samples, truncated after four empty IDs; the sample
# count doubled overflows int. Java first allocates String[2^30], so it reads
# to EOF with a heap of 4 GB or more and runs out of memory with less. Only the
# port's result, EOF, is fixed.
if [ "$1" != java ]; then
  perl -e 'print pack("N", 2055763188), pack("n", 1), "x", pack("N", 0x40000000), "\0\0" x 4' > "$OUT/samples.bref3"
  check bref3-sample-count "java.io.EOFException" \
    ref="$OUT/samples.bref3" gt="$DATA/target.thin.vcf.gz" out="$OUT/samples"
fi

# An imp-step too small to change a cluster's position when added to it.
# CodedSteps.stepStarts then never reaches the next cluster, so Java appends
# the same step start until it runs out of memory. Only the port's refusal is
# fixed.
if [ "$1" != java ]; then
  check imp-step-no-advance "fast-beagle: imp-step=1e-20 is too small to advance" \
    ref="$DATA/ref.vcf.gz" gt="$DATA/target.thin.vcf.gz" out="$OUT/step" imp-step=1e-20
fi

# ImpIbs divides imp-states by round(imp-segment / imp-step), which is 0 here.
check imp-segment "java.lang.ArithmeticException: / by zero" \
  ref="$DATA/ref.vcf.gz" gt="$DATA/target.thin.vcf.gz" out="$OUT/segment" imp-segment=0.01

exit $fail
