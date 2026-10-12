#!/bin/bash
# bgen=plink2 must write the .bgen and .sample that plink2 writes from the
# same VCF: plink2 --vcf <out>.vcf.gz dosage=DS --export bgen-1.2 bits=N
# [--extract-if-info "DR2 >= x"] [--maf x] [--chr-set N], with N from
# bgen-bits (default 8). bgen=phased must write a BGEN that
# tests/check_bgen_phased.py finds consistent with the VCF and that plink2
# --bgen ref-first loads. In both modes tests/check_bgen_info.py must find the
# .info sidecar matching the BGEN and VCF, tests/check_bgen_reader.py must find
# that bgen-reader decodes every BGEN as our own decoder does (one process over
# every file after the runs, since importing bgen-reader costs more than a
# check), and a run that fails anywhere after the writer opens must leave no
# .bgen, .info or .sample.
# The VCF must keep its oracle hash. Invalid
# bgen-bits and bgen-chr-set values must exit 1 with a message, as must every
# refused run.
#
# Every run's .bgen, .sample and .info must also hash as tests/bgen-hashes.txt
# records. BGEN_ORACLE=recorded checks those hashes instead of plink2,
# bgen-reader and the Python checks, and needs neither plink2 nor uv. RECORD=1 rewrites tests/bgen-hashes.txt from a
# live run in which every other check passes.
#
# Usage: PLINK2=<plink2 v2.0.0-a.7.8> [RECORD=1] tests/check-bgen.sh
#        BGEN_ORACLE=recorded tests/check-bgen.sh
# The live checks need uv on PATH; it runs check_bgen_reader.py with its pinned
# bgen-reader.
# BEAGLE overrides the binary (default build/beagle). CASES restricts the run
# to the named cases.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
BEAGLE=${BEAGLE:-$ROOT/build/beagle}
PLINK2_VERSION="PLINK v2.0.0-a.7.8"

HASHES="$ROOT/tests/bgen-hashes.txt"
ORACLE=${BGEN_ORACLE:-live}
case $ORACLE in live|recorded) ;; *) echo "FAIL BGEN_ORACLE must be live or recorded, not $ORACLE"; exit 1 ;; esac
if [ -n "${RECORD:-}" ] && { [ "$ORACLE" != live ] || [ -n "${CASES:-}" ]; }; then
  echo "FAIL RECORD=1 needs the live oracle and every case"; exit 1
fi

if [ "$ORACLE" = live ]; then
  if [ -z "${PLINK2:-}" ] || [ ! -x "$PLINK2" ]; then
    echo "FAIL PLINK2 must name a plink2 binary (got '${PLINK2:-}')"; exit 1
  fi
  version=$("$PLINK2" --version 2>&1)
  if [[ "$version" != "$PLINK2_VERSION"* ]]; then
    echo "FAIL $PLINK2 is '$version', want $PLINK2_VERSION"; exit 1
  fi
  command -v uv > /dev/null || { echo "FAIL uv is not on PATH; it runs tests/check_bgen_reader.py (https://docs.astral.sh/uv/)"; exit 1; }
fi


TABLES=("$ROOT/tests/oracle-cases.txt" "$ROOT/tests/bgen-cases.txt")
# Every case runs with bgen=plink2, which must fail as plink2 does on a
# nonautosome case (plink2 needs sex information for chrX), and with
# bgen=phased. These runs add arguments: case, mode, beagle's arguments, then
# after | the outcome, ok or the message bgen=plink2 must fail with while
# plink2 fails on the VCF too.
RUNS='
imp         plink2  bgen-min-dr2=0.3 bgen-min-maf=0.05            | ok
imp         plink2  bgen-min-dr2=0.30                             | ok
imp-map     plink2  bgen-min-maf=0.2                              | ok
gt          plink2  bgen-min-maf=0.05                             | ok
# gt has no DR2, so its DR2 filter removes every variant.
gt          plink2  bgen-min-dr2=0.3                              | no variants remaining
# Under --chr-set 38, chromosome 23 in imp-chroms is an autosome, and
# chromosome 38 is chrX under --chr-set 37.
imp-chr38   plink2  bgen-chr-set=38                               | ok
imp-chroms  plink2  bgen-chr-set=38                               | ok
imp-chr38   plink2  bgen-chr-set=38 bgen-min-dr2=0.3 bgen-bits=3  | ok
imp-chr38   plink2  bgen-chr-set=37                               | supports autosomes 1-37 only
imp         plink2  bgen-bits=1                                   | ok
imp         plink2  bgen-bits=3                                   | ok
imp         plink2  bgen-bits=12                                  | ok
imp         plink2  bgen-bits=16                                  | ok
gt-phased   plink2  bgen-bits=1                                   | ok
gt-phased   plink2  bgen-bits=3                                   | ok
gt-phased   plink2  bgen-bits=12                                  | ok
gt-phased   plink2  bgen-bits=16                                  | ok
imp-hap22   plink2  bgen-bits=3                                   | ok
imp         plink2  bgen-min-dr2=0.3 bgen-min-maf=0.05 bgen-bits=12 | ok
imp         phased  bgen-bits=16                                  | ok
imp         phased  bgen-bits=3                                   | ok
imp-multi   phased  bgen-bits=1                                   | ok
imp-x       phased  bgen-bits=12                                  | ok
'
# Invalid bgen-bits: the arguments, then the message beagle must print.
BAD_ARGS=(
  "bgen=plink2 bgen-bits=0|Error in \"bgen-bits\" argument: value=0 < 1"
  "bgen=phased bgen-bits=17|Error in \"bgen-bits\" argument: value=17 > 16"
  "bgen=plink2 bgen-bits=abc|abc is not a number"
  "bgen-bits=8|bgen-bits needs bgen=plink2 or bgen=phased"
  "bgen=plink2 bgen-chr-set=0|Error in \"bgen-chr-set\" argument: value=0 < 1"
  "bgen=plink2 bgen-chr-set=96|Error in \"bgen-chr-set\" argument: value=96 > 95"
  "bgen=phased bgen-chr-set=38|bgen-chr-set needs bgen=plink2"
)
THREADS=2

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
fail=0
pass() { echo "PASS $*"; }
bgen_reader() { PYTHONDONTWRITEBYTECODE=1 uv run --python 3.12 --script "$ROOT/tests/check_bgen_reader.py" "$@" 2>&1; }
# The BGENs that passed every other live check, for one bgen_reader call.
READER_FILES=()
failed() { echo "FAIL $*"; fail=1; }

# Sets expect, tags and args to the case's row.
read_case() {  # name
  read -r expect tags args <<< "$(cases "${TABLES[@]}" | awk -v n="$1" '$1 == n {$1 = ""; print; exit}')"
}

# The run's row in tests/bgen-hashes.txt must match the SHA-256 prefixes of its
# .bgen, .sample and .info. Under RECORD=1 the row is collected instead. Sets
# VERDICT on a mismatch.
check_hashes() {  # out key
  local key=$2 f got=() want
  for f in "$1.bgen" "$1.sample" "$1.info"; do got+=("$($SHA "$f" | cut -c1-16)"); done
  if [ -n "${RECORD:-}" ]; then echo "${got[*]} $key" >> "$OUT/hashes"; return; fi
  want=$(awk -v k="$key" '!/^#/ && NF {h = $1 " " $2 " " $3; $1 = $2 = $3 = ""; sub(/^ +/, "")
    if ($0 == k) {print h; exit}}' "$HASHES")
  [ -n "$want" ] || { VERDICT="has no row in tests/bgen-hashes.txt"; return 1; }
  [ "${got[*]}" = "$want" ] || { VERDICT="hashes ${got[*]}, recorded $want"; return 1; }
}

run_plink2() {  # vcf out beagle-args...
  local vcf=$1 out=$2 bits=8 a; shift 2
  local flags=()
  has_tag "$tags" multiallelic && flags+=(--import-max-alleles 2)
  for a in "$@"; do
    case "$a" in
      bgen-min-dr2=*) flags+=(--extract-if-info "DR2 >= ${a#*=}") ;;
      bgen-min-maf=*) flags+=(--maf "${a#*=}") ;;
      bgen-bits=*) bits=${a#*=} ;;
      bgen-chr-set=*) flags+=(--chr-set "${a#*=}") ;;
    esac
  done
  "$PLINK2" --vcf "$vcf" dosage=DS --export bgen-1.2 bits="$bits" ${flags[@]+"${flags[@]}"} --out "$out" > "$out.plink2.log" 2>&1
}

# bgen=plink2 must write the .bgen and .sample that plink2 writes from the VCF.
# shellcheck disable=SC2329  # runs through check_ok as live_$mode
live_plink2() {  # out beagle-args...
  local b=$1 p=$1.plink2; shift
  run_plink2 "$b.vcf.gz" "$p" "$@" || { echo "plink2 exit=$?: $(grep -m1 '^Error' "$p.plink2.log")"; return 1; }
  if ! cmp "$b.bgen" "$p.bgen" || ! cmp "$b.sample" "$p.sample"; then
    echo "bgen or sample differs"; return 1
  fi
  echo "$(head -c 12 "$b.bgen" | tail -c 4 | od -An -tu4 | tr -d ' ') variants"
}

# bgen=phased: check_bgen_phased.py finds the BGEN consistent with the VCF,
# and plink2 loads the BGEN unless the case is nonautosome.
# shellcheck disable=SC2329  # runs through check_ok as live_$mode
live_phased() {  # out beagle-args...
  local b=$1 got
  got=$(python3 "$ROOT/tests/check_bgen_phased.py" "$b.bgen" "$b.vcf.gz") || { echo "$got"; return 1; }
  if ! has_tag "$tags" nonautosome; then
    "$PLINK2" --bgen "$b.bgen" ref-first --sample "$b.sample" --export vcf --out "$b.plink2" > "$b.plink2.log" 2>&1 \
      || { echo "plink2 --bgen: $(grep -m1 '^Error' "$b.plink2.log")"; return 1; }
  fi
  echo "$got"
}

# A run that must succeed. Live checks run only under BGEN_ORACLE=live.
check_ok() {  # label out mode key beagle-args...
  local label=$1 b=$2 mode=$3 key=$4 rows detail=; shift 4
  case_verdict "$expect" "$args" "$b" $THREADS "$BEAGLE" bgen="$mode" "$@" \
    || { failed "$label beagle $VERDICT $(tail -1 "$b.run.log")"; return; }
  if [ "$ORACLE" = live ]; then
    detail=$("live_$mode" "$b" "$@") || { failed "$label $detail"; return; }
    rows=$(python3 "$ROOT/tests/check_bgen_info.py" "$b.bgen" "$b.info" "$b.vcf.gz") || { failed "$label .info: $rows"; return; }
    READER_FILES+=("$b.bgen")
    detail="$detail, .info $rows, "
  fi
  check_hashes "$b" "$key" || { failed "$label $VERDICT"; return; }
  pass "$label ${detail}hashes as recorded"
}

# Both tools must fail: bgen=plink2 with a message naming the cause and no
# .bgen, .info or .sample left behind, plink2 on the VCF that beagle writes
# without bgen=.
check_both_fail() {  # label out message beagle-args...
  local label=$1 b=$2.beagle v=$2.vcfonly p=$2.plink2 message=$3; shift 3
  case_run "$args" "$b" $THREADS "$BEAGLE" bgen=plink2 "$@"
  refused "$b.run.log" $? "$message" "$b.bgen" "$b.info" "$b.sample" || { failed "$label bgen=plink2 $VERDICT"; return; }
  [ "$ORACLE" = live ] || { pass "$label fails: $message"; return; }
  case_verdict "$expect" "$args" "$v" $THREADS "$BEAGLE" || { failed "$label beagle without bgen= $VERDICT"; return; }
  if run_plink2 "$v.vcf.gz" "$p" "$@"; then failed "$label plink2 succeeded"; return; fi
  pass "$label fails in both: $(grep -m1 -A1 '^Error' "$p.plink2.log" | tr '\n' ' ')"
}

check_run() {  # name mode outcome beagle-args...
  local name=$1 mode=$2 outcome=$3 expect tags args out; shift 3
  read_case "$name"
  out="$OUT/$name$(printf '%s' "$*" | tr -c 'a-z0-9.' '_')"
  case "$mode $outcome" in
    "plink2 ok" | "phased ok") check_ok "$name bgen=$mode${*:+ $*}" "$out.$mode" "$mode" "$name $mode${*:+ $*}" "$@" ;;
    plink2\ *) check_both_fail "$name${*:+ $*}" "$out" "$outcome" "$@" ;;
    phased\ *)
      case_run "$args" "$out.phased" $THREADS "$BEAGLE" bgen=phased "$@"
      refused "$out.phased.run.log" $? "$outcome" "$out.phased.bgen" "$out.phased.info" "$out.phased.sample" \
        || { failed "$name bgen=phased $VERDICT"; return; }
      pass "$name bgen=phased fails: $outcome" ;;
    *) failed "$name${*:+ $*}: no check for mode $mode with outcome $outcome" ;;
  esac
}

all_runs() {
  local name tags
  while read -r name _ tags _; do
    if has_tag "$tags" nonautosome; then
      echo "$name plink2 | supports autosomes 1-22 only"
    else
      echo "$name plink2 | ok"
    fi
    if has_tag "$tags" nonfinite; then
      echo "$name phased | cannot encode a non-finite allele probability"
    else
      echo "$name phased | ok"
    fi
  done < <(cases "${TABLES[@]}")
  awk '!/^#/ && NF' <<< "$RUNS"
}

# Beagle must reject the arguments before reading any input.
check_bad_args() {  # "args|message"
  local args=${1%%|*} message=${1#*|} out="$OUT/bad-args"
  # shellcheck disable=SC2086  # args is an argument list; splitting it is the point
  "$BEAGLE" gt="$DATA/test.vcf.gz" out="$out" $args > "$out.stdout" 2> "$out.stderr"
  refused "$out.stderr" $? "$message" || { failed "$args $VERDICT"; return; }
  pass "$args fails: $message"
}

# A run that fails after the BGEN writer opens must remove what it wrote.
check_no_partial() {  # label out message beagle-args...
  local label=$1 out=$2 message=$3; shift 3
  "$BEAGLE" "$@" out="$out" > "$out.stdout" 2> "$out.stderr"
  refused "$out.stderr" $? "$message" "$out.bgen" "$out.info" "$out.sample" || { failed "$label $VERDICT"; return; }
  pass "$label fails and leaves no .bgen, .info or .sample"
}

partial_runs() {
  local stall="$ROOT/tests/fuzz-regressions/window-stall"
  mkdir "$OUT/open-fail.info" "$OUT/vcf-fail.vcf.gz"
  check_no_partial "bgen=phased, .info cannot be opened" "$OUT/open-fail" "Error opening $OUT/open-fail.info" \
    ref="$DATA/ref.vcf.gz" gt="$DATA/target.vcf.gz" bgen=phased
  [ -d "$OUT/open-fail.info" ] || failed "bgen=phased, .info cannot be opened removed a path it did not write"
  check_no_partial "bgen=plink2, VCF cannot be opened" "$OUT/vcf-fail" "Error opening $OUT/vcf-fail.vcf.gz" \
    ref="$DATA/ref.vcf.gz" gt="$DATA/target.vcf.gz" bgen=plink2
  # shellcheck disable=SC2046  # args.txt is an argument list; splitting it is the point
  check_no_partial "bgen=phased, fails after writing records" "$OUT/stall" "does not advance" \
    gt="$stall/targ.vcf" ref="$stall/ref.vcf" $(cat "$stall/args.txt") bgen=phased
}

check_selection "${TABLES[@]}" || fail=1

while IFS='|' read -r run outcome; do
  read -r name mode extra <<< "$run"
  read -r outcome <<< "$outcome"
  # shellcheck disable=SC2086  # extra is an argument list; splitting it is the point
  selected "$name" && check_run "$name" "$mode" "$outcome" $extra
done < <(all_runs)
if [ ${#READER_FILES[@]} -gt 0 ]; then
  if decoded=$(bgen_reader "${READER_FILES[@]}"); then
    pass "bgen-reader decodes all ${#READER_FILES[@]} BGENs as our decoder does"
  else
    failed "bgen-reader: $(tail -1 <<< "$decoded")"
  fi
fi
[ -n "${CASES:-}" ] || for a in "${BAD_ARGS[@]}"; do check_bad_args "$a"; done
[ -n "${CASES:-}" ] || partial_runs
if [ -n "${RECORD:-}" ] && [ "$fail" -eq 0 ]; then
  {
    echo "# The SHA-256 prefixes of the .bgen, .sample and .info that each run of"
    echo "# tests/check-bgen.sh writes, then the run: case, mode and added arguments."
    echo "# Written by RECORD=1 PLINK2=<plink2> tests/check-bgen.sh after plink2,"
    echo "# bgen-reader and our own decoder accept every file."
    LC_ALL=C sort -k4 "$OUT/hashes"
  } > "$HASHES"
  pass "recorded $(wc -l < "$OUT/hashes" | tr -d ' ') runs in tests/bgen-hashes.txt"
fi
exit $fail
