#!/bin/bash
# Builds build/beagle and the BGEN unit tests under AddressSanitizer and
# UndefinedBehaviorSanitizer, in build/san so the normal build is untouched,
# then runs every oracle case as VCF only, with bgen=plink2 and with
# bgen=phased. Every run must exit 0 with the oracle hash and no sanitizer
# report. On Linux LeakSanitizer also runs; it does not support macOS arm64.
# tests/check-tsan.sh runs ThreadSanitizer. Output refusal checks also run
# with leak detection off: util_exit leaves memory allocated by design.
#
# Usage: tests/check-sanitizers.sh
# NTHREADS overrides the thread counts tried (default "1 2").
# CASES restricts the run to the named cases (default all).
# SAN_JOBS sets how many runs go at once (default the CPU count): the runs at
# 1 and 2 threads leave cores idle otherwise.
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
# shellcheck source=cases.sh
source "$ROOT/tests/cases.sh"
SAN="$ROOT/build/san"

mkdir -p "$SAN"
rsync -a --delete --exclude build "$ROOT/Makefile" "$ROOT/java" "$ROOT/src" "$ROOT/third_party" "$ROOT/tests" "$SAN/"
SAN_CFLAGS=${SAN_CFLAGS:-"-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined,float-cast-overflow -fno-sanitize-recover=all"}
# clang cannot vectorise a `#pragma omp simd` loop that holds
# UndefinedBehaviorSanitizer checks, and reports that as -Wpass-failed, which
# -Werror turns into an error. The normal build keeps that error.
SAN_CFLAGS="$SAN_CFLAGS -Wno-pass-failed"
# make does not rebuild objects when only CFLAGS change.
make -C "$SAN" clean > /dev/null
make -C "$SAN" CFLAGS="$SAN_CFLAGS" build/beagle check-bgen-unit check-records check-tracker check-markers check-block-reader check-sliding-window check-hmm-updater build/output/bgen_files_fixture \
  > "$SAN/make.log" 2>&1 || { echo "FAIL sanitizer build or BGEN unit test ($SAN/make.log)"; exit 1; }
echo "PASS sanitizer build and BGEN unit tests"

leaks=0
[ "$(uname -s)" = Linux ] && leaks=1
# Reports go to the run's stderr, which is kept in its log.
export ASAN_OPTIONS="detect_leaks=$leaks"
export UBSAN_OPTIONS="print_stacktrace=1"
REPORT='runtime error:|ERROR: (Address|Leak|UndefinedBehavior)Sanitizer'
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

check_selection "$ROOT/tests/oracle-cases.txt" || exit 1
fail=0
if ASAN_OPTIONS=detect_leaks=0 python3 "$ROOT/tests/check_bgen_files.py" "$SAN/build/output/bgen_files_fixture"; then
  echo "PASS sanitizer BGEN partial-file cleanup"
else
  echo "FAIL sanitizer BGEN partial-file cleanup"; fail=1
fi
if ASAN_OPTIONS=detect_leaks=0 python3 "$ROOT/tests/check_output_failures.py" "$SAN/build/beagle"; then
  echo "PASS sanitizer output refusals"
else
  echo "FAIL sanitizer output refusals"; fail=1
fi
# The runs, one per line: name expect tags nthreads mode args...
runs() {
  while read -r name expect tags args; do
    selected "$name" || continue
    for t in ${NTHREADS:-1 2}; do
      for mode in vcf plink2 phased; do
        [ "$mode" = plink2 ] && has_tag "$tags" nonautosome && continue
        [ "$mode" = phased ] && has_tag "$tags" nonfinite && continue
        echo "$name $expect $tags $t $mode $args"
      done
    done
  done < <(cases)
}
# One run: prints its PASS or FAIL lines and returns 1 on FAIL. xargs runs it
# in a child bash, so it and what it calls are exported below.
# shellcheck disable=SC2329  # runs through xargs
run() {  # name expect tags nthreads mode args...
  local name=$1 expect=$2 t=$4 mode=$5 out bgen=() ok; shift 5
  out="$OUT/$name.t$t.$mode"
  [ "$mode" = vcf ] || bgen=(bgen="$mode")
  case_verdict "$expect" "$*" "$out" "$t" "$SAN/build/beagle" ${bgen[@]+"${bgen[@]}"}
  ok=$?
  if grep -Eq "$REPORT" "$out.run.log"; then
    echo "FAIL $name nthreads=$t $mode: sanitizer report"; grep -E -A12 "$REPORT" "$out.run.log" | head -40; return 1
  elif [ "$ok" -eq 0 ]; then
    echo "PASS $name nthreads=$t $mode"
  else
    echo "FAIL $name nthreads=$t $mode $VERDICT"; tail -5 "$out.run.log"; return 1
  fi
}
export -f run case_verdict case_run vcf_hash want_hash thread_dependent
export OUT SAN REPORT DATA SEED SHA
runs > "$OUT/runs"
# Each run's lines go to its own file, printed in table order once every run ends.
# shellcheck disable=SC2016  # the child bash expands the variables
xargs -P "${SAN_JOBS:-$(getconf _NPROCESSORS_ONLN)}" -L1 \
  bash -c 'run "$@" > "$OUT/$1.t$4.$5.lines" 2>&1 || touch "$OUT/$1.t$4.$5.fail"' _ < "$OUT/runs"
while read -r name _ _ t mode _; do
  cat "$OUT/$name.t$t.$mode.lines"
  [ ! -e "$OUT/$name.t$t.$mode.fail" ] || fail=1
done < "$OUT/runs"

[ $leaks -eq 1 ] && echo "leak checking on" || echo "leak checking off (not supported on $(uname -s) $(uname -m))"
exit $fail
