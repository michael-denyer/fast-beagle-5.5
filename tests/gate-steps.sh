#!/bin/bash
# The gate's checks, in order, run in the checkout <root>. tests/check-local.sh
# runs them natively and on Linux x86_64 in docker, and CI runs them on each
# runner. Prints one pass or FAIL line per check; each check's output goes to
# build/check-<name>.log, or build/check-c-<name>.log in the C tier.
#
# Usage: [GATE_TIER=c] [GATE_GROUP=<group>] [GATE_FUZZ=random] [GATE_LIST=1]
#   [GATE_LIST_GROUPS=1] tests/gate-steps.sh <root>
# Needs Java 21, htslib and uv; the full tier also needs PLINK2 naming the
# pinned plink2 binary (see tests/check-bgen.sh). GATE_TIER=c runs only the checks that compare the C
# binary against recorded results; Java then only builds the bref3 fixtures.
# It skips every check that runs Java or the jar alongside it, the
# fixture-cache check and the TLA+ model. It checks the
# BGEN output against the hashes in tests/bgen-hashes.txt instead of
# plink2, and runs the saved fuzz regressions but no new fuzz
# examples. The default is the full gate. GATE_FUZZ=random fuzzes 1000 new
# examples instead of the fixed 200.
#
# Each check names its group, after an optional tier (full or c) and an
# optional macos, which limits it to macOS and prints a skip line elsewhere.
# GATE_GROUP selects one declared group, so CI can run the groups as parallel
# jobs; the setup checks run in every group. The default, all, runs every check
# in order. GATE_LIST=1 prints the group and name of each check the tier and
# group select, without running it. GATE_LIST_GROUPS=1 prints the groups other
# than setup that have a check in the tier, on any OS, without creating files.
# shellcheck disable=SC2329  # java_build, oracle_trace and trace_threads run through step
set -uo pipefail
cd "$1" || exit 1
SEAMS="T1a T1b T1c T1d T2 T2b T3a T3b0 T3b1 T3b T3c T3d T4a T4b T4c T4d T5a T5b T5c T5d"
# Seams whose content depends on nthreads, rechecked at 1 and 18 threads on the
# cases whose windows are long enough for the thread count to split them: the
# cases with per-thread hashes.
THREAD_SEAMS="T3b0 T3b1 T3b T3c T3d T4a T4b T4c T4d"

fail=0
tier=${GATE_TIER:-full}
case $tier in full|c) ;; *) echo "GATE_TIER must be full or c, not $tier"; exit 2 ;; esac
group=${GATE_GROUP:-all}
bgen_oracle=live logs=build/check-
[ "$tier" = c ] && bgen_oracle=recorded logs=build/check-c-
fuzz_args=(--examples 200)
case ${GATE_FUZZ:-fixed} in
  fixed) ;;
  random) fuzz_args=(--examples 1000 --random) ;;
  *) echo "GATE_FUZZ must be fixed or random, not $GATE_FUZZ"; exit 2 ;;
esac
in_group() { [ "$group" = all ] || [ "$1" = setup ] || [ "$1" = "$group" ]; }
# mode is declare (record each check's tier and group), list or run.
step() {  # [full|c] [macos] group name command...
  local only_tier=any only_macos=
  case $1 in full|c) only_tier=$1; shift ;; esac
  if [ "$1" = macos ]; then only_macos=1; shift; fi
  local owner=$1 name=$2; shift 2
  if [ "$mode" = declare ]; then
    [ "$owner" = setup ] || declared+=("$only_tier $owner")
    return 0
  fi
  in_group "$owner" || return 0
  if [ "$only_tier" != any ] && [ "$only_tier" != "$tier" ]; then
    if [ "$mode" = run ] && [ "$only_tier" = full ]; then echo "  skip  $name (full tier only)"; fi
    return 0
  fi
  if [ -n "$only_macos" ] && [ "$(uname -s)" != Darwin ]; then
    if [ "$mode" = run ]; then echo "  skip  $name (macOS only)"; fi
    return 0
  fi
  if [ "$mode" = list ]; then echo "$owner $name"; return 0; fi
  if "$@" > "$logs$name.log" 2>&1; then
    echo "  pass  $name"
  else
    echo "  FAIL  $name ($logs$name.log)"; fail=1
  fi
}

java_build() {
  mkdir -p build/classes \
    && find java/src -name '*.java' > build/java-sources.txt \
    && javac -nowarn -d build/classes @build/java-sources.txt
}

oracle_trace() {
  mkdir -p build/trace-oracle \
    && tests/check-oracle.sh java -Dbeagle.trace=build/trace-oracle -cp build/java-trace/classes main.Main
}

trace_threads() {
  local t thread_cases
  # shellcheck source=cases.sh
  thread_cases=$(ROOT=$PWD; source tests/cases.sh; cases | while read -r name expect _; do
    if thread_dependent "$expect"; then printf '%s ' "$name"; fi; done) || return 1
  for t in 1 18; do
    # shellcheck disable=SC2086  # the seam list splits into arguments
    NTHREADS=$t CASES="$thread_cases" tests/check-trace.sh $THREAD_SEAMS || return 1
  done
}

checks() {
step setup fixtures tests/fetch-fixtures.sh --ensure
step core gate-tier tests/check-gate-tier.sh
step core log-recording python3 tests/check_log_recording.py
step core make-phase python3 tests/check_make_phase.py
step core lock-exit python3 tests/check_lock_exit.py
step full cases cases python3 tests/check_cases.py
step full core jcompat make check-jcompat
step core tracker make check-tracker
step core oom make check-oom
step core interval make check-interval
step core markers make check-markers
step core block-reader make check-block-reader
step core snv-perms make check-snv-perms
step full java oracle-jar tests/check-oracle.sh java -ea -jar data/beagle.27Feb25.75f.jar
step full java failures-jar tests/check-failures.sh java -ea -jar data/beagle.27Feb25.75f.jar
step full java log-jar tests/check-log.sh java -ea -jar data/beagle.27Feb25.75f.jar
step full java java-build java_build
step full java oracle-source tests/check-oracle.sh java -ea -cp build/classes main.Main
step full java java-trace make java-trace
step full java oracle-trace oracle_trace
step setup c-build make build/beagle
step core oracle-c tests/check-oracle.sh build/beagle
step core gate-planning tests/check-gate-planning.sh
step core failures-c tests/check-failures.sh build/beagle
step core output-failures python3 tests/check_output_failures.py build/beagle
step core log-c tests/check-log.sh build/beagle
step core piece-size make check-piece-size
step core bgen-unit make check-bgen-unit
step core records make check-records
step core bgen-files make check-bgen-files
step core vcf-index make check-vcf-index
step core tbi make check-tbi
step bgen bgen env BGEN_ORACLE="$bgen_oracle" tests/check-bgen.sh
# shellcheck disable=SC2086  # the seam list splits into arguments
step full java trace tests/check-trace.sh $SEAMS
step sanitizers sanitizers tests/check-sanitizers.sh
step macos tsan tsan tests/check-tsan.sh
step full core tla tests/check-tla.sh
step full core fuzz uv run --python 3.12 --script tests/check_fuzz.py "${fuzz_args[@]}"
step c core fuzz-regressions uv run --python 3.12 --script tests/check_fuzz.py --examples 0 --invalid-examples 0
step full java trace-threads trace_threads
}

mode=declare declared=()
checks
groups_in() {  # tier: each group with a check in that tier (any: in any tier), in declaration order
  printf '%s\n' "${declared[@]}" | awk -v t="$1" '(t == "any" || $1 == "any" || $1 == t) && !seen[$2]++ {print $2}'
}
if [ "$group" != all ] && [[ " $(groups_in any | tr '\n' ' ') " != *" $group "* ]]; then
  echo "GATE_GROUP must be all or a group declared in tests/gate-steps.sh, not $group"; exit 2
fi
if [ "${GATE_LIST_GROUPS:-}" = 1 ]; then groups_in "$tier"; exit 0; fi

mode=run
[ "${GATE_LIST:-}" = 1 ] && mode=list
mkdir -p build
checks
exit $fail
