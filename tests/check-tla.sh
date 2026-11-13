#!/bin/bash
# Model-checks the TLA+ specs in tla/ with TLC over a matrix of constants.
# docs/testing.md says what each spec models. Each run checks the spec's
# invariants, no deadlock, and its liveness properties. Needs Java; downloads
# tla2tools.jar TLA_VERSION into data/, checked against TLA_SHA256.
#
# Usage: tests/check-tla.sh [Spec ...]   (default: every spec in SPECS)
# shellcheck disable=SC2329  # the function named after each spec is called through SPECS
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SPECS=(ParallelOrdered BlockReader SlidingWindow FatalExit BgenCleanup VcfIt LineReader Bref3It RefSeqGroup WindowChain VcfIndex BgenWriter PhasingSchedule ClusterLifecycle HmmPass CompHapTracker)
[ $# -gt 0 ] || set -- "${SPECS[@]}"
for spec in "$@"; do
  known=0
  for s in "${SPECS[@]}"; do [ "$s" != "$spec" ] || known=1; done
  [ $known = 1 ] || { echo "FAIL unknown spec $spec"; exit 1; }
done
TLA_VERSION=v1.7.4
TLA_SHA256=936a262061c914694dfd669a543be24573c45d5aa0ff20a8b96b23d01e050e88
JAR="$ROOT/data/tla2tools-$TLA_VERSION.jar"
SHA=$(command -v sha256sum || echo "shasum -a 256")
# As in tests/fetch-fixtures.sh: run java rather than look for it.
java -version > /dev/null 2>&1 || {
  echo "FAIL java cannot run, and TLC needs it: put a JDK first on PATH" >&2
  exit 1
}

mkdir -p "$ROOT/data"
if [ ! -f "$JAR" ]; then
  curl -LSf -o "$JAR.part" "https://github.com/tlaplus/tlaplus/releases/download/$TLA_VERSION/tla2tools.jar" \
    || { echo "FAIL download tla2tools.jar $TLA_VERSION"; exit 1; }
  mv "$JAR.part" "$JAR"
fi
got=$($SHA "$JAR" | cut -d' ' -f1)
[ "$got" = "$TLA_SHA256" ] || { echo "FAIL $JAR has sha256 $got, want $TLA_SHA256"; exit 1; }

OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT
fail=0
run=0

# check <Spec> <label> <invariants> <properties> <<< "<constant assignments>":
# one TLC run of tla/<Spec>.tla.
check() {
  local spec=$1 label=$2
  run=$((run + 1))
  { echo "CONSTANTS"; cat; echo "SPECIFICATION Spec"; echo "INVARIANTS $3"; echo "PROPERTIES $4"; } > "$OUT/MC.cfg"
  cp "$ROOT/tla/$spec.tla" "$OUT/"
  if java -XX:+UseParallelGC -cp "$JAR" tlc2.TLC -workers 4 -metadir "$OUT/states-$run" \
      -config "$OUT/MC.cfg" "$OUT/$spec.tla" > "$OUT/tlc.log" 2>&1 \
      && grep -q "No error has been found" "$OUT/tlc.log"; then
    echo "PASS $spec $label $(sed -n 's/.* \([0-9,]*\) distinct states found, 0 states left on queue.*/\1/p' "$OUT/tlc.log" | tail -1) distinct states found"
  else
    echo "FAIL $spec $label: $(grep -m1 "^Error:" "$OUT/tlc.log")"; grep -v "^\s*$" "$OUT/tlc.log" | tail -60; fail=1
  fi
}

# worker_set <n>: w1,...,w<n>
worker_set() {
  seq -f "w%g" 1 "$1" | paste -sd, -
}

ParallelOrdered() {
  # workers items window
  for model in "1 3 1" "2 4 1" "2 5 2" "3 5 2" "3 6 3" "2 3 5" "3 8 3" "4 7 2"; do
    read -r nw ni win <<< "$model"
    check ParallelOrdered "workers=$nw items=$ni window=$win" \
      "TypeOK WindowBound NoOverwrite ConsumesOwnItem InOrder" "Terminates AllConsumed" <<EOF
  Workers = {$(worker_set "$nw")}
  NItems = $ni
  Window = $win
EOF
  done
}

BlockReader() {
  local invariants="TypeOK SlotsPartition ReaderSlotExclusive InOrder"
  # The code has BLOCK_READER_SLOTS = 4; Slots = 1 is the boundary where a
  # missing broadcast stalls the pipeline.
  for slots in 1 2 3 4; do
    for nb in 0 1 3 5; do
      check BlockReader "slots=$slots batches=$nb close=any" \
        "$invariants" "CloseTerminates WaitEnds Progress" <<EOF
  Slots = $slots
  NBatches = $nb
  MayClose = TRUE
EOF
      check BlockReader "slots=$slots batches=$nb close=never" \
        "$invariants" "WaitEnds SeesSentinel" <<EOF
  Slots = $slots
  NBatches = $nb
  MayClose = FALSE
EOF
    done
  done
}

SlidingWindow() {
  # windows failat
  for model in "1 {}" "1 {1}" "2 {}" "2 {1,2}" "3 {}" "3 {1,2,3}" "4 {}" "4 {1,2,3,4}" "4 {3}"; do
    read -r nw failat <<< "$model"
    check SlidingWindow "windows=$nw failat=$failat" \
      "TypeOK InOrder AheadNotGiven NoLeak ErrorAfterPrefix ErrorNoAhead" \
      "CloseTerminates NextReturns ExitOnlyOnError" <<EOF
  NWindows = $nw
  FailAt = $failat
EOF
  done
}

FatalExit() {
  for nw in 1 2; do
    check FatalExit "workers=$nw" \
      "TypeOK" "Terminates LockReleased MainNotBlocked DeferredRaised" <<EOF
  Workers = {$(worker_set "$nw")}
EOF
  done
}

BgenCleanup() {
  check BgenCleanup "members=3" \
    "TypeOK NoPartialAfterCleanup CompletedSurvive" "Terminates" <<EOF
  Members = {bgen, info, sample}
EOF
}

VcfIt() {
  # The code has VCF_IT_BATCH = 1024; Init chooses every line's kind and the
  # failing read, so each run covers every input of N lines.
  for model in "0 1" "1 1" "3 1" "3 2" "4 2" "4 3" "5 2"; do
    read -r n batch <<< "$model"
    check VcfIt "lines=$n batch=$batch" "TypeOK OnDemandPrefix ReadErrKept" "Completes" <<EOF
  N = $n
  Batch = $batch
EOF
  done
}

LineReader() {
  # Init chooses every input of up to K bytes over {a, CR, LF} and every
  # split into chunks, so one run per mode covers all of them.
  for bgzip in TRUE FALSE; do
    check LineReader "bytes<=6 bgzip=$bgzip" "TypeOK OutPrefix EofReadOnce" "Completes" <<EOF
  K = 6
  BgzipLines = $bgzip
EOF
  done
}

Bref3It() {
  # blocks records: Init chooses up to R records per block, their filter
  # verdicts, and where the file is cut short.
  for model in "0 0" "1 2" "2 2" "3 3"; do
    read -r nb r <<< "$model"
    check Bref3It "blocks=$nb records<=$r" "TypeOK OutPrefix NoReadPastEnd" "Completes" <<EOF
  NB = $nb
  R = $r
EOF
  done
}

RefSeqGroup() {
  # records haplotypes max_nseq. 4 haplotypes with max_nseq 3 is the
  # default for 2 samples; 3 with 2 makes the coder fill sooner; 2 with 3,
  # the default for 1 sample, never fills it.
  for model in "3 3 2" "4 3 2" "2 4 3" "4 2 3"; do
    read -r n h m <<< "$model"
    check RefSeqGroup "records=$n haps=$h maxnseq=$m" \
      "TypeOK OutPrefix PlaceholdersFilled CoderSeqs EmptyCoderFresh GroupOneChrom NoAssert" "Completes" <<EOF
  N = $n
  H = $h
  MaxNSeq = $m
  C = 2
EOF
  done
}

WindowChain() {
  # mode targ-records ref-records chromosomes positions alleles window overlap window-markers
  for model in "targ 3 0 2 3 1 2 1 2" "targ 6 0 2 3 1 2 1 4" "targ 5 0 2 3 1 3 2 4" \
      "ref 3 3 2 2 2 2 1 4" "ref 2 3 2 3 2 2 1 4" "ref 3 4 2 3 1 2 1 4" "ref 2 3 2 3 2 3 2 4" \
      "impute 3 3 2 2 2 2 1 4" "impute 2 3 2 2 2 2 1 2" "impute 2 3 2 3 2 2 1 4" "impute 3 4 1 3 2 2 1 4"; do
    read -r mode nt nr nc p a w o wm <<< "$model"
    local invariants="TypeOK CarrySuffix Paired ChromsBesideCarry StallIsReal NoSilentRepeat"
    [ "$mode" != targ ] || invariants="TypeOK WinInterval Covered OneChrom CarrySuffix BoundaryEmpty DoneAtLast StallIsReal NoSilentRepeat"
    check WindowChain "mode=$mode targ=$nt ref=$nr chroms=$nc pos=$p alleles=$a window=$w overlap=$o markers=$wm" \
      "$invariants" "Terminates" <<EOF
  Mode = "$mode"
  NT = $nt
  NR = $nr
  NC = $nc
  P = $p
  A = $a
  W = $w
  O = $o
  WM = $wm
EOF
  done
}

VcfIndex() {
  # Records whose CHROM starts with '#' are the only lines put can pass to the
  # '#' branch; Init chooses each line's kind.
  for n in 0 1 3 5 7; do
    check VcfIndex "lines=$n" "TypeOK IndexedLikeTabix HeaderEndLikeTabix" "Completes" <<EOF
  N = $n
  Kinds = {"hash", "rec"}
EOF
  done
}

BgenWriter() {
  # BgenWriter.tla instantiates BgenCleanup, and check copies only the spec it runs.
  cp "$ROOT/tla/BgenCleanup.tla" "$OUT/"
  # records chromosomes autosomes
  for model in "0 1 {1}" "2 2 {1}" "3 2 {1}" "2 3 {1,2}"; do
    read -r n nc au <<< "$model"
    check BgenWriter "records=$n chroms=$nc autosomes=$au" \
      "TypeOK WrittenAutosomes CheckedBeforeWrite HeaderCountsWritten CompleteAfterClose Plink2EmptyExits SuccessKeepsAll" \
      "Terminates" <<EOF
  N = $n
  NChroms = $nc
  Autosomes = $au
EOF
  done
}

PhasingSchedule() {
  # The defaults are burnin=3 iterations=12; par.c allows both down to 1, and
  # iterations=1 is where lr_threshold divides by zero after the last run.
  for model in "1 1" "1 2" "2 1" "3 1" "3 3" "4 4" "3 12" "5 6"; do
    read -r burnin its <<< "$model"
    check PhasingSchedule "burnin=$burnin iterations=$its" \
      "TypeOK RunBound FollowsSchedule EarlyExitOnLowRate RunThresholds LastRunOne EmInBurnin Stage2AfterStage1 PhasedRunsNothing" \
      "Terminates" <<EOF
  Burnin = $burnin
  Iterations = $its
EOF
  done
}

ClusterLifecycle() {
  # clusters passes; Init chooses every type sequence set_clusters can produce.
  for model in "3 2" "4 2" "5 2" "5 3" "6 2"; do
    read -r n passes <<< "$model"
    check ClusterLifecycle "clusters=$n passes=$passes" \
      "TypeOK CountsMatch NoExit MaskMatchesRef" "Transitions Terminates" <<EOF
  N = $n
  MaxPasses = $passes
EOF
  done
}

HmmPass() {
  # clusters samples; two samples on one phase_baum2 check begin_sample's reset.
  for model in "3 1" "4 1" "5 1" "3 2"; do
    read -r n ns <<< "$model"
    check HmmPass "clusters=$n samples=$ns" \
      "TypeOK FwdLockstep BwdLockstep IndicesMatch AheadUnchanged SwapsCountFlips RateSums SwapRange NoMarkInBurnin" \
      "LastRunPhasesAll" <<EOF
  N = $n
  S = $ns
EOF
  done
}

CompHapTracker() {
  # haps max_states max_step min_steps calls; min_steps > max_step is INT_MAX.
  # Four slots is the smallest heap whose poll sifts past a node with two children.
  for model in "3 1 2 1 4" "3 2 2 3 5" "3 2 2 1 5" "3 3 2 2 6" "4 3 3 9 5" "5 4 1 3 5"; do
    read -r nh ms st mn ev <<< "$model"
    check CompHapTracker "haps=$nh max_states=$ms steps=$st min_steps=$mn calls=$ev" \
      "TypeOK SizeBound SlotsArePrefix HeapOrdered NilIffSlot KeyBehindStep EvictsLeast TrackedNoChange" \
      "RepairEnds" <<EOF
  H = $nh
  MaxStates = $ms
  MaxStep = $st
  MinSteps = $mn
  MaxEvents = $ev
EOF
  done
}

for spec in "$@"; do "$spec"; done
exit $fail
