#!/bin/bash
# Run build/beagle on the chr20 imputation from tests/bench/fetch-chr20.sh at
# 18 threads and compare its output to the hash that the release jar wrote on
# the same input (docs/perf-baseline.md). The gate's cases span 100 kb; this
# check runs a whole chromosome, so the nightly CI run runs it outside the
# pre-merge gate.
#
# Usage: tests/check-chr20.sh <dir from tests/bench/fetch-chr20.sh>
set -uo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DIR=$1
WANT=38d4ed96cec67a44
SHA=$(command -v sha256sum || echo "shasum -a 256")
OUT=$(mktemp -d)
trap 'rm -rf "$OUT"' EXIT

"$ROOT/build/beagle" ref="$DIR/ref.vcf.gz" gt="$DIR/target.vcf.gz" map="$DIR/chr20.map" \
  out="$OUT/chr20" nthreads=18 > "$OUT/run.log" 2>&1
rc=$?
hash=
[ -f "$OUT/chr20.vcf.gz" ] && hash=$(gzip -dc "$OUT/chr20.vcf.gz" | grep -v '^##source\|^##filedate' | $SHA | cut -c1-16)
if [ "$rc" -eq 0 ] && [ "$hash" = "$WANT" ]; then
  echo "PASS chr20 nthreads=18 exit=0 $hash"
else
  echo "FAIL chr20 nthreads=18 exit=$rc ${hash:-no output}, want exit=0 $WANT"
  tail -5 "$OUT/run.log"
  exit 1
fi
