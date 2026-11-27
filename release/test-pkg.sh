#!/bin/bash
# Smoke-test the actual installer on a disposable macOS runner.
set -euo pipefail
pkg=$1
recipe=$(cd "$(dirname "$0")/conda" && pwd)
sudo installer -pkg "$pkg" -target /
codesign --verify --strict --verbose=2 /usr/local/bin/fast-beagle
brew list htslib > /dev/null 2>&1 || brew install htslib
/usr/local/bin/fast-beagle gt="$recipe/small.vcf" out="$RUNNER_TEMP/small" nthreads=2
test "$(bgzip -dc "$RUNNER_TEMP/small.vcf.gz" | grep -vc '^#')" = 40
echo "${pkg##*/}: installed and run" >> "$GITHUB_STEP_SUMMARY"
