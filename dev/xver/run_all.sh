#!/usr/bin/env bash
## Orchestrator. Runs every producer in its own process, then the comparator.
##
## Producers must be separate PROCESSES, not two calls in one session: a second
## call in the same session inherits mutated RNG state, a populated scratch dir
## and warm caches, so it would understate the noise floor -- and the two reps
## would also collide on bl_pipeline's per-cores image directory.
##
## Usage: bash dev/xver/run_all.sh [scenario-list]
##   e.g. bash dev/xver/run_all.sh S0,S1,S2,S4,S5,S6,S7,S8,S9,S10,S11,S13,S14

set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "$XVER_REPO"

SCEN="${1:-S0}"
echo "== scenarios: $SCEN =="

for ver in main alex; do
  case "$ver" in
    main) LIB="$XVER_LIB_MAIN" ;;
    alex) LIB="$XVER_LIB_ALEX" ;;
  esac
  for rep in 1 2; do
    echo
    echo "### producer $ver rep$rep"
    log="$XVER_OUT/produce-$ver-rep$rep.log"
    # R_LIBS is set so the version library wins; the conda library still supplies
    # the ~366 shared dependencies, identically for both sides.
    if XVER_VERSION="$ver" XVER_REP="$rep" XVER_SCENARIOS="$SCEN" \
       R_LIBS="$LIB" Rscript --vanilla dev/xver/produce.R > "$log" 2>&1; then
      grep -E "^   .*(err|REDUCE-ERR)|^wrote" "$log" | sed 's/^/    /' || true
    else
      echo "    PRODUCER FAILED -- tail of $log:"; tail -25 "$log"; exit 1
    fi
  done
done

echo
echo "### comparator"
# Deliberately no version library on the path: the comparator must not be able to
# load spatialGE, so any payload leak fails loudly instead of deserializing
# against whichever class table happens to be available.
Rscript --vanilla dev/xver/compare.R 2>&1 | tee "$XVER_OUT/compare.log"
