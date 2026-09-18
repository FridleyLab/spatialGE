#!/usr/bin/env bash
## Compares the S3 stdout transcripts, which is the only way to observe delta 4
## (b7e9c73, STclust clobbering its own verbose argument).
##
## Why this is not part of compare.R: sge_worker_msg keeps
## system(sprintf('echo "%s"', msg)) verbatim (fidelity rule 11), which writes to
## the PROCESS's stdout. capture.output() cannot see it, so it never reaches the
## R-level payload -- it only exists in the producer's redirected log. The producer
## brackets each STclust call with an XVMARK line so the sections can be isolated.
##
## Expected result: main prints the same number of lines for verbose=TRUE and
## verbose=FALSE (the argument is ignored); alex prints nothing for verbose=FALSE.
##
## Usage: bash dev/xver/40_stdout_probe.sh [log-suffix]   (default: extra)

set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/env.sh"

SUFFIX="${1:-extra}"
echo "== S3 verbose transcript (delta 4: b7e9c73) =="
printf "%-6s %-16s %s\n" "ver" "verbose" "lines printed"
fail=0
for v in main alex; do
  log="$XVER_OUT/$SUFFIX-$v.log"
  if [ ! -f "$log" ]; then echo "  missing $log -- run produce.R with S3 in XVER_SCENARIOS"; exit 1; fi
  for tag in TRUE FALSE; do
    n=$(sed -n "/XVMARK s3_stclust_verbose$tag/,/XVMARK \/s3_stclust_verbose$tag/p" "$log" \
        | grep -v XVMARK | grep -c . || true)
    printf "  %-4s %-16s %s\n" "$v" "$tag" "$n"
    eval "n_${v}_${tag}=$n"
  done
done

echo
# shellcheck disable=SC2154
if [ "$n_main_TRUE" -gt 0 ] && [ "$n_main_FALSE" -gt 0 ] \
   && [ "$n_alex_TRUE" -gt 0 ] && [ "$n_alex_FALSE" -eq 0 ]; then
  echo "VERIFIED: main ignores verbose=FALSE and prints anyway; alex-refactor honors it."
else
  echo "UNEXPECTED: the pattern main(TRUE>0,FALSE>0) / alex(TRUE>0,FALSE==0) did not hold."
  echo "Do not treat delta 4 as verified without explaining this."
  fail=1
fi
exit $fail
