#!/usr/bin/env bash
## Preconditions for the cross-version comparison. Every check here is one that,
## if silently false, would make a later "no difference" result meaningless
## rather than merely wrong. Fails loudly and early instead.
##
## Usage:
##   source /opt/anaconda3/etc/profile.d/conda.sh && conda activate spatialge-dev
##   bash dev/xver/00_gate.sh

set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "$XVER_REPO"

fail() { echo "GATE FAIL: $*" >&2; exit 1; }
ok()   { echo "  ok  $*"; }

echo "== xver gate =="

## 1. Right R. System R has a different package set and a libcurl that fails
## differently, so a run against it proves nothing about this environment.
command -v Rscript >/dev/null || fail "no Rscript on PATH"
case "$(command -v Rscript)" in
  *spatialge-dev*) ok "Rscript: $(command -v Rscript)" ;;
  *) fail "Rscript is $(command -v Rscript), not the spatialge-dev conda env -- run 'conda activate spatialge-dev' first" ;;
esac
Rscript --vanilla -e 'if(getRversion() < "4.5.0") stop("R too old: ", getRversion())' \
  || fail "R version check failed"
ok "R version $(Rscript --vanilla -e 'cat(as.character(getRversion()))')"

## 2. Clean tree, for anything that can change behavior. The harness installs alex
## from `git archive HEAD`, so an uncommitted change to R/, src/, inst/, NAMESPACE
## or DESCRIPTION would be silently absent from the comparison -- the run would
## measure the last commit while appearing to measure the working tree.
##
## Doc-only edits are a warning, not a failure: NEWS.md and dev/ARCHITECTURE.md
## carry no runtime behavior, and this harness's own output is a set of proposed
## edits to exactly those two files, so blocking on them would make the gate
## unrunnable immediately after a successful run.
BEHAVIORAL="$(git status --porcelain -- R/ src/ inst/ NAMESPACE DESCRIPTION tests/)"
[ -z "$BEHAVIORAL" ] || fail "uncommitted changes that affect behavior -- commit or stash them, or the archived HEAD will not match what you think you are testing:
$BEHAVIORAL"
ok "no uncommitted changes under R/ src/ inst/ tests/ NAMESPACE DESCRIPTION"

OTHER="$(git status --porcelain -- ':!plans/' ':!dev/xver/' ':!R/' ':!src/' ':!inst/' ':!NAMESPACE' ':!DESCRIPTION' ':!tests/')"
if [ -n "$OTHER" ]; then
  echo "  WARN  doc/other files modified (not behavior-affecting, comparison still valid):"
  printf '%s\n' "$OTHER" | sed 's/^/          /'
fi

## 3. The refs resolve, and alex is a strict descendant of main. If main ever
## gains commits alex does not have, "alex vs main" stops being the question
## this harness answers.
SHA_MAIN="$(git rev-parse "$XVER_REF_MAIN")"
SHA_ALEX="$(git rev-parse "$XVER_REF_ALEX")"
[ "$SHA_ALEX" = "$(git rev-parse HEAD)" ] || fail "$XVER_REF_ALEX is not checked out (HEAD is $(git rev-parse --short HEAD))"
git merge-base --is-ancestor "$SHA_MAIN" "$SHA_ALEX" \
  || fail "$XVER_REF_MAIN is not an ancestor of $XVER_REF_ALEX -- the branches have diverged"
ok "main=${SHA_MAIN:0:7}  alex=${SHA_ALEX:0:7}  (ancestor, $(git rev-list --count "$SHA_MAIN..$SHA_ALEX") commits apart)"

## 4. Local refs match the remote. The stated task is to compare the branch *as
## it appears on GitHub*, so a stale local ref would answer a different question.
REMOTE_REFS="$(git ls-remote origin "refs/heads/$XVER_REF_MAIN" "refs/heads/$XVER_REF_ALEX" 2>/dev/null || true)"
if [ -n "$REMOTE_REFS" ]; then
  R_MAIN="$(printf '%s\n' "$REMOTE_REFS" | awk -v r="refs/heads/$XVER_REF_MAIN" '$2==r{print $1}')"
  R_ALEX="$(printf '%s\n' "$REMOTE_REFS" | awk -v r="refs/heads/$XVER_REF_ALEX" '$2==r{print $1}')"
  [ "$R_MAIN"  = "$SHA_MAIN" ] || fail "origin/$XVER_REF_MAIN is ${R_MAIN:0:7}, local is ${SHA_MAIN:0:7}"
  [ "$R_ALEX"  = "$SHA_ALEX" ] || fail "origin/$XVER_REF_ALEX is ${R_ALEX:0:7}, local is ${SHA_ALEX:0:7}"
  ok "both refs match origin"
else
  echo "  WARN  could not reach origin -- comparing local refs only"
fi

## 5. Fidelity rule 12: this variable clamps count_cores() on alex only.
[ -z "${_R_CHECK_LIMIT_CORES_:-}" ] || fail "_R_CHECK_LIMIT_CORES_ is set"
ok "_R_CHECK_LIMIT_CORES_ unset"

## 6. detectCores() must not be NA. count_cores() divides it, and dev/environment
## /README.md records it returning NA in sandboxed shells -- which would make
## every auto-detected core count NA rather than erroring.
Rscript --vanilla -e 'n <- parallel::detectCores(); if(is.na(n)) stop("detectCores() is NA"); cat("detectCores:", n, "\n")' \
  || fail "detectCores() is NA in this shell"

## 7. Every package the harness itself needs, in the shared library. writexl and
## waldo are in alex's Suggests only, so installing main would not reveal a gap
## -- and writexl builds the xlsx fixture that is the only probe for delta 6.
Rscript --vanilla -e '
  need <- c("waldo","writexl","png","jsonlite","Matrix","ggplot2","spdep","sf","sp",
            "ComplexHeatmap","sctransform","dplyr","tibble","tidyr","Rcpp")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
  if(length(miss)) stop("missing packages: ", paste(miss, collapse=", "))
  cat("harness deps: all", length(need), "present\n")' \
  || fail "harness dependency missing"

## 8. hdf5r/rhdf5 and arrow gate only some fixture formats, so a miss is a
## scenario skip rather than a hard stop -- but it must be visible, not silent.
Rscript --vanilla -e '
  for(p in c("hdf5r","rhdf5","arrow")){
    cat(sprintf("  %-8s %s\n", p, if(requireNamespace(p, quietly=TRUE)) "present" else "MISSING (some S1 formats will be skipped)"))
  }'

## 9. TLS. Scenario S12 (real Thrane data) needs a working download; per the plan
## a failure here is a recorded COVERAGE_GAP, not a silent skip, so surface it now.
if Rscript --vanilla -e '
  h <- tryCatch(curlGetHeaders("https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=")[1],
                error=function(e) paste("ERROR:", conditionMessage(e)))
  cat(h)
  if(grepl("^ERROR", h)) quit(status=1)' 2>&1 | sed 's/^/  tls: /'; then
  ok "TLS gate passed (S12 can run)"
else
  echo "  WARN  TLS gate FAILED -- S12 will be recorded as COVERAGE_GAP, not skipped"
  echo "        see the 'Corporate network / TLS interception' section of dev/environment/README.md"
fi

echo "== gate passed =="
