#!/usr/bin/env bash
## Installs both refs under comparison into two separate libraries, from two
## `git archive` exports.
##
## `git archive` rather than `git worktree add`: the export contains no .git, so
## nothing the harness subsequently runs can commit, checkout or reset the real
## repository. A worktree would also write refs and admin files into .git, mark
## `main` as in-use, and need a cleanup step that can leave a half state on abort.
##
## Both refs are exported the same way on purpose. Installing alex from the live
## working tree instead would make the two sides asymmetric in how they were
## built, which is exactly the kind of uncontrolled variable this whole exercise
## exists to eliminate.
##
## Usage: bash dev/xver/10_build_libs.sh

set -euo pipefail
. "$(dirname "${BASH_SOURCE[0]}")/env.sh"
cd "$XVER_REPO"

build_one() {
  local ref="$1" srcdir="$2" libdir="$3" want_ver="$4"
  echo "== $ref -> $libdir =="

  rm -rf "$srcdir" "$libdir"
  mkdir -p "$srcdir" "$libdir"
  git archive "$ref" | tar -x -C "$srcdir"

  ## main has src/*.o and src/spatialGE.so committed. `git archive` stamps every
  ## file with the commit time, so make can decide those objects are current and
  ## relink a .so built against an unknown R/ABI on someone else's machine.
  ## Removing them is correctness, not hygiene; --preclean below is the backstop.
  local stale
  stale="$(find "$srcdir/src" -maxdepth 1 \( -name '*.o' -o -name '*.so' -o -name '*.dylib' \) 2>/dev/null | wc -l | tr -d ' ')"
  if [ "$stale" != "0" ]; then
    echo "  removing $stale committed build artifact(s) from src/"
    find "$srcdir/src" -maxdepth 1 \( -name '*.o' -o -name '*.so' -o -name '*.dylib' \) -delete
  fi

  ## --without-keep.source keeps srcrefs out of function bodies, so nothing that
  ## serializes a closure carries file paths and line numbers that differ between
  ## the two exports for reasons unrelated to behavior. (The flag is spelled
  ## --without-, not --no-; the latter is accepted-and-ignored with a warning.)
  R CMD INSTALL --preclean --no-multiarch --without-keep.source \
    -l "$libdir" "$srcdir" > "$XVER_OUT/install-$(basename "$libdir").log" 2>&1 \
    || { echo "INSTALL FAILED for $ref -- see $XVER_OUT/install-$(basename "$libdir").log" >&2
         tail -40 "$XVER_OUT/install-$(basename "$libdir").log" >&2; exit 1; }

  ## Assert we got the version we asked for, from the library we asked for.
  ## Without this the whole comparison can silently run alex against alex.
  R_LIBS="$libdir" Rscript --vanilla -e "
    lib <- normalizePath('$libdir')
    suppressMessages(library(spatialGE, lib.loc=lib))
    got_lib <- normalizePath(dirname(find.package('spatialGE')))
    if(got_lib != lib) stop('loaded spatialGE from ', got_lib, ' not ', lib)
    v <- as.character(packageVersion('spatialGE'))
    if(v != '$want_ver') stop('expected version $want_ver, got ', v)
    cat('  ok  version', v, 'from', got_lib, '\n')
    cat('  ok ', length(getNamespaceExports('spatialGE')), 'exports\n')
  " || { echo "VERIFY FAILED for $ref" >&2; exit 1; }

  ## Warnings in the install log are worth seeing once: main has never been
  ## compiled in this conda toolchain, so a new warning here is information.
  local nwarn
  nwarn="$(grep -ci 'warning' "$XVER_OUT/install-$(basename "$libdir").log" || true)"
  echo "  install log: $nwarn line(s) matching 'warning'"
}

build_one "$XVER_REF_MAIN" "$XVER_SRC_MAIN" "$XVER_LIB_MAIN" "1.2.2"
build_one "$XVER_REF_ALEX" "$XVER_SRC_ALEX" "$XVER_LIB_ALEX" "2.0.0"

## The comparison is only meaningful if both sides ran against the same
## dependency versions, so record them once and have the comparator assert the
## two payload headers agree.
Rscript --vanilla -e '
  d <- read.dcf(file.path(Sys.getenv("XVER_SRC_ALEX"), "DESCRIPTION"))
  fields <- intersect(c("Imports","Depends","LinkingTo"), colnames(d))
  raw <- unlist(strsplit(paste(d[1, fields], collapse=","), ","))
  pkgs <- trimws(sub("[(].*", "", raw))
  pkgs <- setdiff(pkgs[nzchar(pkgs)], c("R","methods","utils","stats","graphics","grDevices","tools","parallel"))
  vers <- vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error=function(e) NA_character_), character(1))
  out <- data.frame(package=pkgs, version=unname(vers), stringsAsFactors=FALSE)
  fp <- file.path(Sys.getenv("XVER_OUT"), "dep-versions.csv")
  write.csv(out, fp, row.names=FALSE)
  cat("recorded", nrow(out), "dependency versions ->", fp, "\n")
  if(any(is.na(out$version))) cat("  NOTE not installed:", paste(out$package[is.na(out$version)], collapse=", "), "\n")'

echo "== both libraries built =="
