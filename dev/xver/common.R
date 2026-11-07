## Shared setup for every producer process. Sourced, never run directly.
##
## Everything here exists to remove a variable that would otherwise differ
## between the two versions for reasons unrelated to their code. If a producer
## skips this file, a later "the versions differ" result is not trustworthy.

xv_repo <- Sys.getenv("XVER_REPO")
if(!nzchar(xv_repo)) stop("XVER_REPO is unset -- source dev/xver/env.sh first")

## ---- process determinism -------------------------------------------------

## Canonical sorting must not depend on the ambient locale: the fixture
## deliberately contains 'Gene$special' and '1BadGene', whose collation order
## differs between C and en_US.
Sys.setlocale("LC_COLLATE", "C")

## main leaves several mclapply call sites with no mc.cores, falling through to
## getOption("mc.cores", 2L) (fidelity rule 3). Pinning it means the two versions
## fork the same way, so a diff is about their code rather than their topology.
options(mc.cores = 1L)

## Fidelity rule 12: this clamps count_cores() on alex-refactor but not on main.
Sys.unsetenv("_R_CHECK_LIMIT_CORES_")

## pseudobulk_heatmap() ends in ComplexHeatmap::draw(), which needs a device and
## would otherwise litter Rplots.pdf into the working directory.
grDevices::pdf(NULL)

options(stringsAsFactors = FALSE, warn = 1)

## ---- version assertions --------------------------------------------------

## The single most important check in the harness: without it, a mistake in
## R_LIBS silently runs alex against alex and reports a clean comparison.
xv_load_spatialge <- function(lib, want_version){
  lib <- normalizePath(lib, mustWork = TRUE)
  suppressMessages(library(spatialGE, lib.loc = lib, character.only = FALSE))
  got <- normalizePath(dirname(find.package("spatialGE")))
  if(!identical(got, lib)){
    stop("loaded spatialGE from ", got, " but expected ", lib)
  }
  v <- as.character(utils::packageVersion("spatialGE"))
  if(!identical(v, want_version)){
    stop("expected spatialGE ", want_version, " but loaded ", v)
  }
  invisible(list(lib = lib, version = v))
}

## ---- fixtures and pipeline, sourced from the alex working tree ------------

## Both producers deliberately source ONE copy of these, from the working tree,
## never from the installed package. main's tree contains neither file, and
## inst/testdata/ installs to a different path anyway -- so a single shared copy
## is the only way to guarantee both versions run an identical call sequence.
xv_source_driver <- function(){
  source(file.path(xv_repo, "inst", "testdata", "make_fixtures.R"))
  source(file.path(xv_repo, "dev", "baseline", "pipeline.R"))
  source(file.path(xv_repo, "dev", "baseline", "compare_baseline.R"))

  ## bl_pipeline() must resolve the shimmed fixture_write_image() below, which
  ## only happens if it was sourced into the global environment. Assert rather
  ## than assume, because a silent failure here reintroduces the RNG dependency
  ## the shim exists to remove.
  if(!identical(environment(bl_pipeline), globalenv())){
    stop("bl_pipeline was not sourced into globalenv(); the fixture shim will not apply")
  }
  invisible(TRUE)
}

## dev/baseline/pipeline.R writes the load_images PNGs in a loop that sits
## OUTSIDE bl_capture()'s set.seed(), so the pixel values depend on the RNG state
## left behind by the preceding 25 steps. That state can legitimately differ
## between versions, which would show up as a load_images diff that has nothing
## to do with load_images. Seeding from the output path makes each image a pure
## function of where it is written.
##
## Because bl_capture() re-seeds before every step, this loop is the only
## RNG-state-order dependency in the whole pipeline.
xv_shim_fixture_image <- function(){
  assign("fixture_write_image", function(fp, side = 20){
    set.seed(20260115L + sum(utf8ToInt(basename(fp))))
    arr <- array(stats::runif(side * side * 3), dim = c(side, side, 3))
    png::writePNG(arr, fp)
  }, envir = globalenv())
  invisible(TRUE)
}

## ---- isolated scratch space ----------------------------------------------

## Never bare tempdir(): tests/testthat/test-tests.R does
## unlink(tempdir(), recursive=TRUE), and anything sharing a root with that kind
## of call is one stray source() away from losing its fixtures mid-run.
xv_scratch <- function(tag){
  root <- file.path(Sys.getenv("XVER_ROOT"), "scratch", tag)
  unlink(root, recursive = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}

## ---- capture ------------------------------------------------------------

## bl_capture() records only errors, so a regression that adds, drops or
## downgrades a warning is invisible to it -- and delta 8 (ab76c49) is precisely
## a warning path that errored. This wrapper keeps warnings and messages as
## first-class payload.
xv_capture <- function(expr, seed = 20260115L){
  warns <- character(0)
  msgs  <- character(0)
  set.seed(seed)
  val <- withCallingHandlers(
    tryCatch(
      list(ok = TRUE, value = expr),
      error = function(e) list(ok = FALSE, value = conditionMessage(e),
                                error_classes = class(e))
    ),
    warning = function(w){
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m){
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  val$warnings <- warns
  val$messages <- msgs
  val
}
