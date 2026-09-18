## Self-test for dev/xver/reduce.R against one real version, before anything is
## built on top of it. Checks the three properties the harness depends on:
##
##   1. every step reduces without hitting the no-extractor error
##   2. no payload leaks an environment / closure / S4 object (xv_assert_plain)
##   3. reduction is deterministic -- reduce the same object twice, get the same
##      thing, so a later cross-version diff is about the objects and not about
##      the reducer
##
## Usage: XVER_VERSION=main Rscript --vanilla dev/xver/21_test_reduce.R

source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "common.R"))
source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "reduce.R"))

which_ver <- Sys.getenv("XVER_VERSION", "main")
spec <- switch(which_ver,
  main = list(lib = Sys.getenv("XVER_LIB_MAIN"), version = "1.2.2"),
  alex = list(lib = Sys.getenv("XVER_LIB_ALEX"), version = "2.0.0"))
xv_load_spatialge(spec$lib, spec$version)
xv_source_driver()
xv_shim_fixture_image()
setwd(xv_scratch(paste0("testreduce-", which_ver)))

cat(sprintf("== reducer self-test: %s ==\n", which_ver))
res <- bl_pipeline(cores = 1, verbose = FALSE)

fails <- 0L
for(nm in names(res)){
  rec <- res[[nm]]
  ## bl_pipeline uses bl_capture, which has no warnings/messages fields; adapt.
  rec$warnings <- character(0); rec$messages <- character(0)

  r1 <- tryCatch(xv_reduce_step(rec, canon = FALSE),
                 error = function(e) structure(conditionMessage(e), class = "xv_err"))
  if(inherits(r1, "xv_err")){
    cat(sprintf("  REDUCE-FAIL %-22s %s\n", nm, r1)); fails <- fails + 1L; next
  }
  leak <- tryCatch({ xv_assert_plain(r1, nm); NA_character_ },
                   error = function(e) conditionMessage(e))
  if(!is.na(leak)){
    cat(sprintf("  LEAK        %-22s %s\n", nm, leak)); fails <- fails + 1L; next
  }
  rc <- tryCatch(xv_reduce_step(rec, canon = TRUE),
                 error = function(e) structure(conditionMessage(e), class = "xv_err"))
  if(inherits(rc, "xv_err")){
    cat(sprintf("  CANON-FAIL  %-22s %s\n", nm, rc)); fails <- fails + 1L; next
  }
  ## Determinism: a reducer that is not a pure function of its input would make
  ## every cross-version diff untrustworthy.
  r2 <- xv_reduce_step(rec, canon = FALSE)
  d <- waldo::compare(r1, r2, tolerance = NULL)
  if(length(d)){
    cat(sprintf("  NONDETERM   %-22s %d diff(s)\n", nm, length(d))); fails <- fails + 1L; next
  }
  nleaf_raw <- length(xv_flatten_list(r1))
  nleaf_can <- length(xv_flatten_list(rc))
  cat(sprintf("  ok          %-22s raw=%-6d canon=%-6d\n", nm, nleaf_raw, nleaf_can))
}

cat(sprintf("\n%s: %d failure(s)\n", which_ver, fails))
if(fails > 0) quit(status = 1)
