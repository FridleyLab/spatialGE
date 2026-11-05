##
# Golden-baseline regression tests. Recomputes dev/baseline/pipeline.R's
# bl_pipeline() fresh and diffs every result against the corresponding
# tests/golden/<fn>__cores<N>.rds captured by dev/baseline/capture_baseline.R
# (run that script once before enabling these).
#
# Gated on SPATIALGE_GOLDEN=true: golden/ is Rbuildignored/gitignored and
# never ships, so these must not run by default (CRAN, a fresh clone without
# a captured baseline, etc. would otherwise see unconditional failures).
#

.golden_dir <- testthat::test_path("..", "..", "tests", "golden")

.run_golden_comparison <- function(cores){
  source(testthat::test_path("..", "..", "dev", "baseline", "pipeline.R"), local = TRUE)
  source(testthat::test_path("..", "..", "dev", "baseline", "compare_baseline.R"), local = TRUE)

  fresh <- bl_pipeline(cores = cores, verbose = FALSE)
  for(name in names(fresh)){
    fp <- file.path(.golden_dir, sprintf("%s__cores%d.rds", name, cores))
    if(!file.exists(fp)){
      testthat::fail(sprintf("no golden file captured for '%s' (cores=%d) -- re-run dev/baseline/capture_baseline.R", name, cores))
      next
    }
    golden <- readRDS(fp)
    diff <- bl_compare(golden, fresh[[name]])
    testthat::expect_null(diff, label = sprintf("%s (cores=%d) vs captured golden", name, cores))
  }
}

testthat::test_that("golden baseline: fresh cores=1 output matches captured baseline", {
  testthat::skip_if_not(Sys.getenv("SPATIALGE_GOLDEN") == "true", "golden comparison disabled (set SPATIALGE_GOLDEN=true to enable)")
  testthat::skip_if_not(dir.exists(.golden_dir), "no tests/golden/ captured -- run dev/baseline/capture_baseline.R first")

  .run_golden_comparison(cores = 1)
})

testthat::test_that("golden baseline: fresh cores=3 output matches captured baseline", {
  testthat::skip_if_not(Sys.getenv("SPATIALGE_GOLDEN") == "true", "golden comparison disabled (set SPATIALGE_GOLDEN=true to enable)")
  testthat::skip_if_not(dir.exists(.golden_dir), "no tests/golden/ captured -- run dev/baseline/capture_baseline.R first")

  .run_golden_comparison(cores = 3)
})

testthat::test_that("captured baseline is invariant to cores=1 vs cores=3", {
  testthat::skip_if_not(Sys.getenv("SPATIALGE_GOLDEN") == "true", "golden comparison disabled (set SPATIALGE_GOLDEN=true to enable)")
  testthat::skip_if_not(dir.exists(.golden_dir), "no tests/golden/ captured -- run dev/baseline/capture_baseline.R first")

  source(testthat::test_path("..", "..", "dev", "baseline", "compare_baseline.R"), local = TRUE)

  names1 <- list.files(.golden_dir, pattern = "__cores1\\.rds$")
  fn_names <- sub("__cores1\\.rds$", "", names1)
  testthat::expect_gt(length(fn_names), 0)

  for(name in fn_names){
    f1 <- file.path(.golden_dir, sprintf("%s__cores1.rds", name))
    f3 <- file.path(.golden_dir, sprintf("%s__cores3.rds", name))
    if(!file.exists(f3)){
      testthat::fail(sprintf("no cores3 golden file for '%s'", name))
      next
    }
    diff <- bl_compare(readRDS(f1), readRDS(f3))
    testthat::expect_null(diff, label = sprintf("%s: cores1 vs cores3", name))
  }
})
