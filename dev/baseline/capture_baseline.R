## Captures the current, pre-refactor behavior of every exported function
## into tests/golden/<name>__cores<N>.rds (gitignored/Rbuildignored -- these
## never leave the local machine). Run this once now, before any R/ change
## in Phases 4+, and re-run it after each phase to diff against with
## dev/baseline/compare_baseline.R / tests/testthat/test-golden.R.
##
## Usage (run from the package root):
##   conda run -n spatialge-dev Rscript dev/baseline/capture_baseline.R

devtools::load_all(".", quiet = TRUE)
source("inst/testdata/make_fixtures.R")
source("dev/baseline/pipeline.R")

golden_dir <- "tests/golden"
dir.create(golden_dir, showWarnings = FALSE, recursive = TRUE)

for(cores in c(1, 3)){
  message(sprintf("== capturing baseline at cores=%d ==", cores))
  results <- bl_pipeline(cores, verbose = FALSE)
  for(name in names(results)){
    fp <- file.path(golden_dir, sprintf("%s__cores%d.rds", name, cores))
    saveRDS(results[[name]], fp, version = 3)
    ok <- results[[name]]$ok
    message(sprintf("  %-24s %s", name, if(ok) "captured" else paste("captured error:", results[[name]]$value)))
  }
}

message("Done. Golden files written to ", golden_dir)
