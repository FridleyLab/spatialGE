## Phase 9 coverage harness: an objective answer to whether the fixtures and
## test suite reach every branch in R/, instead of relying on "the tests
## pass" as a stand-in for "the tests exercise the code".
##
## covr::package_coverage(type="none", ...) builds and installs an
## instrumented copy of the package itself, so the `code` argument below
## just needs to exercise it -- it does not need devtools::load_all() or any
## of the covr::load_all()-side helpers used elsewhere in dev/.
##
## Two code paths are combined to get broad coverage without needing a
## captured tests/golden/ baseline (which is local-only and gitignored):
##   1. dev/baseline/pipeline.R's bl_pipeline(), run at both cores=1 and
##      cores=3 (mirrors test-golden.R's two-cores comparison, so both the
##      serial and the mc.cores>1 branch of sge_lapply() get exercised) --
##      this threads one STlist through every exported function in R/, in
##      dependency order.
##   2. testthat::test_dir("tests/testthat"), which covers what bl_pipeline()
##      doesn't reach: alternate STlist constructors (list-of-data.frames,
##      Visium/Xenium/CosMx readers, samples-file variants), and the
##      utils_validate.R helpers exercised directly.
##
## Writes a diffable plain-text summary (per-file percentage, then an
## explicit zero-coverage line listing) to dev/profiling/out/coverage.txt,
## consistent with rprof_summary.R's plain-text-over-HTML choice.
##
## Run from the package root:
##   Rscript dev/profiling/coverage.R

out_dir <- "dev/profiling/out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cov <- covr::package_coverage(
  path = ".",
  type = "none",
  code = c(
    'source("inst/testdata/make_fixtures.R"); source("dev/baseline/pipeline.R"); bl_pipeline(cores = 1, verbose = FALSE); bl_pipeline(cores = 3, verbose = FALSE)',
    'testthat::test_dir("tests/testthat", stop_on_failure = FALSE)'
  )
)

out_path <- file.path(out_dir, "coverage.txt")
con <- file(out_path, open = "w")

by_file <- covr::coverage_to_list(cov)$filecoverage
by_file <- by_file[order(by_file)]

writeLines(sprintf("overall: %.2f%%", covr::percent_coverage(cov)), con)
writeLines("", con)
writeLines("by file (ascending, worst first):", con)
for(fn in names(by_file)){
  writeLines(sprintf("  %6.2f%%  %s", by_file[[fn]], fn), con)
}

writeLines("", con)
writeLines("zero-coverage lines:", con)
zero <- covr::zero_coverage(cov)
if(nrow(zero) == 0){
  writeLines("  (none)", con)
} else{
  zero <- zero[order(zero$filename, zero$line), ]
  for(i in seq_len(nrow(zero))){
    writeLines(sprintf("  %s:%d  %s()", zero$filename[i], zero$line[i], zero$functions[i]), con)
  }
}

close(con)
message(sprintf("[coverage] wrote summary to %s", out_path))
