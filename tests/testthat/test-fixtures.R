##
# STlist() construction against the synthetic fixtures in
# inst/testdata/make_fixtures.R (see helper-fixtures.R for how `fixtures`
# gets built). These are all offline/seeded/deterministic -- no skip needed.
#

testthat::test_that("list-of-data.frames input builds a valid STlist", {
  f <- fixtures$list_dfs
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, cores = 1)
  testthat::expect_s4_class(st, "STlist")
  testthat::expect_setequal(names(st@counts), fixtures$base$samples)
})

testthat::test_that("list-of-data.frames + csv samples file builds a valid STlist", {
  f <- fixtures$list_dfs_samplesfile_csv
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("list-of-data.frames + xlsx samples file: known detect_input_samples() bug", {
  # detect_input_samples()'s xlsx branch matches on the regex "//.xlsx$" (two
  # slashes), which can never match a real path ending in ".xlsx" -- so an
  # xlsx samples file always falls through to comma/tab content-sniffing on
  # its raw (binary) bytes and fails. Tracked for a fix in Phase 11; this test
  # documents the current, expected failure rather than asserting success.
  f <- fixtures$list_dfs_samplesfile_xlsx
  testthat::expect_error(
    STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples),
    "not comma or tab-delimited"
  )
})

testthat::test_that("comma-delimited matrix files build a valid STlist", {
  f <- fixtures$delim_csv
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("tab-delimited matrix files build a valid STlist", {
  f <- fixtures$delim_tsv
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("samplesfile_matrices (paths-only samples file): known STlist() bugs", {
  # Two independent bugs currently make this input unusable:
  #  1. STlist()'s cores-resolution `if(is.null(cores) && !is.null(rnacounts))`
  #     is FALSE whenever rnacounts is NULL (as it must be for this format),
  #     so it falls to `cores = as.integer(cores)` on a NULL cores, giving
  #     integer(0); `is.na(integer(0))` in the next `if()` then errors with
  #     "argument is of length zero" before samples-file-specific logic runs.
  #  2. process_sample_names_from_file(), the function that would resolve the
  #     counts_path/coords_path columns from the samples file, is defined in
  #     R/STList.R but never called from anywhere -- this input path is dead
  #     code regardless of bug #1. Tracked for Phase 11.
  f <- fixtures$samplesfile_matrices
  testthat::expect_error(STlist(samples = f$samples), "argument is of length zero")
})

testthat::test_that("Visium MEX directories build a valid STlist", {
  f <- fixtures$visium_mex
  st <- STlist(rnacounts = f$rnacounts, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("CosMx exprMat/metadata csv build a valid STlist", {
  f <- fixtures$cosmx
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("Visium 10x h5 directories build a valid STlist", {
  f <- fixtures$visium_h5
  st <- STlist(rnacounts = f$rnacounts, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("Xenium 10x h5 directories build a valid STlist", {
  f <- fixtures$xenium_h5
  st <- STlist(rnacounts = f$rnacounts, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})

testthat::test_that("Xenium MEX directories build a valid STlist", {
  f <- fixtures$xenium_mex
  st <- STlist(rnacounts = f$rnacounts, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
})
