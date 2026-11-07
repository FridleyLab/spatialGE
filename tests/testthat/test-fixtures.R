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

testthat::test_that("list-of-data.frames + xlsx samples file builds a valid STlist", {
  f <- fixtures$list_dfs_samplesfile_xlsx
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, samples = f$samples, cores = 1)
  testthat::expect_s4_class(st, "STlist")
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

testthat::test_that("samplesfile_matrices (paths-only samples file): unsupported input", {
  # process_sample_names_from_file(), the function that would have resolved
  # counts_path/coords_path columns from a paths-only samples file, was dead
  # code (never called from anywhere) and was removed in Phase 4 -- so this
  # input format (samples file with no rnacounts/spotcoords) has never
  # actually been supported. STlist() now fails with a clear validation error
  # (rnacounts is required) instead of the internal "argument is of length
  # zero" it previously raised due to a separate cores-resolution bug.
  f <- fixtures$samplesfile_matrices
  testthat::expect_error(STlist(samples = f$samples), "rnacounts argument must be provided")
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
