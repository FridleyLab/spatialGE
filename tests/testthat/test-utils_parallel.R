##
# Tests for the parallel backend, including the non-forking (PSOCK) path added
# so Windows gets parallelism at all.
#
# WHAT THESE CAN AND CANNOT VERIFY. PSOCK works on every platform, so the
# mechanism is exercised here on whatever CI runs -- forced on with
# options(spatialGE.psock=TRUE). What that does NOT verify is Windows itself:
# path handling, the absence of fork, and the shell used by the progress messages
# are all untested until this runs on a Windows runner. Treat a green result here
# as "the non-forking path is correct", not "Windows works".
#

testthat::test_that("resolve_cores consults the override before any platform rule", {
  ## The override used to sit AFTER an unconditional return(1L) for Windows, so
  ## options(spatialGE.cores=) and SPATIALGE_CORES were silently ignored on the
  ## one platform whose parallel behavior differs most.
  withr::local_options(list(spatialGE.cores = 3))
  testthat::expect_equal(resolve_cores(10, cores = NULL), 3L)
  testthat::expect_equal(resolve_cores(10, cores = 8), 3L)
})

testthat::test_that("resolve_cores no longer forces a single core", {
  withr::local_options(list(spatialGE.cores = NULL))
  withr::local_envvar(list(SPATIALGE_CORES = NA))

  ## An explicit request must survive. Previously Windows got 1L regardless,
  ## making the whole package serial there.
  testthat::expect_equal(resolve_cores(10, cores = 4, coerce = "as_integer"), 4L)

  ## coerce='ceiling' deliberately returns a double, and coerce='none' passes the
  ## value through untouched -- fidelity rule 5 makes both load-bearing, since
  ## non-mclapply consumers rely on a non-integral value surviving.
  testthat::expect_identical(resolve_cores(10, cores = 2.5, coerce = "ceiling"), 3)
  testthat::expect_identical(resolve_cores(10, cores = 2.5, coerce = "none"), 2.5)
})

testthat::test_that("sge_use_psock follows the platform but is overridable", {
  withr::local_options(list(spatialGE.psock = NULL))
  testthat::expect_equal(sge_use_psock(), .Platform$OS.type == "windows")

  withr::local_options(list(spatialGE.psock = TRUE))
  testthat::expect_true(sge_use_psock())
})

testthat::test_that("the PSOCK path returns the same results as the serial path", {
  withr::local_options(list(spatialGE.psock = TRUE))

  serial <- lapply(1:6, function(i) i^2)
  parallelised <- sge_lapply(1:6, function(i) i^2, cores = 2)
  testthat::expect_identical(parallelised, serial)

  ## Extra arguments must forward.
  testthat::expect_identical(
    sge_lapply(1:3, function(i, k) i * k, k = 10, cores = 2),
    lapply(1:3, function(i) i * 10))
})

testthat::test_that("PSOCK workers can reach package internals", {
  withr::local_options(list(spatialGE.psock = TRUE))

  ## The failure this guards against is specific and was hit while building it:
  ## a PSOCK worker is a fresh R process, so FUN's references to package
  ## internals fail with "could not find function" unless the worker loads the
  ## package. Detecting a development load needs the .__DEVTOOLS__ namespace
  ## marker -- checking for a DESCRIPTION file reports "installed" for a
  ## dev-loaded package, because its namespace path IS the source directory.
  mats <- list(matrix(c(0, 1, 2, 3, 0, 0), nrow = 2),
               matrix(c(4, 0, 6, 0, 8, 1), nrow = 2))
  res <- sge_lapply(mats, function(m){
    sparse_row_sds(methods::as(Matrix::Matrix(m * 1.0, sparse = TRUE), "CsparseMatrix"))
  }, cores = 2)

  testthat::expect_length(res, 2)
  testthat::expect_true(all(vapply(res, function(v) all(is.finite(v)), logical(1))))
})

testthat::test_that("PSOCK falls back to serial below two cores or two elements", {
  withr::local_options(list(spatialGE.psock = TRUE))

  ## Mirrors mclapply's own short-circuit, so a single-sample STlist does not pay
  ## cluster startup for nothing.
  testthat::expect_identical(sge_lapply(1:4, function(i) i, cores = 1), as.list(1:4))
  testthat::expect_identical(sge_lapply(list(5), function(i) i, cores = 4), list(5))
})

testthat::test_that("the cores override forces serial execution on both backends", {
  ## This is what makes profvis able to see past the parallel layer, so it must
  ## hold for the PSOCK path too and not only for mclapply.
  withr::local_options(list(spatialGE.cores = 1, spatialGE.psock = TRUE))
  testthat::expect_identical(sge_lapply(1:4, function(i) i * 2, cores = 4),
                             lapply(1:4, function(i) i * 2))
})

testthat::test_that("SThet produces identical results through PSOCK and fork", {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  genes <- c("Gene1", "GradGene1", "DomainGeneA")

  ## cores = 2, not 3. R CMD check --as-cran sets _R_CHECK_LIMIT_CORES_, which
  ## makes parallel:::.check_ncores() ERROR above two simultaneous processes.
  ## count_cores() clamps its AUTO-DETECTED value under that variable but
  ## deliberately leaves an explicit cores= alone (fidelity rule 12), so passing 3
  ## here failed the check -- in the test, not the package. Two workers still
  ## exercise real forking and real PSOCK dispatch, which is all this needs.
  forked <- SThet(st, genes = genes, method = "moran", cores = 2, verbose = FALSE)

  psock <- withr::with_options(
    list(spatialGE.psock = TRUE),
    SThet(st, genes = genes, method = "moran", cores = 2, verbose = FALSE))

  ## End-to-end equality on a real analysis, not just on a toy closure: this is
  ## what demonstrates that the sliced payloads serialize correctly, which is the
  ## property the non-forking backend depends on.
  for(s in names(forked@gene_meta)){
    testthat::expect_equal(psock@gene_meta[[s]]$moran_i,
                           forked@gene_meta[[s]]$moran_i,
                           tolerance = 0, label = paste("moran_i via PSOCK", s))
  }
})
