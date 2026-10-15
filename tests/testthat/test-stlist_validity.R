##
# STlist class contract: dim() and validity.
#
# `dim()` used to return a *named list* of per-sample dimensions, which violates
# the base R contract (dim() must return an integer vector or NULL). It made
# nrow()/ncol() meaningless on an STlist and was reported by R CMD check as
# "dim() returns class 'list'". Nothing in the package or its tests called it, so
# nothing caught it.
#
# The class also had no setValidity(), so a malformed STlist could be constructed
# freely and failed later, far from the cause.

.sv_fixture <- function() {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  suppressWarnings(suppressMessages(
    STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
           cores = 1, verbose = FALSE)
  ))
}

testthat::test_that("dim() returns an integer vector for a single-sample STlist", {
  st <- .sv_fixture()
  one <- new("STlist",
             counts = st@counts[1],
             spatial_meta = st@spatial_meta[1],
             gene_meta = list(),
             sample_meta = st@sample_meta,
             tr_counts = list(),
             gene_krige = list(),
             misc = st@misc)

  d <- dim(one)

  testthat::expect_type(d, "integer")
  testthat::expect_length(d, 2)
  testthat::expect_identical(d, c(nrow(st@counts[[1]]), ncol(st@counts[[1]])))
})

testthat::test_that("dim() returns NULL (not a list) for a multi-sample STlist", {
  ## Regression: this used to be a named list, e.g. list(alpha = c(G, n), ...).
  st <- .sv_fixture()
  testthat::expect_gt(length(st@counts), 1)

  d <- dim(st)

  testthat::expect_null(d)
  ## The base contract, stated directly: nrow()/ncol() must not return a list.
  testthat::expect_false(is.list(nrow(st)))
  testthat::expect_null(nrow(st))
})

testthat::test_that("a normally-built STlist is valid", {
  st <- .sv_fixture()

  testthat::expect_true(validObject(st))
})

testthat::test_that("an STlist stays valid through transform_data()", {
  ## Zero-count spots are dropped from both @counts and @spatial_meta, so the
  ## two must remain in step.
  st <- .sv_fixture()
  suppressWarnings(suppressMessages(
    st <- transform_data(st, method = "log", cores = 1)
  ))

  testthat::expect_true(validObject(st))
})

testthat::test_that("setValidity rejects a spatial_meta/counts size mismatch", {
  st <- .sv_fixture()
  s <- names(st@counts)[1]
  st@spatial_meta[[s]] <- st@spatial_meta[[s]][-1, , drop = FALSE]

  testthat::expect_error(validObject(st), "rows but the count matrix has")
})

testthat::test_that("setValidity rejects a sample-keyed slot for an unknown sample", {
  st <- .sv_fixture()
  st@gene_meta <- list(ghost = data.frame(gene = "x"))

  testthat::expect_error(validObject(st), "samples not in @counts")
})

testthat::test_that("setValidity rejects a count matrix without spot names", {
  st <- .sv_fixture()
  colnames(st@counts[[1]]) <- NULL

  testthat::expect_error(validObject(st), "spot/cell \\(column\\) names")
})

testthat::test_that("setValidity rejects a spatial_meta that is not a data.frame", {
  st <- .sv_fixture()
  s <- names(st@counts)[1]
  st@spatial_meta[[s]] <- as.matrix(st@spatial_meta[[s]])

  testthat::expect_error(validObject(st), "must be a data.frame")
})

testthat::test_that("setValidity rejects spatial_meta without libname", {
  st <- .sv_fixture()
  s <- names(st@counts)[1]
  st@spatial_meta[[s]][["libname"]] <- NULL

  testthat::expect_error(validObject(st), "libname")
})
