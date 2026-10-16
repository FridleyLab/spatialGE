##
# Error-catalog consistency.
#
# Every user-facing validation error in the package is supposed to be raised
# through raise_err(), so that it (a) carries the same "spatialGE exception:"
# signature a caller can pattern-match on, and (b) has its wording stored in
# inst/err.csv rather than inline at the call site. That was not true before:
# ~30 sites used a bare stop(), five of them for the identical condition
# ("Could not recognize number of cores requested"), and three sites hardcoded a
# "spatialGE exception:" string by hand instead of using the catalog.
#
# These tests cover the shared machinery and pin a representative error from each
# family, so a new bare stop() in one of them is caught.

.ec_fixture <- function(transformed = TRUE) {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- suppressWarnings(suppressMessages(
    STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
           cores = 1, verbose = FALSE)
  ))
  if (transformed) {
    st <- suppressWarnings(suppressMessages(
      transform_data(st, method = "log", cores = 1)
    ))
  }
  st
}

# The error a call raises, or NA if it does not raise.
.err_of <- function(expr) {
  tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
}

testthat::test_that("every code in err.csv is unique and resolvable", {
  db <- suppressWarnings(utils::read.csv(
    system.file("err.csv", package = "spatialGE"),
    header = FALSE, quote = "'"))

  testthat::expect_gt(nrow(db), 50)
  ## Duplicate codes would make raise_err() return a length-2 vector, which the
  ## guard turns into an "unknown error code" report -- so catch them here.
  testthat::expect_identical(anyDuplicated(db[[1]]), 0L)
  testthat::expect_false(any(is.na(db[[2]])))
  ## read.csv(quote = "'") means an apostrophe inside a message would terminate
  ## the field early and silently corrupt that entry.
  testthat::expect_false(any(grepl("'", db[[2]], fixed = TRUE)))
})

testthat::test_that("raise_err() names an unknown code instead of failing obscurely", {
  ## Previously the lookup returned character(0) and stop(character(0)) failed
  ## with "argument is of length zero", naming neither the code nor the cause.
  msg <- .err_of(raise_err(err_code = "error9999"))
  testthat::expect_match(msg, "unknown error code error9999")

  msg_null <- .err_of(raise_err())
  testthat::expect_match(msg_null, "unknown error code <NULL>")
})

testthat::test_that("raise_err() interpolates formatted messages", {
  msg <- .err_of(raise_err(err_code = "error0073", samplename = NULL,
                           "GeneA, GeneB"))
  testthat::expect_match(msg, "spatialGE exception")
  testthat::expect_match(msg, "GeneA, GeneB")

  ## samplename is appended after the formatted text.
  msg2 <- .err_of(raise_err(err_code = "error0070", samplename = "alpha"))
  testthat::expect_match(msg2, "Sample: alpha\\.")
})

testthat::test_that("user-facing validation errors carry the spatialGE signature", {
  st <- .ec_fixture()

  ## One representative site from each family that used a bare stop() before.
  sites <- list(
    "SThet(genes=NULL)"          = function() SThet(st, genes = NULL),
    "STdiff(pairwise bad)"       = function() STdiff(st, annot = NULL),
    "STclust(w out of range)"    = function() STclust(st, ws = 5),
    "gene_interpolation(NULL)"   = function() gene_interpolation(st, genes = NULL)
  )

  for (nm in names(sites)) {
    msg <- .err_of(sites[[nm]]())
    testthat::expect_false(is.na(msg), label = paste(nm, "did not error"))
    testthat::expect_match(msg, "^spatialGE exception:", label = nm)
  }
})

testthat::test_that("resolve_cores() rejects an unusable core count in one place", {
  ## Previously the caller checked, and did so inconsistently: some sites raised
  ## error0032/error0036, five used a bare stop(), and transform_data()'s
  ## gene_meta path did not check at all, passing NA straight to sge_lapply().
  for (bad in list("abc", NA, 0, -1)) {
    msg <- .err_of(resolve_cores(10, bad))
    testthat::expect_match(msg, "spatialGE exception")
    testthat::expect_match(msg, "number of cores")
  }

  ## Valid requests must still resolve.
  testthat::expect_equal(resolve_cores(10, cores = 2), 2L)
})

testthat::test_that("resolve_cores() preserves the caller's coercion", {
  ## Fidelity rule: coerce='ceiling' and 'none' deliberately keep a non-integral
  ## value, because consumers that are not mclapply depend on it surviving. The
  ## centralised guard must validate without coercing.
  testthat::expect_identical(resolve_cores(10, cores = 2.5, coerce = "ceiling"), 3)
  testthat::expect_identical(resolve_cores(10, cores = 2.5, coerce = "none"), 2.5)
})

testthat::test_that("a bad cores= value fails cleanly at a call site, not as NA", {
  st <- .ec_fixture(transformed = FALSE)

  ## transform_data()'s gene_meta path was the unguarded one.
  msg <- .err_of(transform_data(st, method = "log", cores = "abc"))
  testthat::expect_match(msg, "^spatialGE exception:.*number of cores")
})
