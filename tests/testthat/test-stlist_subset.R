##
# STlist subsetting and length.
#
# Before this, an STlist could not be subset at all: `st[1]` failed with
# "object of type 'S4' is not subsettable", so there was no way to pull one
# sample out of a multi-sample object and carry on working with it -- the only
# route was to rebuild an STlist from the raw files. `length()` did not exist
# either, so seq_along()/seq_len(length(x)) idioms silently used length 0.
#
# Subsetting has to keep four sample-keyed slots in step with @counts
# (@spatial_meta, @gene_meta, @tr_counts and the per-sample level of
# @gene_krige) plus align @sample_meta rows. Getting that wrong is exactly the
# kind of thing setValidity() now catches, so every result here is also checked
# with validObject().

.ss_fixture <- function(transformed = FALSE) {
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

testthat::test_that("length() reports the number of spatial arrays", {
  st <- .ss_fixture()

  testthat::expect_identical(length(st), length(names(st@counts)))
  ## The idiom this enables.
  testthat::expect_length(seq_len(length(st)), length(names(st@counts)))
})

testthat::test_that("x[i] subsets samples by name, position and logical", {
  st <- .ss_fixture()
  sn <- names(st@counts)
  testthat::expect_gt(length(sn), 1)

  by_name <- st[sn[2]]
  testthat::expect_s4_class(by_name, "STlist")
  testthat::expect_identical(names(by_name@counts), sn[2])
  testthat::expect_true(validObject(by_name))

  ## A single-sample result is now genuinely matrix-like.
  testthat::expect_identical(dim(by_name), c(nrow(st@counts[[sn[2]]]), ncol(st@counts[[sn[2]]])))

  testthat::expect_identical(names(st[c(sn[3], sn[1])]@counts), c(sn[3], sn[1]))
  testthat::expect_identical(names(st[seq_len(2)]@counts), sn[seq_len(2)])
  testthat::expect_identical(names(st[-1]@counts), sn[-1])

  logical_i <- rep(FALSE, length(sn)); logical_i[2] <- TRUE
  testthat::expect_identical(names(st[logical_i]@counts), sn[2])
})

testthat::test_that("subsetting keeps every sample-keyed slot in step", {
  st <- .ss_fixture(transformed = TRUE)
  sn <- names(st@counts)

  sub <- st[c(sn[3], sn[1])]

  testthat::expect_true(validObject(sub))
  for (slot_name in c("counts", "spatial_meta", "gene_meta", "tr_counts")) {
    testthat::expect_identical(names(methods::slot(sub, slot_name)), c(sn[3], sn[1]),
                               label = paste0("@", slot_name))
  }
  ## Values must be untouched, not merely renamed.
  testthat::expect_equal(sub@tr_counts[[1]], st@tr_counts[[sn[3]]], tolerance = 0)
  testthat::expect_identical(nrow(sub@gene_meta[[1]]), nrow(st@gene_meta[[sn[3]]]))
})

testthat::test_that("subsetting aligns @sample_meta with the selected samples", {
  st <- .ss_fixture()
  sn <- names(st@counts)

  ## The identifying column is not fixed: STlist() uses 'sample_name' for list
  ## input, while a clinical CSV keeps its own header (the Thrane data uses
  ## 'sampleID'). Match on the first column, which is what show() treats as the
  ## sample identifier.
  testthat::expect_identical(nrow(st@sample_meta), length(sn))

  sub <- st[c(sn[3], sn[1])]

  testthat::expect_identical(nrow(sub@sample_meta), 2L)
  testthat::expect_identical(as.character(sub@sample_meta[[1]]), c(sn[3], sn[1]))
})

testthat::test_that("an unknown sample name is an error that lists the alternatives", {
  st <- .ss_fixture()
  sn <- names(st@counts)

  testthat::expect_error(
    st["NoSuchSample"],
    "spatialGE exception.*unknown sample name"
  )
  ## And it names what is available, rather than just refusing.
  msg <- tryCatch(st["NoSuchSample"], error = function(e) conditionMessage(e))
  testthat::expect_match(msg, sn[1], fixed = TRUE)
})

testthat::test_that("invalid selections are rejected rather than silently recycled", {
  st <- .ss_fixture()

  ## A short logical vector would recycle and select the wrong samples.
  testthat::expect_error(st[c(TRUE, FALSE)], "spatialGE exception")
  ## Out of range must not silently return a shorter object.
  testthat::expect_error(st[99], "spatialGE exception")
  ## Selecting nothing is not a useful STlist.
  testthat::expect_error(st[rep(FALSE, length(st))], "spatialGE exception")
  ## An unsupported index type.
  testthat::expect_error(st[list(1)], "spatialGE exception")
})

testthat::test_that("two-dimensional subsetting explains why it is unsupported", {
  st <- .ss_fixture()

  ## An STlist holds several arrays with independent gene/spot sets, so there is
  ## no shared second dimension. Say that instead of failing generically.
  testthat::expect_error(st[1, 1], "one-dimensional")
  ## [[ is not a sample accessor either.
  testthat::expect_error(st[[1]], "spatialGE exception")
})

testthat::test_that("validObject() accepts an object with kriging results", {
  ## Regression: setValidity() originally listed @gene_krige among the
  ## sample-keyed slots and compared its names to the sample names. But
  ## @gene_krige is keyed by GENE first (x@gene_krige[[gene]][[sample]]), so
  ## every kriged object was reported invalid. It went unnoticed because slot
  ## assignment via `@<-` does not run validity, so gene_interpolation()
  ## succeeded and only a later validObject() (or new()) failed.
  st <- .ss_fixture(transformed = TRUE)
  sn <- names(st@counts)

  st@gene_krige <- list(
    GeneA = stats::setNames(lapply(sn, function(s) list(success_or_not = "ok")), sn)
  )

  testthat::expect_true(validObject(st))
})

testthat::test_that("subsetting filters kriging results to the kept samples", {
  st <- .ss_fixture(transformed = TRUE)
  sn <- names(st@counts)

  st@gene_krige <- list(
    GeneA = stats::setNames(lapply(sn, function(s) list(krige_out = s)), sn),
    GeneB = stats::setNames(lapply(sn, function(s) list(krige_out = s)), sn)
  )

  sub <- st[sn[1]]

  testthat::expect_true(validObject(sub))
  ## Genes survive (kriging is per gene) ...
  testthat::expect_setequal(names(sub@gene_krige), c("GeneA", "GeneB"))
  ## ... but only the kept sample remains within each.
  for (g in names(sub@gene_krige)) {
    testthat::expect_identical(names(sub@gene_krige[[g]]), sn[1])
  }
  testthat::expect_identical(sub@gene_krige$GeneA[[sn[1]]]$krige_out, sn[1])
})
