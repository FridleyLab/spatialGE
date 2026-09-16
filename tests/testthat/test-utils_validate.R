##
# Unit tests for the R/utils_validate.R helpers (assert_stlist, assert_transformed,
# resolve_samples) introduced to de-duplicate validation logic scattered across
# ~19 files. These tests exercise the helpers directly against small, purpose-built
# STlist stand-ins (rather than the full inst/testdata/ fixtures) so that every
# parameter combination used at a real call site -- including known, intentionally
# preserved bugs -- can be pinned down exactly. See plans/this-is-an-r-velvety-music.md
# Phase 6 for the call-site catalogue these tests are checked against.
#

make_stlist_stub <- function(spatial_meta = list(), tr_counts = list(), counts = list()) {
  methods::new("STlist",
               counts = counts,
               spatial_meta = spatial_meta,
               gene_meta = list(),
               sample_meta = tibble::tibble(),
               tr_counts = tr_counts,
               gene_krige = list(),
               misc = list())
}


# assert_stlist -----------------------------------------------------------------

testthat::test_that("assert_stlist passes silently for a valid STlist", {
  st <- make_stlist_stub()
  testthat::expect_no_error(assert_stlist(st))
})

testthat::test_that("assert_stlist default message/check_null+check_class matches original sites", {
  testthat::expect_error(assert_stlist(NULL), "The input must be a STlist.", fixed = TRUE)
  testthat::expect_error(assert_stlist(list()), "The input must be a STlist.", fixed = TRUE)
})

testthat::test_that("assert_stlist honors a custom message (compare_SThet.R / summarize_STlist.R sites)", {
  testthat::expect_error(
    assert_stlist(NULL, message = "The input must be a STList."),
    "The input must be a STList.",
    fixed = TRUE
  )
})

testthat::test_that("assert_stlist(check_null=FALSE) reproduces filter_data.R's A10 gap: a non-null, non-STlist input is NOT rejected", {
  # filter_data.R historically only checked is.null(x) and never checked the
  # class -- a pre-existing gap, preserved verbatim rather than fixed here.
  # Tracked for a fix in Phase 11.
  testthat::expect_no_error(assert_stlist(list(), message = "Please, specify an STList to filter.", check_class = FALSE))
  testthat::expect_error(
    assert_stlist(NULL, message = "Please, specify an STList to filter.", check_class = FALSE),
    "Please, specify an STList to filter.",
    fixed = TRUE
  )
})

testthat::test_that("assert_stlist(check_null=FALSE) reproduces summarize_STlist.R's A4 site: only the class is checked", {
  st <- make_stlist_stub()
  testthat::expect_no_error(assert_stlist(st, message = "The input is not an STlist", check_null = FALSE))
  testthat::expect_error(
    assert_stlist(list(), message = "The input is not an STlist", check_null = FALSE),
    "The input is not an STlist",
    fixed = TRUE
  )
})

testthat::test_that("assert_stlist dispatches to raise_err() when err_code is supplied (pseudobulk_samples.R site)", {
  testthat::expect_error(
    assert_stlist(NULL, err_code = "error0018"),
    "spatialGE exception: pseudobulk_samples: The input must be an STList.",
    fixed = TRUE
  )
})


# assert_transformed -------------------------------------------------------------

testthat::test_that("assert_transformed passes silently once tr_counts is populated", {
  st <- make_stlist_stub(tr_counts = list(a = 1))
  testthat::expect_no_error(assert_transformed(st))
})

testthat::test_that("assert_transformed default message matches original sites", {
  st <- make_stlist_stub()
  testthat::expect_error(assert_transformed(st), "There are no transformed counts in this STlist.", fixed = TRUE)
})

testthat::test_that("assert_transformed honors a custom message (gene_interpolation.R / plot_counts.R sites)", {
  st <- make_stlist_stub()
  testthat::expect_error(
    assert_transformed(st, message = "There are not transformed counts in this STList."),
    "There are not transformed counts in this STList.",
    fixed = TRUE
  )
})

testthat::test_that("assert_transformed dispatches to raise_err() when err_code is supplied (STclust.R site)", {
  st <- make_stlist_stub()
  testthat::expect_error(
    assert_transformed(st, err_code = "error0007"),
    "spatialGE exception: STclust: Could not find normalized counts. Please use transform data before using STclust.",
    fixed = TRUE
  )
})


# resolve_samples -----------------------------------------------------------------

testthat::test_that("resolve_samples returns the slot's names, in order, when samples is NULL", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  testthat::expect_identical(resolve_samples(st, NULL, slot = "spatial_meta"), c("a", "b", "c"))
})

testthat::test_that("resolve_samples(numeric_handling='na_omit') drops out-of-range indices (STclust.R/STdiff.R/STenrich.R/STgradient.R sites)", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  out <- resolve_samples(st, c(1, 5), slot = "spatial_meta", numeric_handling = "na_omit")
  testthat::expect_identical(out, "a")
})

testthat::test_that("resolve_samples(numeric_handling='plain') keeps a resulting NA for an out-of-range index (plot_counts.R/distribution_plots.R/plot_image.R/... sites)", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  out <- resolve_samples(st, c(1, 5), slot = "spatial_meta", numeric_handling = "plain")
  testthat::expect_identical(out, c("a", NA_character_))
})

testthat::test_that("resolve_samples(char_transform='exact') filters to only the names present, in the caller's order (STclust.R/STdiff.R/STenrich.R sites)", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  out <- resolve_samples(st, c("c", "a", "zzz"), slot = "spatial_meta", char_transform = "exact", post_check = "none")
  testthat::expect_identical(out, c("c", "a"))
})

testthat::test_that("resolve_samples(char_transform='none') never filters or validates character samples (SThet.R/gene_interpolation.R/STplot_interpolation.R/... sites)", {
  # These sites never checked a character `samples` argument against the
  # available names at all -- a pre-existing gap, preserved verbatim. Tracked
  # for a fix in Phase 11.
  st <- make_stlist_stub(tr_counts = list(a = 1, b = 2))
  out <- resolve_samples(st, c("a", "zzz"), slot = "tr_counts", char_transform = "none", post_check = "none")
  testthat::expect_identical(out, c("a", "zzz"))
})

testthat::test_that("resolve_samples(char_transform='validate_all') stops without filtering when a name is missing (filter_data.R's C8 site)", {
  st <- make_stlist_stub(counts = list(a = 1, b = 2))
  testthat::expect_identical(
    resolve_samples(st, c("a", "b"), slot = "counts", numeric_handling = "plain",
                    char_transform = "validate_all", post_check = "none"),
    c("a", "b")
  )
  testthat::expect_error(
    resolve_samples(st, c("a", "zzz"), slot = "counts", numeric_handling = "plain",
                    char_transform = "validate_all", post_check = "none",
                    message = "Samples in rm_tissue are not present in this STlist."),
    "Samples in rm_tissue are not present in this STlist.",
    fixed = TRUE
  )
})

testthat::test_that("resolve_samples(post_check='any_in') stops when no requested sample matches (STclust.R/STdiff.R sites)", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  testthat::expect_error(
    resolve_samples(st, "zzz", slot = "spatial_meta", numeric_handling = "na_omit",
                    char_transform = "exact", post_check = "any_in",
                    message = "None of the requested samples are present in the STlist."),
    "None of the requested samples are present in the STlist.",
    fixed = TRUE
  )
})

testthat::test_that("resolve_samples(post_check='any_in') dispatches to raise_err() when err_code is supplied (STenrich.R/STgradient.R sites)", {
  st <- make_stlist_stub(spatial_meta = list(a = 1, b = 2, c = 3))
  testthat::expect_error(
    resolve_samples(st, "zzz", slot = "spatial_meta", numeric_handling = "na_omit",
                    char_transform = "exact", post_check = "any_in", err_code = "error0020"),
    "spatialGE exception: STenrich: None of the requested samples are present in the STlist.",
    fixed = TRUE
  )
})

testthat::test_that("resolve_samples(post_check='grep') known bug: a substring match falsely validates a nonexistent sample (plot_spatial_geneset.R/plot_spatial_meta.R/plot_spatial_expression.R/plot_image.R sites)", {
  # These four sites check sample existence with
  # `grep(paste0(samples, collapse='|'), available)`, an unanchored
  # substring/regex existence check rather than an exact %in% match. It also
  # never filters `samples` -- it only checks that *some* loose match exists,
  # then returns the caller's `samples` unchanged. Here "S1" is not itself an
  # available sample name, but it is a substring of both "S10" and "S11", so
  # the check passes and "S1" is returned as if it were valid. Tracked for a
  # fix in Phase 11.
  st <- make_stlist_stub(spatial_meta = list(S10 = 1, S11 = 2))
  out <- resolve_samples(st, "S1", slot = "spatial_meta", numeric_handling = "plain",
                         char_transform = "none", post_check = "grep",
                         message = "The requested samples are not present in the STlist.")
  testthat::expect_identical(out, "S1")
})

testthat::test_that("resolve_samples(post_check='grep') still rejects a sample with no substring match at all", {
  st <- make_stlist_stub(spatial_meta = list(S10 = 1, S11 = 2))
  testthat::expect_error(
    resolve_samples(st, "zzz", slot = "spatial_meta", numeric_handling = "plain",
                    char_transform = "none", post_check = "grep",
                    message = "The requested samples are not present in the STlist."),
    "The requested samples are not present in the STlist.",
    fixed = TRUE
  )
})
