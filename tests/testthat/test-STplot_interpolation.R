##
# STplot_interpolation(): the default (no-argument) call, and its kriging
# prerequisite.
#
# Regression coverage for a defect that made the documented default unusable:
# `stop("Please, enter one or more genes to plot.")` ran *before* the branch that
# resolves `genes='top'`, so calling the function without `genes` -- exactly what
# ?STplot_interpolation shows, and what its own @param documents as "If 'top'"
# -- always failed. The gitignored golden tier only ever called it with explicit
# gene names, so nothing caught it.
#
# The synthetic fixtures make these offline and deterministic.

.sti_fixture <- function() {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  transform_data(st, method = "log", cores = 1)
}

.sti_kriged <- function(top_n = 3) {
  st <- .sti_fixture()
  ## Kriging on the small fixture emits variogram-fit warnings; this file tests
  ## STplot_interpolation(), not gene_interpolation().
  suppressWarnings(
    suppressMessages(gene_interpolation(st, genes = "top", top_n = top_n,
                                        cores = 1, verbose = FALSE))
  )
}

testthat::test_that("STplot_interpolation() errors clearly before kriging has run", {
  st <- .sti_fixture()

  ## Used to fail later on an empty color scale, or produce +/-Inf.
  testthat::expect_error(
    STplot_interpolation(st),
    "spatialGE exception.*STplot_interpolation.*gene_interpolation"
  )
})

testthat::test_that("STplot_interpolation() works with no `genes` argument", {
  ## The headline regression: this call used to always error.
  st <- .sti_kriged()

  plots <- STplot_interpolation(st)

  testthat::expect_true(length(plots) > 0)
  testthat::expect_s3_class(plots[[1]], "ggplot")
})

testthat::test_that("genes=NULL and genes='top' select the same plots", {
  st <- .sti_kriged()

  default_plots <- STplot_interpolation(st, top_n = 3)
  top_plots <- STplot_interpolation(st, genes = "top", top_n = 3)

  testthat::expect_equal(length(default_plots), length(top_plots))
})

testthat::test_that("STplot_interpolation() still accepts explicit gene names", {
  st <- .sti_kriged()
  genes <- utils::head(names(st@gene_krige), 2)

  plots <- STplot_interpolation(st, genes = genes)

  testthat::expect_true(length(plots) > 0)
})

testthat::test_that("STplot_interpolation() warns and plots the rest when a named gene was not kriged", {
  st <- .sti_kriged()
  present <- utils::head(names(st@gene_krige), 1)

  testthat::expect_warning(
    plots <- STplot_interpolation(st, genes = c(present, "NoSuchKrigedGene")),
    "No kriging results for: NoSuchKrigedGene"
  )
  testthat::expect_true(length(plots) > 0)
})

testthat::test_that("STplot_interpolation() errors clearly when no named gene was kriged", {
  st <- .sti_kriged()

  testthat::expect_error(
    STplot_interpolation(st, genes = "NoSuchKrigedGene"),
    "spatialGE exception.*STplot_interpolation.*gene_interpolation"
  )
})

testthat::test_that("a default 'top' request does not warn about genes that were never kriged", {
  ## top_n here is independent of the top_n used at kriging time, so asking for
  ## more genes than were interpolated is expected rather than a user error.
  st <- .sti_kriged(top_n = 2)

  testthat::expect_no_warning(STplot_interpolation(st, top_n = 50))
})
