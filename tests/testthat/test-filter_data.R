##
# filter_stage_genes_expr(), STEP 4 of filter_data()'s fixed filtering order.
# Regression test for the bug where the @gene_meta branch indexed with a
# malformed 3-subscript expression and rownames() (tibbles have no real
# rownames, so this returns row-number strings) instead of the `gene`
# column, making that branch a silent no-op. filter_data()'s own STEP 5
# (filter_stage_counts) happens to re-sync @gene_meta from @tr_counts
# afterward regardless, which masks the bug end-to-end -- so this test
# calls filter_stage_genes_expr() directly rather than through filter_data().
#

testthat::test_that("filter_stage_genes_expr removes matching genes from gene_meta, not just counts", {
  f <- fixtures$list_dfs
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, cores = 1)
  st <- transform_data(st, method = "log", cores = 1)

  target_gene <- "GradGene1"
  st_filtered <- filter_stage_genes_expr(st, names(st@counts), paste0("^", target_gene, "$"))

  for (s in names(st_filtered@gene_meta)) {
    testthat::expect_false(target_gene %in% st_filtered@gene_meta[[s]][["gene"]])
    testthat::expect_false(target_gene %in% rownames(st_filtered@counts[[s]]))
    testthat::expect_equal(nrow(st_filtered@gene_meta[[s]]), nrow(st@gene_meta[[s]]) - 1)
  }
})


##
# filter_data() reporting.
#
# Every count threshold defaults to 0 (and the *_max* thresholds to NULL), so a
# call that supplies none of them filters nothing -- which previously looked like
# a successful filter. verbose=TRUE now says so explicitly, and otherwise reports
# what was removed.

.fd_fixture <- function() {
  f <- fixtures$list_dfs
  suppressWarnings(suppressMessages(
    STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords,
           cores = 1, verbose = FALSE)
  ))
}

testthat::test_that("filter_data() says so when a call filters nothing", {
  st <- .fd_fixture()

  testthat::expect_message(
    out <- filter_data(st, verbose = TRUE),
    "nothing was filtered"
  )
  ## And it really is a pass-through.
  testthat::expect_equal(lapply(out@counts, dim), lapply(st@counts, dim))
})

testthat::test_that("filter_data() reports how much it removed", {
  st <- .fd_fixture()
  ## A threshold high enough to drop at least one spot from the fixture.
  cutoff <- stats::quantile(Matrix::colSums(st@counts[[1]]), 0.9)

  testthat::expect_message(
    out <- filter_data(st, spot_minreads = cutoff, verbose = TRUE),
    "removed .* sample\\(s\\), .* spot\\(s\\)/cell\\(s\\) and .* gene\\(s\\)"
  )
  testthat::expect_lt(ncol(out@counts[[1]]), ncol(st@counts[[1]]))
})

testthat::test_that("filter_data(verbose=FALSE) is silent", {
  st <- .fd_fixture()

  testthat::expect_silent(suppressWarnings(filter_data(st, verbose = FALSE)))
})

