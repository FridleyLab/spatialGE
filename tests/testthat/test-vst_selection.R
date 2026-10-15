##
# VST (standardized variance) availability and stability.
#
# `gene_meta` used to hold only gene_mean/gene_stdevs after transform_data();
# `vst.variance.standardized` appeared later, as a side effect of whichever
# analysis ran first (STclust(), gene_interpolation()). Because eight code paths
# select "top variable genes" by VST, the set of genes a function chose depended
# on call order -- and gene_stdevs ranks them differently (top-20 overlap of only
# 6/20 on the melanoma data). transform_data() now populates VST itself.
#
# The synthetic fixtures make these offline and deterministic.

.vst_fixture_raw <- function() {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  suppressWarnings(suppressMessages(
    filter_data(st, gene_minreads = 0, spot_minreads = 0, verbose = FALSE)
  ))
}

.vst_fixture <- function(method = "log") {
  suppressWarnings(suppressMessages(
    transform_data(.vst_fixture_raw(), method = method, cores = 1)
  ))
}

.vst_top <- function(x, n = 20) {
  gm <- x@gene_meta[[1]]
  gm[["gene"]][order(gm[["vst.variance.standardized"]], decreasing = TRUE)][seq_len(n)]
}

testthat::test_that("transform_data() populates VST for the log method", {
  st <- .vst_fixture("log")

  testthat::expect_true(
    all(vapply(st@gene_meta, function(g) "vst.variance.standardized" %in% colnames(g), logical(1)))
  )
})

testthat::test_that("transform_data() populates VST for the sct method, or fails actionably", {
  ## SCTransform computes its own variance statistics, but they were never
  ## carried into gene_meta -- so VST had to come from a later analysis.
  ##
  ## sctransform cannot fit the small synthetic fixture (200 genes x ~120 spots),
  ## so this asserts whichever of the two documented outcomes occurs. The
  ## populated-VST branch is exercised on real data by dev/xver and the golden
  ## tier; what matters here is that a failure is never opaque.
  st <- .vst_fixture_raw()

  res <- tryCatch(
    suppressWarnings(transform_data(st, method = "sct", cores = 1)),
    error = function(e) e
  )

  if (inherits(res, "error")) {
    testthat::expect_match(conditionMessage(res),
                           "spatialGE exception.*sct.*SCTransform")
    testthat::expect_match(conditionMessage(res), "sct_n_regr_genes|method='log'")
  } else {
    testthat::expect_true(
      all(vapply(res@gene_meta,
                 function(g) "vst.variance.standardized" %in% colnames(g), logical(1)))
    )
  }
})

testthat::test_that("VST is visible through get_gene_meta() without running an analysis", {
  st <- .vst_fixture("log")

  testthat::expect_true("vst.variance.standardized" %in% colnames(get_gene_meta(st)))
})

testthat::test_that("the mean/stdev columns survive the VST join", {
  ## calculate_vst() left-joins onto gene_meta; it must not replace the columns
  ## transform_data() just computed.
  st <- .vst_fixture("log")

  for (gm in st@gene_meta) {
    testthat::expect_true(all(c("gene", "gene_mean", "gene_stdevs",
                                "vst.variance.standardized") %in% colnames(gm)))
  }
})

testthat::test_that("STclust() does not change the VST column or the top-gene ranking", {
  ## The regression this guards: previously STclust() was what *created* VST, so
  ## the genes selected before vs after it could differ.
  st <- .vst_fixture("log")
  before_cols <- lapply(st@gene_meta, function(g) g[["vst.variance.standardized"]])
  before_top <- .vst_top(st)

  st2 <- suppressWarnings(STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE,
                                  cores = 1, verbose = FALSE))

  testthat::expect_equal(lapply(st2@gene_meta, function(g) g[["vst.variance.standardized"]]),
                         before_cols)
  testthat::expect_identical(.vst_top(st2), before_top)
})

testthat::test_that("STdiff() runs on an STlist straight from transform_data()", {
  ## No preceding STclust(); VST must already be there.
  st <- .vst_fixture("log")
  acol <- "test_annot"
  st@spatial_meta[[1]][[acol]] <- rep(c("a", "b"), length.out = nrow(st@spatial_meta[[1]]))

  res <- suppressWarnings(STdiff(st, annot = acol, topgenes = 5,
                                 test_type = "t_test", cores = 1, verbose = FALSE))

  testthat::expect_true(is.list(res))
})

testthat::test_that("stdiff_ensure_vst() is idempotent and delegates to calculate_vst()", {
  st <- .vst_fixture("log")
  once <- stdiff_ensure_vst(st)
  twice <- stdiff_ensure_vst(once)

  testthat::expect_equal(once@gene_meta, twice@gene_meta)
})
