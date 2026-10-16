##
# VST (standardized variance) availability and stability.
#
# `gene_meta` used to hold only gene_mean/gene_stdevs after transform_data();
# `vst.variance.standardized` appeared later, as a side effect of whichever
# analysis ran first (STclust(), gene_interpolation()). transform_data() now
# populates VST itself, so get_gene_meta() shows it regardless of call order.
#
# Two corrections to an earlier version of this comment, both measured: the
# gene_stdevs-vs-VST top-20 overlap is 7-13 across the eight tnbc_bassiouni
# samples (8-11 on this fixture), not 6/20; and the set of genes a function
# *selects* never depended on call order, because STclust(), STgradient(),
# gene_interpolation() and STdiff() each computed VST themselves before ranking.
# What the change fixed is where the column comes from, not what is done with it.
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

## An STlist with `nvar` varying genes plus `nconst` constant ones, on a 10x10
## grid. Deliberately tiny: the variance-mean loess cannot be fit below 4
## non-constant genes.
.vst_tiny_fixture <- function(nvar, nconst = 0) {
  set.seed(20260115)
  xy <- expand.grid(x = 1:10, y = 1:10)
  libs <- paste0("spot", seq_len(nrow(xy)))
  rows <- list()
  for (i in seq_len(nvar)) rows[[length(rows) + 1L]] <- stats::rpois(nrow(xy), 5 + i)
  for (i in seq_len(nconst)) rows[[length(rows) + 1L]] <- rep(3, nrow(xy))
  mtx <- do.call(rbind, rows)
  expr <- data.frame(genename = paste0("gene_", seq_len(nrow(mtx))), mtx)
  colnames(expr) <- c("genename", libs)
  coords <- data.frame(libname = libs, y = xy$y, x = xy$x)
  STlist(rnacounts = list(tiny = expr), spotcoords = list(tiny = coords),
         cores = 1, verbose = FALSE)
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

testthat::test_that("transform_data() warns and skips VST when too few genes vary", {
  ## Populating VST here means transform_data() now fits a loess ACROSS GENES,
  ## which has a minimum data requirement the rest of the function does not: at
  ## span=0.3 it needs 4 non-constant genes, and below that stats::loess stops
  ## with 'span is too small'. That made transform_data() fail on an STlist it
  ## handled in 1.2.2. The simulated 3-gene example in
  ## basic_functions_vignette.Rmd is exactly this case, which is how it was
  ## found -- the vignette could no longer be built.
  st <- .vst_tiny_fixture(nvar = 3)

  w <- testthat::capture_warnings(st <- transform_data(st, cores = 1))
  testthat::expect_true(any(grepl("non-zero variance", w)))
  testthat::expect_true(any(grepl("tiny", w)))   # names the sample

  ## The transformation itself, which is what the caller asked for, succeeded.
  testthat::expect_false(rlang::is_empty(st@tr_counts))
  testthat::expect_identical(dim(st@tr_counts[["tiny"]]), dim(st@counts[["tiny"]]))
  testthat::expect_true(all(c("gene_mean", "gene_stdevs") %in% colnames(st@gene_meta[["tiny"]])))

  ## ...and the column that could not be fit is absent rather than zero-filled.
  ## calculate_vst() only computes when it is missing, so a later analysis
  ## retries instead of inheriting a meaningless value.
  testthat::expect_false("vst.variance.standardized" %in% colnames(st@gene_meta[["tiny"]]))

  ## SThet() does not use VST, so the vignette's actual goal still works.
  st <- SThet(st, genes = paste0("gene_", 1:3), method = "moran", cores = 1, verbose = FALSE)
  testthat::expect_true("moran_i" %in% colnames(get_gene_meta(st, sthet_only = TRUE)))
})

testthat::test_that("it is the non-constant gene count that gates the VST fit", {
  ## A guard on nrow() would have missed this: the loess is fit on the
  ## non-constant genes only, so a matrix with plenty of genes of which just
  ## three vary fails in exactly the same way.
  st <- .vst_tiny_fixture(nvar = 3, nconst = 97)
  testthat::expect_identical(nrow(st@counts[["tiny"]]), 100L)

  w <- testthat::capture_warnings(st <- transform_data(st, cores = 1))
  testthat::expect_true(any(grepl("non-zero variance", w)))
  testthat::expect_false("vst.variance.standardized" %in% colnames(st@gene_meta[["tiny"]]))

  ## One more varying gene is enough to fit it.
  ok <- suppressWarnings(transform_data(.vst_tiny_fixture(nvar = 4), cores = 1))
  testthat::expect_true("vst.variance.standardized" %in% colnames(ok@gene_meta[["tiny"]]))
})

testthat::test_that("calculate_vst() on a matrix errors when the trend cannot be fit", {
  ## The matrix form exists to rank genes for a caller's topgenes, so there is no
  ## useful degraded answer: error when nothing matches, warn when only some does.
  st <- .vst_tiny_fixture(nvar = 3)

  testthat::expect_error(calculate_vst(x = st@counts[["tiny"]]),
                         "spatialGE exception.*too few to rank genes by variability")

  ## And it still works on a matrix that can be fit, returning the gene/statistic
  ## frame STgradient ranks (a plain data.frame, not a tibble).
  fit <- calculate_vst(x = .vst_tiny_fixture(nvar = 6)@counts[["tiny"]])
  testthat::expect_identical(colnames(fit), c("gene", "vst.variance.standardized"))
  testthat::expect_identical(nrow(fit), 6L)
})
