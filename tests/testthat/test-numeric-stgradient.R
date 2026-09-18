##
# Numeric assertions for STgradient()'s distance-vs-expression correlations.
#
# See test-numeric-sthet.R's header for why this tier exists.
#
# ORDER SENSITIVITY. Uses distsumm='min' throughout. Per fidelity rule 13,
# distsumm='avg' averages distances that then tie exactly, so a ~1e-16
# summation-order change flips a tie and therefore a rank -- moving Spearman's r
# by up to 4.3e-3 (dev/xver/REPORT.md). Pinning an 'avg' value would be pinning
# noise. distsumm='min' agrees to ~1e-15 under reordering.
#
# ref is given explicitly as "1" rather than derived from unique(labels)[1]:
# unique() returns labels in order of appearance, so deriving the reference that
# way makes the whole result depend on spot order even with distsumm='min'.
#

stgradient_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
          cores = 1, verbose = FALSE)
}

STG_ANNOT <- "stclust_spw0.025_dsplFalse"

testthat::test_that("STgradient distsumm='min' matches pinned values", {
  sg <- STgradient(stgradient_fixture(), annot = STG_ANNOT, ref = "1",
                   topgenes = 25, distsumm = "min", cores = 1, verbose = FALSE)
  d <- sg[["alpha"]]
  d <- d[order(d$gene), ]

  pick <- function(gene) d[d$gene == gene, ]

  ## Gene$special carries a shell metacharacter deliberately, so pinning it
  ## also guards the check.names=FALSE path through this function.
  r <- pick("Gene$special")
  testthat::expect_equal(r$min_spearman_r, -0.111076773725, tolerance = 1e-6)
  testthat::expect_equal(r$min_spearman_r_pval, 0.43307224540, tolerance = 1e-6)

  r <- pick("Gene103")
  testthat::expect_equal(r$min_spearman_r, -0.167277810082, tolerance = 1e-6)
  testthat::expect_equal(r$min_spearman_r_pval, 0.23589558403, tolerance = 1e-6)

  r <- pick("Gene107")
  testthat::expect_equal(r$min_spearman_r, 0.163784617042, tolerance = 1e-6)
  testthat::expect_equal(r$min_spearman_r_pval, 0.24596248203, tolerance = 1e-6)
})

testthat::test_that("STgradient recovers the fixture's designed expression gradient", {
  ## topgenes=200 so every GradGene is included -- the variance-based top-25
  ## admits only 2 of them, which is too few to test a trend.
  sg <- STgradient(stgradient_fixture(), annot = STG_ANNOT, ref = "1",
                   topgenes = 200, distsumm = "min", cores = 1, verbose = FALSE)
  d <- sg[["alpha"]]

  grad <- d[grepl("^GradGene", d$gene), ]
  rand <- d[grepl("^Gene[0-9]+$", d$gene), ]
  testthat::expect_equal(nrow(grad), 10L)
  testthat::expect_gt(nrow(rand), 100)

  ## The fixture gives GradGene1..10 a lambda rising linearly in xpos, and the
  ## clustering aligns with that axis, so distance from the reference domain
  ## runs against expression for all ten. Directional consistency is the real
  ## signal: under a null, ten independent signs agreeing has probability
  ## ~0.001. Allowing one exception keeps it robust without weakening it much.
  n_negative <- sum(grad$min_spearman_r < 0, na.rm = TRUE)
  testthat::expect_gte(n_negative, 9L)

  ## And the effect is stronger than for genes with no spatial term. Measured
  ## medians are 0.203 vs 0.089; compared as medians rather than means so a
  ## single outlier gene cannot carry it.
  testthat::expect_gt(stats::median(abs(grad$min_spearman_r), na.rm = TRUE),
                      stats::median(abs(rand$min_spearman_r), na.rm = TRUE))
})

testthat::test_that("STgradient adjusted p-values are valid and monotone", {
  sg <- STgradient(stgradient_fixture(), annot = STG_ANNOT, ref = "1",
                   topgenes = 50, distsumm = "min", cores = 1, verbose = FALSE)
  d <- sg[["alpha"]]

  ok <- !is.na(d$min_spearman_r_pval) & !is.na(d$min_spearman_r_pval_adj)
  testthat::expect_gt(sum(ok), 20)

  testthat::expect_true(all(d$min_spearman_r_pval[ok] >= 0))
  testthat::expect_true(all(d$min_spearman_r_pval[ok] <= 1))
  testthat::expect_true(all(d$min_spearman_r_pval_adj[ok] >= d$min_spearman_r_pval[ok] - 1e-12))

  o <- order(d$min_spearman_r_pval[ok])
  testthat::expect_false(any(diff(d$min_spearman_r_pval_adj[ok][o]) < -1e-12))

  ## Spearman's r is a correlation and must stay in [-1, 1] -- worth pinning
  ## because the planned C++ rewrite computes the ranks itself.
  r_ok <- !is.na(d$min_spearman_r)
  testthat::expect_true(all(abs(d$min_spearman_r[r_ok]) <= 1 + 1e-12))
})

testthat::test_that("STgradient distsumm='min' is invariant to spot order", {
  sc <- stgradient_fixture()
  base <- STgradient(sc, annot = STG_ANNOT, ref = "1", topgenes = 25,
                     distsumm = "min", cores = 1, verbose = FALSE)

  sc2 <- sc
  perm <- rev(seq_len(ncol(sc2@tr_counts[["alpha"]])))
  sc2@counts[["alpha"]] <- sc2@counts[["alpha"]][, perm, drop = FALSE]
  sc2@tr_counts[["alpha"]] <- sc2@tr_counts[["alpha"]][, perm, drop = FALSE]
  sm <- sc2@spatial_meta[["alpha"]][perm, , drop = FALSE]
  ## Reset row names: base `[` on a plain data.frame leaves them
  ## non-contiguous, which downstream column_to_rownames() calls reject.
  rownames(sm) <- NULL
  sc2@spatial_meta[["alpha"]] <- sm

  permuted <- STgradient(sc2, annot = STG_ANNOT, ref = "1", topgenes = 25,
                         distsumm = "min", cores = 1, verbose = FALSE)

  b <- base[["alpha"]]; p <- permuted[["alpha"]]
  common <- intersect(b$gene, p$gene)
  testthat::expect_gt(length(common), 20)

  ## This is what licenses the pinned values above: 'min' takes a minimum over
  ## distances, which no summation-order change can perturb, unlike 'avg'.
  for(g in common){
    testthat::expect_equal(p$min_spearman_r[p$gene == g],
                           b$min_spearman_r[b$gene == g],
                           tolerance = 1e-10, label = paste("order-invariant", g))
  }
})
