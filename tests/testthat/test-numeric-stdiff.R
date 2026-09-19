##
# Numeric assertions for STdiff()'s non-spatial differential expression.
#
# See test-numeric-sthet.R's header for why this tier exists. Only test_type
# 't_test' is covered: the 'mm' path is spaMM::fitme with a Matern correlation
# structure, measured at 20.8 s and 27.4 GB for a single fit
# (dev/profiling/FINDINGS.md), which is far outside a default test budget.
#
# The strongest assertions here are invariants rather than pinned values. A
# two-cluster comparison must be antisymmetric in log fold-change and identical
# in p-value between the two directions; that has to hold under any
# implementation, including the vectorized rewrite planned for
# non_spatial_de(), whereas pinned values only detect that something moved.
#

stdiff_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
          cores = 1, verbose = FALSE)
}

run_stdiff <- function(sc){
  STdiff(sc, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 20,
         test_type = "t_test", cores = 1)
}

testthat::test_that("STdiff t_test is antisymmetric between the two comparison directions", {
  res <- run_stdiff(stdiff_fixture())
  d <- res[[1]]

  ## alpha has exactly 2 clusters, so each gene yields one row per direction.
  a <- d[d$sample == "alpha", ]
  testthat::expect_gt(nrow(a), 0)

  split_by_gene <- split(a, a$gene)
  two_way <- Filter(function(g) nrow(g) == 2, split_by_gene)
  testthat::expect_gt(length(two_way), 10)

  for(g in two_way){
    ## log2 fold-change must be exactly negated between directions...
    testthat::expect_equal(g$avg_log2fc[1], -g$avg_log2fc[2], tolerance = 1e-10,
                           label = paste("antisymmetric log2fc", g$gene[1]))
    ## ...and a two-sided t-test's p-value must be direction-independent.
    testthat::expect_equal(g$ttest_p_val[1], g$ttest_p_val[2], tolerance = 1e-12,
                           label = paste("symmetric p", g$gene[1]))
  }
})

testthat::test_that("STdiff t_test matches pinned values", {
  res <- run_stdiff(stdiff_fixture())
  d <- res[[1]]
  d <- d[order(d$sample, d$gene, d$cluster_1), ]

  pick <- function(gene, cl){
    r <- d[d$sample == "alpha" & d$gene == gene & d$cluster_1 == cl, ]
    testthat::expect_equal(nrow(r), 1L, label = paste("one row for", gene, cl))
    r
  }

  ## Gene$special carries a shell metacharacter on purpose (it proves
  ## check.names=FALSE survives end to end), so pinning it also guards the
  ## name-handling path.
  r <- pick("Gene$special", "1")
  testthat::expect_equal(r$avg_log2fc, -0.198678027, tolerance = 1e-6)
  testthat::expect_equal(r$ttest_p_val, 0.24862712, tolerance = 1e-6)
  testthat::expect_equal(r$adj_p_val, 0.38250327, tolerance = 1e-6)

  r <- pick("Gene103", "1")
  testthat::expect_equal(r$avg_log2fc, -0.085882348, tolerance = 1e-6)
  testthat::expect_equal(r$ttest_p_val, 0.62287448, tolerance = 1e-6)
  testthat::expect_equal(r$adj_p_val, 0.73279351, tolerance = 1e-6)
})

testthat::test_that("STdiff adjusted p-values are monotone in the raw p-values", {
  res <- run_stdiff(stdiff_fixture())
  d <- res[[1]]

  ## Benjamini-Hochberg is monotone within the set it adjusts over, and
  ## adjustment can only increase a p-value. Both are cheap to check and catch
  ## an adjustment applied over the wrong grouping -- which pinned values would
  ## not reveal, since they would simply all shift together.
  testthat::expect_true(all(d$adj_p_val >= d$ttest_p_val - 1e-12))
  testthat::expect_true(all(d$adj_p_val <= 1 + 1e-12))
  testthat::expect_true(all(d$ttest_p_val >= 0))

  for(s in unique(d$sample)){
    ds <- d[d$sample == s, ]
    o <- order(ds$ttest_p_val)
    adj <- ds$adj_p_val[o]
    testthat::expect_false(any(diff(adj) < -1e-12), label = paste("monotone adj_p", s))
  }
})

testthat::test_that("STdiff returns the documented column contract", {
  res <- run_stdiff(stdiff_fixture())
  d <- res[[1]]

  ## Column position and type are API (fidelity rule 8).
  for(col in c("sample", "gene", "avg_log2fc", "cluster_1", "ttest_p_val", "adj_p_val")){
    testthat::expect_true(col %in% colnames(d), label = paste("column", col))
  }
  testthat::expect_type(d$avg_log2fc, "double")
  testthat::expect_type(d$ttest_p_val, "double")
  testthat::expect_type(d$adj_p_val, "double")
})

testthat::test_that("STdiff wilcoxon matches the t_test log-fold-change contract", {
  ## The wilcoxon path shares stdiff_mean_test() with t_test. It was uncovered
  ## before the item-2 rewrite, which replaced the per-combination dplyr passes
  ## with base indexing; the two paths must agree on avg_log2fc exactly (same
  ## group means) and both must produce a bounded, direction-independent p-value.
  ## STclust()'s internal joins emit dplyr many-to-many warnings on this
  ## fixture, and wilcox.test warns "cannot compute exact p-value with ties";
  ## neither is what this test is about.
  sc <- suppressWarnings(stdiff_fixture())
  res_t <- suppressWarnings(STdiff(sc, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 20,
                                   test_type = "t_test", cores = 1))
  res_w <- suppressWarnings(STdiff(sc, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 20,
                                   test_type = "wilcoxon", cores = 1))

  dt <- res_t[[1]]; dw <- res_w[[1]]
  key_t <- paste(dt$sample, dt$gene, dt$cluster_1)
  key_w <- paste(dw$sample, dw$gene, dw$cluster_1)
  testthat::expect_setequal(key_t, key_w)

  m <- match(key_t, key_w)
  testthat::expect_equal(dt$avg_log2fc, dw$avg_log2fc[m], tolerance = 1e-12,
                         label = "same log-fold-change")
  testthat::expect_true(all(dw$wilcoxon_p_val >= 0 & dw$wilcoxon_p_val <= 1))
  testthat::expect_true(all(dw$adj_p_val >= dw$wilcoxon_p_val - 1e-12))

  ## Two-sided wilcoxon is direction-independent, like the t-test.
  a <- dw[dw$sample == "alpha", ]
  two_way <- Filter(function(g) nrow(g) == 2, split(a, a$gene))
  testthat::expect_gt(length(two_way), 5)
  for(g in two_way){
    testthat::expect_equal(g$wilcoxon_p_val[1], g$wilcoxon_p_val[2], tolerance = 1e-12,
                           label = paste("symmetric wilcoxon p", g$gene[1]))
  }
})
