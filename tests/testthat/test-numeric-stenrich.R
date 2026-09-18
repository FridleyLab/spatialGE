##
# Numeric assertions for STenrich()'s spatial gene-set enrichment.
#
# See test-numeric-sthet.R's header for why this tier exists.
#
# ORDER SENSITIVITY -- read before pinning anything here. Per fidelity rule 13,
# STenrich p-values move at FULL RESOLUTION under spot reordering:
# computeSubsampleSums() (src/spatialge_helpers.cpp) draws Rcpp::sample() over
# row POSITIONS of coords_mat, so a fixed seed selects different spots when the
# order changes. p-values are therefore reproducible only for a fixed input
# order, which the seeded fixture does provide -- but unlike SThet and
# STgradient there is no order-insensitive variant to fall back on.
#
# Consequently the emphasis here is inverted relative to the other files: the
# load-bearing assertions are the quantization, bounds and a negative control,
# and the pinned p-values are secondary and explicitly documented as
# order-dependent. Do NOT add a permute-and-rerun test expecting invariance --
# it would fail by design, and dev/xver/30_order_probe.R already characterizes
# the effect.
#

stenrich_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
          cores = 1, verbose = FALSE)
}

## STenrich scores whether the spots highly expressing a gene set's MEAN are
## spatially closer than random spots, so what it detects is set COHERENCE, not
## whether the individual genes are spatially structured. These four sets
## separate those ideas, and the measured medians across the 3 samples are:
##
##   gradient     3 genes, all ramping the same way   0.020
##   domainA_only 1 gene, one domain                  0.045
##   domainAB     2 genes, opposite domains           0.545
##   flat         3 genes, no spatial term            0.325
##
## domainAB scoring WORSE than the flat control is correct, not a defect:
## DomainGeneA is high on the left and DomainGeneB high on the right, so their
## per-spot mean is nearly uniform and the signal cancels. Asserting that
## explicitly is more valuable than asserting a p-value, because it pins what
## the statistic means.
STENRICH_SETS <- list(
  gradient     = c("GradGene1", "GradGene2", "GradGene3"),
  domainA_only = c("DomainGeneA"),
  domainAB     = c("DomainGeneA", "DomainGeneB"),
  flat         = c("Gene10", "Gene11", "Gene12")
)

run_stenrich <- function(sc, reps = 100){
  set.seed(20260115)
  STenrich(sc, gene_sets = STENRICH_SETS, reps = reps, min_units = 1,
           min_genes = 1, seed = 123, cores = 1, verbose = FALSE)
}

testthat::test_that("STenrich p-values are properly bounded and quantized by reps", {
  res <- run_stenrich(stenrich_fixture(), reps = 100)
  d <- res[["alpha"]]

  testthat::expect_gt(nrow(d), 0)
  testthat::expect_true(all(d$p_value >= 0 & d$p_value <= 1))
  testthat::expect_true(all(d$adj_p_value >= 0 & d$adj_p_value <= 1))

  ## A permutation p-value over `reps` draws can only take multiples of 1/reps.
  ## This is the cheapest check that the permutation count is actually honored
  ## -- a silently ignored reps= would show up here immediately.
  scaled <- d$p_value * 100
  testthat::expect_true(all(abs(scaled - round(scaled)) < 1e-9))

  ## Adjustment can only raise a p-value.
  testthat::expect_true(all(d$adj_p_value >= d$p_value - 1e-12))
})

stenrich_pooled <- function(reps = 200){
  res <- run_stenrich(stenrich_fixture(), reps = reps)
  do.call(rbind, lapply(names(res), function(s){
    x <- res[[s]]
    if(is.null(x) || !nrow(x)) return(NULL)
    x$samp <- s
    x
  }))
}

testthat::test_that("STenrich detects a coherent spatial gene set above a flat control", {
  ## Pooled across samples: a single sample's p-value at reps=200 is coarse, and
  ## the claim under test is about the statistic, not one draw.
  d <- stenrich_pooled()
  testthat::expect_gt(nrow(d), 6)

  p_med <- tapply(d$p_value, d$gene_set, function(v) stats::median(v, na.rm = TRUE))
  p_max <- tapply(d$p_value, d$gene_set, function(v) max(v, na.rm = TRUE))

  ## GradGene1-3 all ramp with xpos, so their mean ramps too. Significant in
  ## every sample -- asserted on the max, which is the strong form.
  testthat::expect_lt(p_max[["gradient"]], 0.05)

  ## And clearly separated from genes with no spatial term.
  testthat::expect_lt(p_med[["gradient"]], p_med[["flat"]])

  ## The flat control must NOT look significant, or the comparison above proves
  ## nothing about specificity.
  testthat::expect_gt(p_med[["flat"]], 0.10)
})

testthat::test_that("STenrich cancels anti-correlated markers, scoring them below a single marker", {
  d <- stenrich_pooled()
  p_med <- tapply(d$p_value, d$gene_set, function(v) stats::median(v, na.rm = TRUE))

  ## DomainGeneA alone is spatially clustered and scores accordingly. Adding
  ## DomainGeneB -- high in the complementary domain -- flattens the set mean and
  ## destroys the signal. This is the assertion that pins what STenrich measures:
  ## coherence of the set's mean, not structure in its members.
  testthat::expect_lt(p_med[["domainA_only"]], 0.10)
  testthat::expect_gt(p_med[["domainAB"]], p_med[["domainA_only"]])
  testthat::expect_gt(p_med[["domainAB"]], 0.30)
})

testthat::test_that("STenrich handles a gene set with exactly one gene present", {
  ## Regression test. stenrich_score_one_combo() subset with
  ## delayed_mat[genes, ], which drops to a dimensionless vector at length 1, so
  ## DelayedMatrixStats::colMeans2() aborted with "Argument 'dim.' must be an
  ## integer vector of length two". min_genes=1 explicitly permits single-gene
  ## sets -- testing one marker gene for spatial clustering is an ordinary thing
  ## to want -- so this path was reachable and always failed.
  sc <- stenrich_fixture()

  set.seed(20260115)
  testthat::expect_no_error(
    res <- STenrich(sc, gene_sets = list(solo = c("DomainGeneA")), reps = 50,
                    min_units = 1, min_genes = 1, seed = 123, cores = 1, verbose = FALSE)
  )

  d <- res[["alpha"]]
  testthat::expect_equal(nrow(d), 1L)
  testthat::expect_equal(d$size_test, 1L)
  testthat::expect_false(is.na(d$p_value))
  testthat::expect_true(d$p_value >= 0 && d$p_value <= 1)
})

testthat::test_that("STenrich reports set sizes consistent with the requested sets", {
  res <- run_stenrich(stenrich_fixture(), reps = 100)
  d <- res[["alpha"]]

  for(i in seq_len(nrow(d))){
    set_name <- d$gene_set[i]
    ## size_gene_set counts the genes of the set present in the data;
    ## size_test counts those actually used. Neither may exceed the request.
    testthat::expect_lte(d$size_gene_set[i], length(STENRICH_SETS[[set_name]]))
    testthat::expect_lte(d$size_test[i], d$size_gene_set[i])
    testthat::expect_gte(d$size_test[i], 1L)
  }

  ## Column contract (fidelity rule 8).
  for(col in c("sample_name", "gene_set", "size_test", "size_gene_set",
               "prop_size_test", "p_value", "adj_p_value")){
    testthat::expect_true(col %in% colnames(d), label = paste("column", col))
  }
})

testthat::test_that("STenrich matches pinned p-values for the fixture's spot order", {
  res <- run_stenrich(stenrich_fixture(), reps = 100)
  d <- res[["alpha"]]

  ## ORDER-DEPENDENT PINS. Valid only for the seeded fixture's spot order; see
  ## the file header. If these fail after a change that touches spot ordering,
  ## run dev/xver/30_order_probe.R before treating it as a regression.
  pin <- function(set, p, adj){
    r <- d[d$gene_set == set, ]
    testthat::expect_equal(nrow(r), 1L, label = paste("one row for", set))
    testthat::expect_equal(r$p_value, p, tolerance = 1e-9, label = paste("p_value", set))
    testthat::expect_equal(r$adj_p_value, adj, tolerance = 1e-9, label = paste("adj_p_value", set))
  }

  pin("domainA_only", 0.00, 0.00)
  pin("gradient",     0.04, 0.08)
  pin("domainAB",     0.06, 0.08)
  pin("flat",         0.37, 0.37)
})
