##
# Regression tests for pseudobulk_samples()'s zero-variance-gene handling, and
# numeric assertions on the PCA it produces.
#
# These run on the committed 6-sample pseudobulk fixture (>= 4 samples, which
# pseudobulk_samples() requires), so unlike the golden gate they are
# self-contained, platform-portable, and run by default on CI. Numeric
# assertions carry explicit tolerances for that reason -- the PCA goes through
# prcomp/LAPACK, which is not bit-reproducible across BLAS implementations.
#

pb_stlist <- function(cores = 1){
  base <- fixture_pbulk_base_data()
  fx <- fixture_write_pbulk(base, file.path(tempdir(), "test-pbulk"))
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               samples = fx$samples, cores = cores, verbose = FALSE)
  transform_data(st, method = "log", cores = cores)
}

testthat::test_that("pseudobulk_samples survives genes with zero variance across samples", {
  st <- pb_stlist()

  ## ZeroGene (all-zero) and ConstGene (constant 7) both have an identical
  ## pseudobulk total in every sample. Before the fix, vargenes retained them
  ## whenever max_var_genes >= the gene count, scale() produced an all-NaN
  ## column, and prcomp() aborted with "cannot rescale a constant/zero column
  ## to unit variance" -- making the function unusable on legitimate input.
  testthat::expect_true(all(c("ZeroGene", "ConstGene") %in% rownames(st@counts[[1]])))

  testthat::expect_warning(
    pb <- pseudobulk_samples(st, calc_umap = FALSE),
    "zero variance across samples"
  )

  testthat::expect_s4_class(pb, "STlist")
  testthat::expect_false(is.null(pb@misc[["pbulk_pca"]]))
  testthat::expect_true(all(is.finite(as.matrix(pb@misc[["pbulk_pca"]]))))
})

testthat::test_that("pseudobulk_samples drops the zero-variance genes rather than keeping NaN columns", {
  st <- pb_stlist()
  pb <- suppressWarnings(pseudobulk_samples(st, calc_umap = FALSE))

  mtx <- pb@misc[["scaled_pbulk_mtx"]]
  testthat::expect_false(any(is.nan(mtx)))

  ## ZeroGene is all-zero, so it is zero-variance before AND after the
  ## library-size normalization pseudobulk applies -- it is the only gene that
  ## actually reaches scale() with sd == 0.
  testthat::expect_false("ZeroGene" %in% colnames(mtx))

  ## ConstGene is deliberately NOT dropped, and that distinction matters: it is
  ## constant in raw counts, but pseudobulk divides by library size, and
  ## groupB's 3x lift on the signal genes changes its library size. So
  ## ConstGene's normalized pseudobulk genuinely differs between groups and is
  ## legitimately variable. A fix that dropped it would be over-eager.
  testthat::expect_true("ConstGene" %in% colnames(mtx))

  testthat::expect_equal(nrow(mtx), 6L)
  testthat::expect_equal(ncol(mtx), length(fixture_gene_names()) - 1L)
})

testthat::test_that("pseudobulk PCA separates the two fixture groups on PC1", {
  st <- pb_stlist()
  pb <- suppressWarnings(pseudobulk_samples(st, calc_umap = FALSE))

  pca <- pb@misc[["pbulk_pca"]]
  testthat::expect_equal(nrow(pca), 6L)

  pc1 <- stats::setNames(pca$PC1, rownames(pca))
  a <- pc1[c("pbA1", "pbA2", "pbA3")]
  b <- pc1[c("pbB1", "pbB2", "pbB3")]

  ## The fixture puts a 3x lift on 42 signal genes in groupB, so PC1 is the
  ## group axis. Assert separation by sign agreement within group and a gap
  ## between groups, not on absolute coordinates -- prcomp's eigenvector sign
  ## is arbitrary and can flip between LAPACK builds.
  testthat::expect_true(all(sign(a) == sign(a[1])))
  testthat::expect_true(all(sign(b) == sign(b[1])))
  testthat::expect_true(sign(a[1]) != sign(b[1]))
  testthat::expect_gt(min(abs(range(a) - range(b))), 20)

  ## PC1 carries essentially all the variance given how strong the group
  ## effect is. Loose bound so a BLAS difference cannot trip it.
  var_expl <- pb@misc[["pbulk_pca_var"]]
  testthat::expect_gt(var_expl[["PC1"]], 0.90)
  testthat::expect_lt(var_expl[["PC2"]], 0.05)
})

testthat::test_that("pseudobulk_samples still refuses fewer than four samples", {
  st <- pb_stlist()
  st3 <- filter_data(st, spot_minreads = 0, gene_minreads = 0,
                     samples = c("pbA1", "pbA2", "pbA3"))
  st3@counts <- st3@counts[c("pbA1", "pbA2", "pbA3")]
  st3@tr_counts <- st3@tr_counts[c("pbA1", "pbA2", "pbA3")]

  testthat::expect_error(pseudobulk_samples(st3, calc_umap = FALSE), "four samples")
})

testthat::test_that("pseudobulk_heatmap returns a drawable HeatmapList with the requested gene count", {
  st <- pb_stlist()
  pb <- suppressWarnings(pseudobulk_samples(st, calc_umap = FALSE))

  hm <- pseudobulk_heatmap(pb, hm_display_genes = 5)
  testthat::expect_s4_class(hm, "HeatmapList")

  ## The displayed matrix is at @ht_list[[1]]@matrix, not @matrix --
  ## pseudobulk_heatmap() ends in ComplexHeatmap::draw().
  mtx <- hm@ht_list[[1]]@matrix
  testthat::expect_equal(dim(mtx), c(5L, 6L))
  testthat::expect_setequal(colnames(mtx), fixture_pbulk_sample_names())
})
