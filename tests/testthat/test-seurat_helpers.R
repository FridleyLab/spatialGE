##
# Regression tests for Seurat_FindVariableFeatures()'s input coercion.
#
# The VST statistics this produces land in @gene_meta as vst.* columns and drive
# topgenes selection in STclust, STgradient, gene_interpolation and
# STdiff_prepare, so the coercion path is load-bearing even though every current
# caller happens to hand it a dgCMatrix already.
#

sfvf_dense <- function(){
  set.seed(20260115)
  m <- matrix(as.double(rpois(200, 4)), nrow = 20,
              dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
  m[1, ] <- 0      # all-zero gene -> variance 0, exercises the not.const branch
  m[2, ] <- 5      # constant gene -> variance 0
  m
}

testthat::test_that("Seurat_FindVariableFeatures accepts a base matrix", {
  ## Before the fix this died with 'could not find function "as.sparse"'. The
  ## two coercion branches could not work together: the first turned a non-Matrix
  ## input into a dgeMatrix, the second then saw a non-dgCMatrix and called
  ## as.sparse(), which does not exist in this package or any of its imports. It
  ## stayed hidden because 'as.sparse' was listed in utils::globalVariables(),
  ## suppressing the R CMD check NOTE that would have caught it.
  m <- sfvf_dense()

  testthat::expect_no_error(hvf <- Seurat_FindVariableFeatures(m))
  testthat::expect_equal(nrow(hvf), 20L)
  testthat::expect_true(all(c("vst.mean", "vst.variance", "vst.variance.expected",
                              "vst.variance.standardized") %in% colnames(hvf)))
  testthat::expect_true(all(is.finite(hvf$vst.mean)))
})

testthat::test_that("Seurat_FindVariableFeatures accepts a dense Matrix and a data.frame", {
  m <- sfvf_dense()

  ## dgeMatrix: a dense Matrix. This is exactly what the old first branch
  ## produced and the old second branch then choked on.
  dge <- methods::as(m, "denseMatrix")
  testthat::expect_false(inherits(dge, "dgCMatrix"))
  testthat::expect_no_error(from_dge <- Seurat_FindVariableFeatures(dge))

  testthat::expect_no_error(from_df <- Seurat_FindVariableFeatures(as.data.frame(m)))

  ## All three routes must agree exactly with the sparse one -- coercion may not
  ## change the statistic.
  from_sparse <- Seurat_FindVariableFeatures(methods::as(m, "CsparseMatrix"))
  testthat::expect_equal(from_dge, from_sparse, tolerance = 1e-12)
  testthat::expect_equal(from_df, from_sparse, tolerance = 1e-12)
})

testthat::test_that("Seurat_FindVariableFeatures coerces to a general double CSC matrix", {
  ## The C++ takes Eigen::SparseMatrix<double>, so the coercion must land on
  ## dgCMatrix specifically -- not lgCMatrix and not dsCMatrix. A plain
  ## as(x, 'CsparseMatrix') preserves logical and symmetric structure, which is
  ## why SeuratObject::as.sparse() routes through dMatrix then generalMatrix, and
  ## why this reproduces both steps.
  ## Matrices are 20 rows, not a handful: stats::loess with span=0.3 raises
  ## "span is too small" on a couple of points, which would mask the coercion
  ## behavior under test.
  M <- sfvf_dense()

  logical_sparse <- methods::as(Matrix::Matrix(M > 0, sparse = TRUE), "CsparseMatrix")
  testthat::expect_true(inherits(logical_sparse, "lgCMatrix"))
  testthat::expect_no_error(Seurat_FindVariableFeatures(logical_sparse))

  set.seed(20260115)
  A <- matrix(as.double(rpois(400, 4)), nrow = 20)
  S <- A + t(A)                                    # force exact symmetry
  dimnames(S) <- list(paste0("g", 1:20), paste0("g", 1:20))
  sym <- methods::as(methods::as(S, "dMatrix"), "symmetricMatrix")
  testthat::expect_false(inherits(sym, "generalMatrix"))
  testthat::expect_no_error(Seurat_FindVariableFeatures(sym))

  ## Integer input must become double rather than failing inside RcppEigen.
  ## Drawn at random rather than patterned: rep(1:10, 20) fills column-major and
  ## makes every ROW constant, which leaves loess with zero non-constant rows and
  ## fails with "invalid 'x'" for reasons unrelated to the coercion.
  set.seed(20260116)
  int_mat <- matrix(as.integer(rpois(200, 4)), nrow = 20,
                    dimnames = list(paste0("g", 1:20), NULL))
  testthat::expect_type(int_mat[1, 1], "integer")
  testthat::expect_no_error(Seurat_FindVariableFeatures(int_mat))
})

testthat::test_that("Seurat_FindVariableFeatures leaves an existing dgCMatrix untouched", {
  m <- methods::as(sfvf_dense(), "CsparseMatrix")

  ## The no-op path every current caller takes. Asserted so a future change to
  ## the coercion cannot silently start densifying the common case -- which
  ## would be a memory regression rather than a wrong answer, and so would not
  ## show up in any value comparison.
  hvf <- Seurat_FindVariableFeatures(m)
  testthat::expect_equal(nrow(hvf), nrow(m))
  testthat::expect_equal(hvf$vst.mean, as.vector(Matrix::rowMeans(m)), tolerance = 1e-12)
})
