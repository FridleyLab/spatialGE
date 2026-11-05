##
# Tests for sparse_row_sds(), which replaces apply(m, 1, sd) -- a silent full
# densification, since base apply() calls as.matrix() on any non-array.
#
# The equivalence asserted here is "agrees with apply(m, 1, sd) to floating-point
# tolerance", NOT bit-identity. The summation order genuinely differs. Where the
# two disagree, sparse_row_sds() is the more accurate of the two: it adds the
# structural zeros' contribution as a single (n - nnz) * mu^2 term rather than
# summing n - nnz copies of it, so the dense path accumulates more rounding over
# a long row. Measured max relative disagreement is 9.6e-16 at 120 columns and
# 5.3e-13 at 50,000.
#

sp <- function(M){
  if(is.null(rownames(M))) rownames(M) <- paste0("g", seq_len(nrow(M)))
  methods::as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix")
}

testthat::test_that("sparse_row_sds agrees with apply(m, 1, sd) on a realistic matrix", {
  set.seed(20260115)
  M <- matrix(as.double(rpois(200 * 60, 3)), nrow = 200)
  M[M < 2] <- 0                      # force genuine sparsity
  m <- sp(M)

  testthat::expect_equal(unname(sparse_row_sds(m)), unname(apply(m, 1, sd)),
                         tolerance = 1e-12)
  testthat::expect_named(sparse_row_sds(m), rownames(m))
})

testthat::test_that("sparse_row_sds matches apply() on the degenerate rows that matter", {
  ## Each of these is a row shape the fixture or real data actually produces, and
  ## each is where a one-pass sum-of-squares shortcut would misbehave -- an
  ## all-zero or constant row can drive sqrt() to a negative radicand.
  cases <- list(
    all_zero      = rbind(c(0, 0, 0, 0), c(1, 2, 3, 4)),
    constant      = rbind(c(7, 7, 7, 7), c(0, 0, 1, 0)),
    one_nonzero   = rbind(c(0, 0, 5, 0), c(0, 0, 0, 0)),
    negatives     = rbind(c(-1, 0, 1, 0), c(-5, -5, -5, -5)),
    large_values  = rbind(c(1e8, 1e8 + 1, 1e8, 1e8 + 2), c(0, 0, 0, 1))
  )

  for(nm in names(cases)){
    m <- sp(cases[[nm]])
    testthat::expect_equal(unname(sparse_row_sds(m)), unname(apply(m, 1, sd)),
                           tolerance = 1e-8, label = paste("row sds:", nm))
  }
})

testthat::test_that("sparse_row_sds returns NA for a single column, like stats::sd", {
  m <- sp(matrix(c(1, 2, 3), ncol = 1))

  ## sd of one observation is NA. A naive (n-1) denominator would divide by zero
  ## and yield NaN or Inf instead, which would then propagate into @gene_meta.
  res <- sparse_row_sds(m)
  testthat::expect_true(all(is.na(res)))
  testthat::expect_equal(unname(res), unname(apply(m, 1, sd)))
})

testthat::test_that("sparse_row_sds allocates far less than the apply() it replaces", {
  ## The point of the function: apply(m, 1, sd) coerces to a dense array first.
  ##
  ## Measured as a RATIO against apply() on the same matrix, not against an
  ## absolute byte bound. gc()'s "max used" is a session-wide high-water mark, so
  ## an absolute reading picks up whatever the session already did and is
  ## meaningless here; the difference between two adjacent measurements is not.
  ## Timing is not used either -- too noisy on a shared CI runner.
  set.seed(20260115)
  m <- methods::as(Matrix::rsparsematrix(2000, 2000, density = 0.01), "CsparseMatrix")

  measure <- function(expr){
    gc(reset = TRUE, full = TRUE)
    before <- sum(gc(full = TRUE)[, 6])
    force(expr)
    after <- sum(gc(full = TRUE)[, 6])
    after - before
  }

  dense_alloc <- measure(apply(m, 1, sd))
  sparse_alloc <- measure(sparse_row_sds(m))

  ## Asserted purely as a ratio, with no absolute bound: gc()'s column 6 is not
  ## reliably in MB across R builds, so an absolute threshold is a unit bug
  ## waiting to happen. A ratio compares like with like and needs no units.
  ## Measured ~150x (14,405,808 vs 96,114 in gc's own units) at 2000 x 2000 with
  ## 1% density; requiring only 10x leaves ample headroom while still failing
  ## outright if the densification returns.
  testthat::expect_gt(dense_alloc, 0)
  testthat::expect_lt(sparse_alloc, dense_alloc / 10)
})

testthat::test_that("sparse_row_sds accepts a non-CSC input by coercing it", {
  M <- rbind(c(0, 1, 2), c(3, 0, 0))
  rownames(M) <- c("a", "b")

  base_mat <- sparse_row_sds(M)
  triplet <- sparse_row_sds(methods::as(Matrix::Matrix(M, sparse = TRUE), "TsparseMatrix"))
  csc <- sparse_row_sds(sp(M))

  testthat::expect_equal(base_mat, csc, tolerance = 1e-12)
  testthat::expect_equal(triplet, csc, tolerance = 1e-12)
})

testthat::test_that("gene_stdevs in @gene_meta is produced by sparse_row_sds", {
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)

  ## Pins the wiring, so a revert to apply(m, 1, sd) is caught as a change here
  ## rather than only as an unexplained golden diff.
  for(s in names(st@gene_meta)){
    gm <- st@gene_meta[[s]]
    expected <- sparse_row_sds(sample_tr_counts(st, s))
    testthat::expect_equal(gm$gene_stdevs, unname(expected[gm$gene]),
                           tolerance = 0, label = paste("gene_stdevs wiring", s))
  }
})

testthat::test_that("log_transf_one_sample returns a well-formed dgCMatrix", {
  set.seed(20260115)
  M <- matrix(as.double(rpois(200, 3)), nrow = 20,
              dimnames = list(paste0("g", 1:20), paste0("c", 1:10)))
  M[M < 2] <- 0
  m <- methods::as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix")

  res <- log_transf_one_sample(m, "s1", 10000)$counts

  testthat::expect_s4_class(res, "dgCMatrix")

  ## @x must be an UNNAMED numeric vector. Matrix::colSums() returns a vector
  ## named by colnames, so scaling by libsizes[col] silently propagates those
  ## names into @x -- values stay correct, so an as.matrix() comparison reports
  ## the two as identical and only a slot-level check catches it.
  testthat::expect_null(names(res@x))
  testthat::expect_type(res@x, "double")

  ## Sparsity preserved, not destroyed and rebuilt: log1p(0) == 0, so the stored
  ## pattern must be unchanged.
  testthat::expect_equal(length(res@x), length(m@x))
  testthat::expect_equal(res@i, m@i)
  testthat::expect_equal(res@p, m@p)
  testthat::expect_identical(dimnames(res), dimnames(m))

  ## And bit-identical to the dense sweep/multiply/log1p sequence it replaced.
  d <- as.matrix(m)
  l <- colSums(d)
  expected <- log1p(sweep(d, 2, l, "/") * 10000)
  testthat::expect_identical(as.matrix(res), expected)
})

testthat::test_that("expandSparse(transpose=TRUE) is identical to transpose-then-frame", {
  ## Equivalence guard for the item-5 rewrite: the transposing callers used to do
  ## expandSparse(m) %>% t() %>% as.data.frame(), which densified, framed,
  ## transposed the frame, and framed again. transpose=TRUE densifies the sparse
  ## transpose once and must return the same data frame, names included.
  set.seed(20260115)
  M <- matrix(as.double(rpois(5 * 4, 3)), nrow = 5,
              dimnames = list(c("Gene1", "Gene$special", "Gene3", "Gene4", "Gene5"),
                              c("cellA", "cellB", "cell C", "cellD")))
  m <- methods::as(Matrix::Matrix(M, sparse = TRUE), "CsparseMatrix")

  old_way <- as.data.frame(t(data.frame(as.matrix(m), check.names = FALSE)))
  new_way <- expandSparse(m, transpose = TRUE)
  testthat::expect_identical(new_way, old_way)

  ## The default path is unchanged.
  testthat::expect_identical(expandSparse(m), data.frame(as.matrix(m), check.names = FALSE))
})
