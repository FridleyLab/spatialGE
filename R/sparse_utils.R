# Sparse-matrix conversion helpers shared by STlist() and its readers.

##
# makeSparse: takes a dataframe input and makes it sparse saving space
# wrapped from Matrix package
# @param dataframe, a dataframe or tibble with gene names in the first column and
# different spots in the following columns
# @return a sparsed data matrix
#
makeSparse = function(dataframe){

  # To prevent NOTES in R CMD check
  . = NULL

  if(any(class(dataframe) == 'dgCMatrix')){
    numdat = dataframe
  } else if(!is.matrix(dataframe)){
    # Force data frame if tibble (tibbles do not support rownames)
    if(tibble::is_tibble(dataframe)){
      dataframe = as.data.frame(dataframe)
    }
    genecol = colnames(dataframe)[1]

    numdat = dataframe
    rownames(numdat) = dataframe[[genecol]]
    numdat = numdat[, -which(colnames(numdat) == genecol)]
    numdat = as.matrix(numdat)
    numdat = as(numdat, "sparseMatrix")
  } else {
    numdat = dataframe %>% as(., "sparseMatrix")
  }
  return(numdat)
}

##
# expandSparse: takes a sparsed matrix and returns a dataframe
# @param sparsedMatrix, a sparse matrix made from the Matrix
# package
# @return a data frame
#
expandSparse = function(sparsedMatrix){
  NonSparse = data.frame(as.matrix(sparsedMatrix), check.names=F)
  return(NonSparse)
}

##
# sparse_row_sds: per-row standard deviation without densifying.
#
# Replaces `apply(m, 1, sd)`, which is a silent full densification: base
# apply() calls as.matrix() on anything that is not already an array, so a
# genes x cells dgCMatrix becomes a dense double matrix first. That is 16 MB for
# a 500 x 4096 sample and 80 MB for 200 x 50000, allocated per sample and, under
# mclapply, per fork.
#
# Two-pass, matching stats::sd's own approach rather than the one-pass
# sqrt((ss - s^2/n)/(n-1)) shortcut, which can produce NaN from a negative
# radicand on a near-constant row. The structural zeros contribute
# (n - nnz) * mu^2 to the centered sum of squares, added back analytically
# instead of being visited.
#
# NOT bit-identical to apply(m, 1, sd): the summation order differs, giving
# ~1e-16 relative disagreement at small n (it happens to agree exactly at
# 4096 columns). That is acceptable here because gene_stdevs, its only
# consumer, is a reported value with no downstream reader -- `topgenes`
# selection uses vst.variance.standardized from Seurat_FindVariableFeatures(),
# not this -- so the shift cannot cascade into gene selection. Verify that claim
# still holds before reusing this anywhere new.
#
# @param m a dgCMatrix (genes in rows).
# @return a named numeric vector of per-row standard deviations.
#
sparse_row_sds = function(m){
  if(!inherits(m, 'dgCMatrix')){
    m = methods::as(m, 'CsparseMatrix')
  }

  n = ncol(m)
  # sd of a single observation is NA, matching stats::sd.
  if(n < 2){
    return(stats::setNames(rep(NA_real_, nrow(m)), rownames(m)))
  }

  mu = Matrix::rowMeans(m)

  # Centered sum of squares over the stored entries only. m@i holds 0-based row
  # indices of each stored value, so mu[m@i + 1L] aligns each entry with its
  # row mean.
  centered = m
  centered@x = (centered@x - mu[centered@i + 1L])^2
  ss = Matrix::rowSums(centered)

  # ...plus the structural zeros, each contributing mu^2.
  nnz = base::tabulate(m@i + 1L, nbins = nrow(m))
  ss = ss + (n - nnz) * mu^2

  stats::setNames(sqrt(ss / (n - 1)), rownames(m))
}
