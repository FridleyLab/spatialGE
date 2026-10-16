##
# UTILITY AND HELPER FUNCTIONS MODIFIED FROM THE SEURAT PACKAGE
#
# Copyright (c) 2021 Seurat authors
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons
# to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

##
# @title Seurat_FindVariableFeatures
# @description Extracted from the Seurat package
#
#' @importFrom stats loess setNames
#
#
Seurat_FindVariableFeatures = function(object=NULL){
  # Fixed arguments
  selection.method="vst"
  loess.span=0.3
  clip.max='auto'

  # Coerce to a general, double, CSC matrix. SparseRowVar2()/SparseRowVarStd()
  # take an Eigen::SparseMatrix<double> (src/seurat_helpers.cpp), so RcppEigen
  # needs exactly that -- a dense dgeMatrix, a logical lgCMatrix and a symmetric
  # dsCMatrix are all unacceptable.
  #
  # PROVENANCE. This file is vendored Seurat code, and the original line here was
  # `object <- as.sparse(x = object)`, i.e. SeuratObject::as.sparse(). That
  # function exists and is exported by SeuratObject, but SeuratObject is only in
  # Suggests and is never imported, so the call was never resolvable from this
  # package's namespace -- the vendored code came across without its dependency.
  # It was not lost between 1.2.2 and 2.0.0: `main` carries the identical broken
  # call, and as.sparse has never been defined anywhere in this repository. Any
  # input that was not already a dgCMatrix failed with 'could not find function
  # "as.sparse"'. It stayed invisible because every current caller happens to
  # pass a dgCMatrix, and because 'as.sparse' was listed in
  # utils::globalVariables() in R/globals.R, suppressing the one R CMD check NOTE
  # that would have flagged it.
  #
  # Rather than take a hard SeuratObject dependency for this, the coercion is
  # reproduced with Matrix. SeuratObject::as.sparse() is
  # as(as(as(x, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix'), and both inner
  # steps are load-bearing: 'dMatrix' forces double (a logical sparse input
  # otherwise stays lgCMatrix) and 'generalMatrix' forces general (a symmetric
  # input otherwise stays dsCMatrix). The added as.matrix() branch covers
  # data.frame input, which SeuratObject's chain errors on.
  if(!inherits(x=object, what='dgCMatrix')){
    if(!inherits(x=object, what='Matrix')){
      object <- as.matrix(x = object)
    }
    object <- methods::as(methods::as(methods::as(object, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
  }
  if (selection.method == "vst") {
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = object))
    }
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0

    # The loess below is fit across GENES, on the non-constant ones only, so it
    # has a minimum data requirement that the rest of this function does not. At
    # span=0.3 and loess' default degree=2 it needs 4 such genes; with 1-3 it
    # stops with 'span is too small', an internal stats message naming neither
    # the sample nor the cause. Note the count that matters is of non-constant
    # genes, not of genes: a 100-gene matrix in which only 3 vary fails too.
    #
    # Returning NULL rather than raising lets calculate_vst() decide, because the
    # right answer differs by caller -- transform_data() can still deliver the
    # normalization a caller actually asked for, while a topgenes selection
    # genuinely cannot proceed.
    if(sum(not.const) < 4){
      return(NULL)
    }

    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    # use c function to get variance after feature standardization
    hvf.info$variance.standardized <- SparseRowVarStd(
      mat = object,
      mu = hvf.info$mean,
      sd = sqrt(hvf.info$variance.expected),
      vmax = clip.max
    )
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }
  return(hvf.info)
}

