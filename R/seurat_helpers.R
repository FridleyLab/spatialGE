##
# Utility and helper functions from Seurat
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
# @title FindVariableFeatures
# @description Extracted from the Seurat package
#

Seurat_FindVariableFeatures = function(object=NULL, verbose=F){
  # Fixed arguments
  selection.method="vst"
  loess.span=0.3
  clip.max='auto'
  #mean.function=FastExpMean
  #dispersion.function=FastLogVMR
  #num.bin=20
  #binning.method="equal_width"

  if (!inherits(x=object, 'Matrix')) {
    object <- as(object = as.matrix(x = object), Class = 'Matrix')
  }
  if (!inherits(x = object, what = 'dgCMatrix')) {
    object <- as.sparse(x = object)
  }
  if (selection.method == "vst") {
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = object))
    }
    hvf.info <- data.frame(mean = rowMeans(x = object))
    hvf.info$variance <- SparseRowVar2(
      mat = object,
      mu = hvf.info$mean,
      display_progress = verbose
    )
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
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
      vmax = clip.max,
      display_progress = verbose
    )
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }
  #  else {
  #   if (!inherits(x = mean.function, what = 'function')) {
  #     stop("'mean.function' must be a function")
  #   }
  #   if (!inherits(x = dispersion.function, what = 'function')) {
  #     stop("'dispersion.function' must be a function")
  #   }
  #   feature.mean <- mean.function(object, verbose)
  #   feature.dispersion <- dispersion.function(object, verbose)
  #   names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = object)
  #   feature.dispersion[is.na(x = feature.dispersion)] <- 0
  #   feature.mean[is.na(x = feature.mean)] <- 0
  #   data.x.breaks <- switch(
  #     EXPR = binning.method,
  #     'equal_width' = num.bin,
  #     'equal_frequency' = c(
  #       -1,
  #       quantile(
  #         x = feature.mean[feature.mean > 0],
  #         probs = seq.int(from = 0, to = 1, length.out = num.bin)
  #       )
  #     ),
  #     stop("Unknown binning method: ", binning.method)
  #   )
  #   data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
  #   names(x = data.x.bin) <- names(x = feature.mean)
  #   mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
  #   sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
  #   feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
  #     sd.y[as.numeric(x = data.x.bin)]
  #   names(x = feature.dispersion.scaled) <- names(x = feature.mean)
  #   hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
  #   rownames(x = hvf.info) <- rownames(x = object)
  #   colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
  # }
  return(hvf.info)
}

##
# @title SparseRowVar2
# @description Extracted from the Seurat package
#

SparseRowVar2 = function(mat, mu, display_progress) {
  .Call('_Seurat_SparseRowVar2', PACKAGE='Seurat', mat, mu, display_progress)
}

##
# @title SparseRowVarStd
# @description Extracted from the Seurat package
#
SparseRowVarStd <- function(mat, mu, sd, vmax, display_progress) {
    .Call('_Seurat_SparseRowVarStd', PACKAGE = 'Seurat', mat, mu, sd, vmax, display_progress)
}

