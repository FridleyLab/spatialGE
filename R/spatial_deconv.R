##
#' @title spatial_deconv
#' @description Applies deconvolution methods to ST data.
#' @details
#' This function applies deconvolution methods (xCell by default) to the stored
#' normalized matrices in order to obtain cell scores for each of the library/spots.
#' The results are stored as untransformed, and square-root transformed scores.
#' The method ESTIMATE is automatically applied regardless of the selected method.
#'
#' @param x, a STList with normalized count matrices. If 'none', only ESTIMATE is
#' applied.
#' @return x, an updated STList with deconvolution scores.
#' @export
#
#
spatial_deconv <- function(x=NULL, method='xcell'){

  if(is.null(x@cell_deconv[['ESTIMATE']])){
    cat('This will take some time...\n')
    x <- spatial_purity(x)
    x <- cluster_purity(x)
  }

  if(method == 'none'){
    cat('Only ESTIMATE was applied to the spatial arrays.')
  }

  if(method == 'xcell'){
    require('xCell') # Needs to be 'required' because of databases loaded by packages.
    x <- spatial_xcell(x)
  }

  return(x)

}
