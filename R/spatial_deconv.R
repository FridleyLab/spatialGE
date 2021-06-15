##
#' @title spatial_deconv: Deconvolution of ST data
#' @description Applies gene expression deconvolution methods to data from spatial
#' arrays. Produces ESTIMATE tumor/stroma classes via model-based clustering.
#' @details
#' This function applies deconvolution methods (xCell and ESTIMATE, for now) to the stored
#' normalized matrices in order to obtain cell scores for each of the library/spots.
#' The results are stored witin the STList. The method ESTIMATE is automatically
#' applied regardless of the selected method.
#'
#' @param x, an STList with normalized count matrices. If 'none', only ESTIMATE is
#' applied.
#' @param method, the name of the method to deconvolute data. As of now, only xCell is supported.
#' @return x, an updated STList with deconvolution scores.
#' @export
#
#
spatial_deconv <- function(x=NULL, method='xcell'){

  method=tolower(method)

  if(is.null(x@cell_deconv[['ESTIMATE']])){
    cat('This will take some time...\n')
    x <- spatial_purity(x)
    x <- cluster_purity(x)
  }

  if(method == 'none' | method == 'estimate'){
    cat('Only ESTIMATE was applied to the spatial arrays.')
  }

  if(method == 'xcell'){
    require('xCell') # Needs to be 'required' because of databases loaded by packages.
    x <- spatial_xcell(x)
  }

  return(x)

}
