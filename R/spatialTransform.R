##
#' @title spatialTransform:
#' @description Applies
#' @details
#' This function takes
#'
#'
#'
#'
#'
#'
#' The function works on parallel using "forking" (not in Windows OS).
#'
#' @param x an STList with raw count matrices.
#' @param method one of 'log', 'logtmm', or 'voom'
#' @param cores integer, the number of cores to use in calculations
#' @return x an updated STList with transformed counts.
#'
#' @examples
#' #
#' #
#'
#' @export
#
#
spatialTransform = function(x=NULL, method='log', scale_f=1000000, cores=1){

  # Detect transformation method
  if(method == 'voom'){
    x = voom_norm(x)
  } else if(method == 'log'){
    x = log_transf(x, scale_f=scale_f, cores=cores)
  }  else if(method == 'logtmm'){
    x = log_tmm_transf(x, cores=cores)
  }

  return(x)

}

