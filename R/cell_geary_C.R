##
#' @title cell_geary_C
#' @description Calculates Geary's C from cell deconvolution scores in ST data.
#' @details
#' This function takes an STList and a vector with cell names (from xCell, for
#' example) and returns Geary's C for each element of the vector.
#'
#' @param x, an STList with deconvoluted cell data.
#' @param cells, a vector of cell names present in the deconvolution matrix.
#' @param subj, an integer indicating the spatial array for which the statistic
#' will be calculated.
#' @param method, the deconvolutionn method to used estimate the cell scores.
#' @return x, a STList including the values corresponding to Geary's C for each
#' cell in cells.
#' @export
#
#
cell_geary_C <- function(x=NULL, cells=NULL, subj=NULL, method='xcell') {

  method=tolower(method)

  # Get requested list of deconvoluted matrices.
  if(method == 'xcell'){
    deconv_list <- x@cell_deconv$xCell
  } else if(method == 'ssgsea'){
    deconv_list <- x@cell_deconv$ssGSEA
  } else{
    stop('Please, specify a deconvolution method to plot.')
  }

  # Test if normalized cell data are available.
  # if(rlang::is_empty(x@cell_deconv$xCell[[subj]]$sqrt_scores)){
  if(rlang::is_empty(deconv_list[[subj]]$sqrt_scores)){
    stop(paste0("There are no ", method, " normalized cell data in STList."))
  }

  # Create distance matrix based on the coordinates of each sampled location.
  subj_dists <- as.matrix(dist(x@coords[[subj]][2:3]))
  subj_dists_inv <- 1/subj_dists
  diag(subj_dists_inv) <- 0

  for(cell in cells){
    # Test if cell name exists in normalized data.
    #if(!any(x@cell_deconv$xCell[[subj]]$sqrt_scores[[1]] == cell)){
    if(!any(deconv_list[[subj]]$sqrt_scores[[1]] == cell)){
      stop(paste(cell, "is not a cell name in the normalized data."))
    }

    # Extract cell data (deconvoluted matrix) for a given cell.
    # cell_data <- unlist(x@cell_deconv$xCell[[subj]]$sqrt_scores[x@cell_deconv$xCell[[subj]]$sqrt_scores[[1]] == cell, -1])
    cell_data <- unlist(deconv_list[[subj]]$sqrt_scores[deconv_list[[subj]]$sqrt_scores[[1]] == cell, -1])

    # Estimate statistic.
    geary_est <- spdep::geary.test(cell_data, spdep::mat2listw(subj_dists_inv))

    # Test if list to store spatial heterogeneity statistics, and create one if
    # needed.
    if(!is.null(x@cell_het[[cell]])){
      if(length(x@cell_het[[cell]]) < subj){
        x@cell_het[[cell]][[subj]] <- list(morans_I=NULL,
                                           gearys_C=NULL,
                                           getis_ord_Gi=NULL)
      }
    }else{
      x@cell_het[[cell]] <- list()
      x@cell_het[[cell]][[subj]] <- list(morans_I=NULL,
                                         gearys_C=NULL,
                                         getis_ord_Gi=NULL)
    }

    # Store statistic in object.
    x@cell_het[[cell]][[subj]]$gearys_C <- geary_est

  }
  return(x)

}
