##
#' @title cell_getis_Gi
#' @description Calculates Getis-Ord Gi C from cell deconvolution scores in ST data.
#' @details
#' This function takes an STList and a vector with cell names (from xCell, for
#' example) and returns Getis-Ord Gi for each element of the vector.
#'
#' @param x, an STList with deconvoluted cell data.
#' @param cells, a vector of cell names present in the deconvolution matrix.
#' @param subj, an integer indicating the spatial array for which the statistic
#' will be calculated.
#' @return x, a STList including the values corresponding to Getis-Ord Gi for each
#' cell in cells.
#' @export
#
#
cell_getis_Gi <- function(x=NULL, cells=NULL, subj=NULL) {

  # Test if normalized cell data are available.
  if(is_empty(x@cell_deconv$xCell[[subj]]$sqrt_scores)){
    stop("There are no normalized cell data in STList.")
  }

  # Create distance matrix based on the coordinates of each sampled location.
  subj_dists <- as.matrix(dist(x@coords[[subj]][2:3]))
  subj_dists_inv <- 1/subj_dists
  diag(subj_dists_inv) <- 0

  for(cell in cells){
    # Test if cell name exists in normalized data.
    if(!any(x@cell_deconv$xCell[[subj]]$sqrt_scores[[1]] == cell)){
      stop(paste(cell, "is not a cell name in the normalized data."))
    }

    # Extract cell data (deconvoluted matrix) for a given cell.
    cell_data <- unlist(x@cell_deconv$xCell[[subj]]$sqrt_scores[x@cell_deconv$xCell[[subj]]$sqrt_scores[[1]] == cell, -1])

    # Estimate statistic.
    getis_est <- spdep::globalG.test(cell_data, spdep::mat2listw(subj_dists_inv, style='B'))

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
    x@cell_het[[cell]][[subj]]$getis_ord_Gi <- getis_est

  }
  return(x)

}
