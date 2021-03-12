##
# This function takes an STList and a vector with cell names (from xCell, for example)
# and returns Geary's C for each element of the vector.
#
# @param x, a STList with cell data.
# @return x, a STList including the values corresponding to Geary's C for each cell
#
#
# Load 'tidyverse' for tibble manipulation.
# Load 'spdep' for estimation of spatial heterogeneity.
require('tidyverse')
require('spdep')

cell_geary_C <- function(x=NULL, cells=NULL) {

  # Test if normalized cell data are available.
  if(is_empty(x@cell_deconv$transf_deconv_matrix)){
    stop("There are no normalized cell data in STList.")
  }

  # Create distance matrix based on the coordinates of each sampled location.
  subj_dists <- as.matrix(dist(x@coords[2:3]))
  subj_dists_inv <- 1/subj_dists
  diag(subj_dists_inv) <- 0

  for(cell in cells){
    # Test if cell name exists in normalized data.
    if(!any(x@cell_deconv$transf_deconv_matrix[[1]] == cell)){
      stop(paste(gene, "is not a cell name in the normalized data."))
    }

    # Extract cell data (deconvoluted matrix) for a given cell.
    cell_data <- unlist(x@cell_deconv$transf_deconv_matrix[
      x@cell_deconv$transf_deconv_matrix[[1]] == cell, -1])

    # Estimate statistic.
    geary_est <- geary.test(cell_data, mat2listw(subj_dists_inv))

    # Test if list to store spatial heterogeneity statistics, and create one if
    # needed.
    if(is.null(x@cell_het[[cell]])){
      x@cell_het[[cell]] <- list(morans_I=NULL,
                                 gearys_C=NULL,
                                 getis_ord_Gi=NULL)
    }

    # Store statistic in object.
    x@cell_het[[cell]]$gearys_C <- geary_est

  }
  return(x)

}
