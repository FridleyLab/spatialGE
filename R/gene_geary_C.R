##
# This function takes an STList and a vector with HUGO gene names and returns
# Geary's C for each element of the vector.
#
# @param x, a STList with normalized counts.
# @return x, a STList including the values corresponding to Gearys's C for each gene.
#
#
# Load 'tidyverse' for tibble manipulation.
# Load 'spdep' for estimation of spatial heterogeneity.
require('tidyverse')
require('spdep')

gene_geary_C <- function(x=NULL, genes=NULL) {

  # Test if voom normalized counts are available.
  if(is_empty(x@voom_counts)){
    stop("There are no normalized counts in STList.")
  }

  # Create distance matrix based on the coordinates of each sampled location.
  subj_dists <- as.matrix(dist(x@coords[2:3]))
  subj_dists_inv <- 1/subj_dists
  diag(subj_dists_inv) <- 0

  for(gene in genes){
    # Test if gene name exists in normalized count matrix.
    if(!any(x@voom_counts[[1]] == gene)){
      stop(paste(gene, "is not a gene in the normalized count matrix."))
    }

    # Extract expression data (voom counts) for a given gene.
    gene_expr <- unlist(x@voom_counts[x@voom_counts[[1]] == gene, -1])

    # Estimate statistic.
    geary_est <- geary.test(gene_expr, mat2listw(subj_dists_inv))

    # Test if list to store spatial heterogeneity statistics, and create one if
    # needed.
    if(is.null(x@gene_het[[gene]])){
      x@gene_het[[gene]] <- list(morans_I=NULL,
                                 gearys_C=NULL,
                                 getis_ord_Gi=NULL)
    }

    # Store statistic in object.
    x@gene_het[[gene]]$gearys_C <- geary_est

  }
  return(x)

}
