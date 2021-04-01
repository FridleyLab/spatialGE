##
# This function takes an STList and a vector with HUGO gene names and returns
# Morans' I for each element of the vector.
#
# @param x, a STList with normalized counts.
# @return x, a STList including the values corresponding to Moran's I for each gene.
#
#
gene_moran_I <- function(x=NULL, genes=NULL, subj=NULL) {

  # Test if voom normalized counts are available.
  if(is_empty(x@voom_counts[[subj]])){
    stop("There are no normalized counts in STList.")
  }

  # Create distance matrix based on the coordinates of each sampled location.
  subj_dists <- as.matrix(dist(x@coords[[subj]][2:3]))
  subj_dists_inv <- 1/subj_dists
  diag(subj_dists_inv) <- 0

  for(gene in genes){
    # Test if gene name exists in normalized count matrix.
    if(!any(x@voom_counts[[subj]][[1]] == gene)){
      stop(paste(gene, "is not a gene in the normalized count matrix."))
    }

    # Extract expression data (voom counts) for a given gene.
    gene_expr <- unlist(x@voom_counts[[subj]][x@voom_counts[[subj]][[1]] == gene, -1])

    # Estimate statistic.
    moran_est <- spdep::moran.test(gene_expr, spdep::mat2listw(subj_dists_inv))

    # Test if list to store spatial heterogeneity statistics, and create one if
    # needed.
    if(!is.null(x@gene_het[[gene]])){
      if(length(x@gene_het[[gene]]) < subj){
        x@gene_het[[gene]][[subj]] <- list(morans_I=NULL,
                                          gearys_C=NULL,
                                          getis_ord_Gi=NULL)
      }
    }else{
      x@gene_het[[gene]] <- list()
      x@gene_het[[gene]][[subj]] <- list(morans_I=NULL,
                                        gearys_C=NULL,
                                        getis_ord_Gi=NULL)
    }

    # Store statistic in object.
    x@gene_het[[gene]][[subj]]$morans_I <- moran_est

  }
  return(x)

}
