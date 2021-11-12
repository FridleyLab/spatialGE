##
#' @title STclust: Detect clusters on ST arrays
#' @description Perform unsupervised spatially-informed clustering on the spots of an
#' spatial array.
#' @details
#' The function takes an STList and performs PCA on the gene expression data, then
#' calculates euclidean distances of the PCs and spatial distances. Weighted averages
#' between the two distances are calculated to perform hierarchical clustering. The
#' user can define how much weight the spatial distance should have (values
#' between 0.05 to 0.25 work reasonably well). The method aims to detect compartments
#' within a tissue.
#'
#' @param x an STList with normalized expression data.
#' @param weights a double [0-1] indicating the weight to be applied to spatial distances.
#' Defaults to w=0.025.
#' @param pcs the number of principal components (PCs) to retain.
#' @param vperc the minimum percentage of explained variance explained by the PCs to
#' decide how many of them to retain.
#' @param method the linkage method applied to hierarchical clustering. It is passed
#' to `hclust` and defaults to 'complete'.
#' @param ks the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied.
#' @param topgenes the number of high spot-to-spot standard deviation to retain before PCA.
#' @param deepSplit, a logical or integer [1-4], to be passed to `cutreeDynamic` and
#' control cluster resolution.
#' @return x, the STList with cluster assignments.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # Using Dynamic Tree Cuts:
#' melanoma <- cluster_STspot(melanoma, ks='dtc', weights=0.1)
#'
#' # Using a range of ks:
#' # melanoma <- cluster_STspot(melanoma, ks=c(2:6), weights=0.1)
#'
#' @export
#
#
STclust = function(x=NULL, weights=0.025, pcs=NULL, vperc=NULL, method='complete', ks='dtc', topgenes=3000, deepSplit=F) {

  require('magrittr')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Check to ensure number of ks is acceptable
  if(is.numeric(ks)){
    if(any(ks < 2)){
      stop('Refusing to generate < 2 clusters.')
    }
  }

  # Check to see what PC filter was requested
  if(is.null(vperc) && is.null(pcs)){
    cat('Using vperc = 0.8\n')
    vperc <- 0.8
  } else if(!is.null(vperc) && !is.null(pcs)){
    cat('Both number of retained PCs and explanined variance specified. Using only specified vperc\n')
    pcs <- NULL
  }

  grp_list = list()
  for(i in 1:length(x@tr_counts)){
    # Get transformed counts
    trcounts_df = x@tr_counts[[i]][order(x@gene_var[[i]][['gene_stdevs']], decreasing=T), ][1:topgenes, ]

    # Convert counts and coordinate data to matrices
    A = t(as.matrix(trcounts_df[, -1]))
    rownames(A) = colnames(trcounts_df[, -1])
    B = as.matrix(x@coords[[i]][,-1])
    rownames(B) = x@coords[[i]][[1]]

    # Perform PCA
    count_pcs = prcomp(A, scale=TRUE)
    count_pcs_var = summary(count_pcs)
    count_pcs_var = as.vector(count_pcs_var$importance[3, ])

    # Define the number of PCs to retain
    if(is.null(pcs)){
      vperc_mask <- count_pcs_var <= vperc
      if(sum(vperc_mask) == 0){
        pcs = 1
      } else{
        pcs = sum(vperc_mask)
      }
    }

    # Use all PCs if maximum variance is lower than vperc
    if(max(count_pcs_var) > vperc){
      pcs = length(count_pcs_var)
    }

    A = as.matrix(count_pcs$x[, 1:pcs])

    # Get distance matrices
    #dA <- wordspace::dist.matrix(A)
    #dA <- wordspace::dist.matrix(A, method='euclidean')
    dA <- wordspace::dist.matrix(A, method='manhattan')
    dB <- dist(B, upper=T, diag=T)
    dAm <- as.matrix(dA)
    dBm <- as.matrix(dB)

    # Scale matrices
    dAm <- dAm/max(dAm)
    dBm <- dBm/max(dBm)

    for(w in 1:length(weights)){
      weight_d <- weights[w]
      weight_g <- 1-weight_d

      # Create vector of weights for Reduce
      weightList <- c(weight_g, weight_d)
      dMxs <- list(dAm, dBm)

      # Apply wieight element-wise
      weight_m <- Reduce(`+`,Map(`*`, dMxs, weightList))
      # COnvert result to matrix
      weight_m <- matrix(weight_m, nrow(dAm))
      weight_m_d <- as.dist(weight_m)

      # Performe hclust
      hierclusters <- hclust(weight_m_d, method=method)

      # Apply dtc or split to k
      if(is.character(ks)){
        if(ks == 'dtc'){
          x@misc[['STclust_cuttype']] <- 'dtc'
          grp_df <- dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weight_m, deepSplit=deepSplit, verbose=F)
          grp_df <- tibble::tibble(colnames(dAm), as.factor(grp_df))
          names(grp_df) <- c('libname', 'WCluster')

          grp_df$WCluster[grp_df$WCluster == 0] <- NA

          grp_df <- dplyr::full_join(x@coords[[i]], grp_df, by='libname')

          grp_list[[paste0(names(x@tr_counts[i]), "_spw", weights[w])]] <- grp_df
        }
      }else if(is.numeric(ks)){
        x@misc[['STclust_cuttype']] <- 'multiK'
        for(k in ks){
          singlek <-  cutree(hierclusters, k=k)
          singlek <- tibble::tibble(colnames(dAm), as.factor(singlek))
          names(singlek) <- c('libname', 'WCluster')
          singlek <- dplyr::full_join(x@coords[[i]], singlek, by='libname')
          grp_list[[paste0(names(x@tr_counts[i]), "_spw", weights[w], '_k', k)]] <- singlek
        }
      } else{
        stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
      }
    }
  }
  x@st_clusters <- grp_list

  return(x)
}
