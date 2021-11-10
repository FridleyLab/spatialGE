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
#' Defaults to w=0.1.
#' @param pcs the number of principal components (PCs) to retain.
#' @param vperc the minimum percentage of explained variance explained by the PCs to
#' decide how many of them to retain.
#' @param method the linkage method applied to hierarchical clustering. It is passed
#' to `hclust` and defaults to 'ward.D'.
#' @param ks the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied.
#' @param spotfilter the number of genes with more than zero counts for a spot to be
#' included in the clustering analysis. The lower this number, the slower the analysis. If
#' spotfilter=0, all spots are retained.
#' @param topgenes the number of high spot-to-spot standard deviation to retain before PCA.
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
STclust <- function(x=NULL, weights=0.1, pcs=NULL, vperc=NULL, method='ward.D', ks='dtc', spotfilter=0, topgenes=2000) {

  require('magrittr')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  if(is.numeric(ks)){
    if(any(ks < 2)){
      stop('Refusing to generate < 2 clusters.')
    }
  }


  if(is.null(vperc) && is.null(pcs)){
    cat('Using vperc = 0.8\n')
    vperc <- 0.8
  }

  if(!is.null(vperc) && !is.null(pcs)){
    cat('Both number of retained PCs and explanined variance specified. Using specified vperc\n')
    vperc <- 0.8
    pcs <- NULL
  }

  grp_list <- list()
  grp_list[['clust_dfs']] <- list()
  for(i in 1:length(x@counts)){

    counts_df <- x@counts[[i]]
    #voom_df <- x@voom_counts[[i]][x@voom_counts[[i]]$gene %in% xCell::xCell.data$genes, ]
    voom_df <- x@voom_counts[[i]][order(x@gene_stdev[[i]]$gene_stdevs, decreasing=T), ][1:topgenes, ]
    #voom_df <- x@voom_counts[[i]]

    spotlibs <- c()
    for(dfcol in 2:ncol(counts_df)){
      if(sum(counts_df[, dfcol] > 0) > spotfilter){
        spotlibs <- append(spotlibs, names(counts_df)[dfcol])
      }
    }

    voom_df <- voom_df[, names(voom_df) %in% spotlibs]
    coords_df <- x@coords[[i]][x@coords[[i]]$libname %in% spotlibs, ]
    #if(purity){
    #  purity_df <- x@cell_deconv$ESTIMATE[[i]]$purity_clusters[x@cell_deconv$ESTIMATE[[i]]$purity_clusters$X1 %in% spotlibs, ]
    #}

    A <- t(as.matrix(voom_df))
    rownames(A) <- spotlibs

    B <- as.matrix(coords_df[,-1])
    rownames(B) <- spotlibs

    voom_pcs <- prcomp(A, scale=TRUE)
    voom_pcs_var <- summary(voom_pcs)
    voom_pcs_var <- as.vector(voom_pcs_var$importance[3, ])

    if(!is.null(vperc)){
      vperc_mask <- voom_pcs_var < vperc
      pcs <- sum(vperc_mask) + 1
    }

    A <- as.matrix(voom_pcs$x[, 1:pcs])

    #dA <- wordspace::dist.matrix(A)
    #dA <- wordspace::dist.matrix(A, method='euclidean')
    dA <- wordspace::dist.matrix(A, method='manhattan')
    dB <- dist(B, upper = T, diag = T)

    dAm <- as.matrix(dA)
    dBm <- as.matrix(dB)

    # Scale matrices
    dAm <- dAm/max(dAm)
    dBm <- dBm/max(dBm)

    for(w in 1:length(weights)){

      weight_d <- weights[w]
      weight_g <- 1-weight_d

      weightList <- c(weight_g, weight_d)
      dMxs <- list(dAm, dBm)

      weight_m <- Reduce(`+`,Map(`*`, dMxs, weightList))

      # Old, slow element-wise weighted average calculation.
      # weight_m2 <- c()
      # for(i in 1:nrow(dAm)) {
      #   for(j in 1:ncol(dAm)) {
      #     weight_m2 <- append(weight_m2, (weight_g*dAm[i, j]) + (weight_d*dBm[i, j]))
      #   }
      # }
      #

      weight_m <- matrix(weight_m, nrow(dAm))
      weight_m_d <- as.dist(weight_m)

      hierclusters <- hclust(weight_m_d, method=method)

      # Coordinates not in clustering analysis
      #coords_abs_df <- x@coords[[i]][!(x@coords[[i]]$X1 %in% spotlibs), ]


      if(is.character(ks)){
        if(ks == 'dtc'){
        grp_list[['type']] <- 'dtc'
        grp_df <- dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weight_m, deepSplit=F, verbose=F)
        grp_df <- tibble::tibble(colnames(dAm), as.factor(grp_df))
        names(grp_df) <- c('libname', 'WCluster')

        grp_df$WCluster[grp_df$WCluster == 0] <- NA

        grp_df <- dplyr::full_join(x@coords[[i]], grp_df, by='libname')
        #grp_df <- dplyr::full_join(grp_df, coords_abs_df, by='X1')
        #if(purity){
        #  grp_df <- dplyr::full_join(grp_df, purity_df, by='X1')
        #  grp_df$cluster <- as.factor(grp_df$cluster)
        #  names(grp_df)[5] <- "EstCluster"
        #}
        grp_list$clust_dfs[[paste0("sub", i, "_spw", weights[w])]] <- grp_df
        }
      }else if(is.numeric(ks)){
        grp_list[['type']] <- 'multiK'
        #grp_list$clust_dfs[[paste0("sub", i, "_spw", weights[w])]] <- list()
        for(k in ks){
          singlek <-  cutree(hierclusters, k=k)
          singlek <- tibble::tibble(colnames(dAm), as.factor(singlek))
          names(singlek) <- c('libname', 'WCluster')
          singlek <- dplyr::full_join(x@coords[[i]], singlek, by='libname')
          #singlek <- dplyr::full_join(singlek, coords_abs_df, by='X1')
          #if(purity){
          #  singlek <- dplyr::full_join(singlek, purity_df, by='X1')
          #  singlek$cluster <- as.factor(singlek$cluster)
          #  names(singlek)[5] <- "EstCluster"
          #}
          #grp_list$clust_dfs[[paste0("sub", i, "_spw", weights[w])]][[paste0('k', k)]] <- singlek
          grp_list$clust_dfs[[paste0("sub", i, "_spw", weights[w], '_k', k)]] <- singlek
        }
      } else{
        stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
      }

    }

  }
  x@st_clusters <- grp_list

  return(x)
}
