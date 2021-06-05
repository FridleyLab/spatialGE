##
#' @title cluster_STspot: Detect spot clusters on ST arrays
#' @description Perform spatially-informed hierarchical clustering of spots within a
#' spatial array using a genetic distance matrix weighted by spatial distances.
#' @details
#' The function takes an STList and calculates euclidean distances of normalized
#' gene expression and spatial distances. Then, it calculates the weighted average
#' between the two distances to perform hierarchical clustering. The user can define
#' how much weight the spatial distance should have (values between 0.1 to 025
#' work reasonably well).
#'
#' @param x, an STList with normalized expression data.
#' @param weights, a double [0-1] indicating the weight to be applied to spatial
#' distances.
#' @param method, the linkage method applied to hierarchical clustering. It is passed
#' to `hclust` and defaults to 'ward.D'.
#' @param ks, the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied.
#' @param spotfilter, the number of genes with more than zero counts for a spot to be
#' included in the clustering analysis. The lower this number, the slower the analysis.
#' @export
#
#
cluster_STspot <- function(x=NULL, weights=0.1, method='ward.D', ks='dtc', spotfilter=500) {

  require('magrittr')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  grp_list <- list()
  grp_list[['clust_dfs']] <- list()
  for(i in 1:length(x@counts)){

    counts_df <- x@counts[[i]]
    voom_df <- x@voom_counts[[i]][x@voom_counts[[i]]$gene %in% xCell::xCell.data$genes, ]

    spotlibs <- c()
    for(dfcol in 2:ncol(counts_df)){
      if(sum(counts_df[, dfcol] > 0) > spotfilter){
        spotlibs <- append(spotlibs, names(counts_df)[dfcol])
      }
    }

    voom_df <- voom_df[, names(voom_df) %in% spotlibs]
    coords_df <- x@coords[[i]][x@coords[[i]]$X1 %in% spotlibs, ]
    #if(purity){
    #  purity_df <- x@cell_deconv$ESTIMATE[[i]]$purity_clusters[x@cell_deconv$ESTIMATE[[i]]$purity_clusters$X1 %in% spotlibs, ]
    #}

    A <- t(as.matrix(voom_df))
    rownames(A) <- spotlibs

    B <- as.matrix(coords_df[,-1])
    rownames(B) <- spotlibs

    dA <- wordspace::dist.matrix(A)
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
        grp_list['type'] <- 'dtc'
        grp_df <- dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weight_m, deepSplit=F, verbose=F)
        grp_df <- tibble::tibble(colnames(dAm), as.factor(grp_df))
        names(grp_df) <- c('X1', 'WCluster')
        grp_df <- dplyr::full_join(x@coords[[i]], grp_df, by='X1')
        #grp_df <- dplyr::full_join(grp_df, coords_abs_df, by='X1')
        #if(purity){
        #  grp_df <- dplyr::full_join(grp_df, purity_df, by='X1')
        #  grp_df$cluster <- as.factor(grp_df$cluster)
        #  names(grp_df)[5] <- "EstCluster"
        #}
        grp_list$clust_dfs[[i]] <- grp_df
        }
      }else if(is.numeric(ks)){
        grp_list['type'] <- 'multiK'
        grp_list$clust_dfs[[i]] <- list()
        for(k in ks){
          singlek <-  cutree(hierclusters, k=k)
          singlek <- tibble::tibble(colnames(dAm), as.factor(singlek))
          names(singlek) <- c('X1', 'WCluster')
          singlek <- dplyr::full_join(x@coords[[i]], singlek, by='X1')
          #singlek <- dplyr::full_join(singlek, coords_abs_df, by='X1')
          #if(purity){
          #  singlek <- dplyr::full_join(singlek, purity_df, by='X1')
          #  singlek$cluster <- as.factor(singlek$cluster)
          #  names(singlek)[5] <- "EstCluster"
          #}
          grp_list$clust_dfs[[i]][[paste0('k', k)]] <- singlek
        }
      } else{
        stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
      }

    }

  }
  x@st_clusters <- grp_list

  return(x)
}
