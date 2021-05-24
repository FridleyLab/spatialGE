##
#' @title cluster_purity
#' @description Perform model-based clustering to define stromal/tumor clusters from
#' ESTIMATE purity scores.
#' @details
#' Takes ESTIMATE tumor purity scores and performs model-based clustering. The function
#' only evaluates k=2, to separate between likely tumor and stromal spots.
#'
#' @param x, an STList.
#' @param who, an integer indicating the spatial array to be analyzed.
#
#
cluster_purity <- function(x=NULL) {

  require('mclust') #Mclust can't call mclust if the library is not loaded.

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Iterate over spatial arrays.
  for(i in 1:length(x@cell_deconv$ESTIMATE)){

    # Extract purity scores.
    estimate_df <- x@cell_deconv$ESTIMATE[[i]]$estimate_purity

    # Transpose values to cluster.
    estimate_mx <- t(as.matrix(estimate_df[1, -1]))

    # Perform clustering.
    BIC <- mclust::mclustBIC(estimate_mx, G=c(2))
    clust_mod <- mclust::Mclust(estimate_mx, x=BIC, G=c(2))

    #factoextra::fviz_mclust(clust_mod, "BIC", palette = "jco")

    # Store tumor/stroma classifications.
    clusters <- tibble::as_tibble_col(clust_mod$classification, column_name='cluster')
    clusters$cluster <- gsub('1', 'tumor', clusters$cluster)
    clusters$cluster <- gsub('2', 'stroma', clusters$cluster)
    clusters <- tibble::add_column(clusters, X1=names(clust_mod$classification), .before=1)

    x@cell_deconv$ESTIMATE[[i]][['purity_clusters']] <- clusters

  }

  return(x)

}

