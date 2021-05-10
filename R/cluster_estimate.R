##
#' @title cluster_purity
#' @description Perform model-based clustering (using BIC) to define stromal/tumor
#' clusters from ESTIMATE purity scores.
#' @details
#' yada yada yada
#'
#' @param x, an STList.
#' @param who, an integer indicating the spatial array to be analyzed.
#' @export
#
#
cluster_purity <- function(x=NULL) {

  require('mclust')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

# plot_list <- list()

  for(i in 1:length(x@cell_deconv$ESTIMATE)){

#    for(j in 2:3){
      estimate_df <- x@cell_deconv$ESTIMATE[[i]]$estimate_purity

      estimate_mx <- t(as.matrix(estimate_df[1, -1]))

      BIC <- mclust::mclustBIC(estimate_mx, G=c(2))

      clust_mod <- mclust::Mclust(estimate_mx, x=BIC, G=c(2))

      #factoextra::fviz_mclust(clust_mod, "BIC", palette = "jco")

      clusters <- tibble::as_tibble_col(clust_mod$classification, column_name='cluster')
      clusters <- tibble::add_column(clusters, X1=names(clust_mod$classification), .before=1)
#      clusters <- dplyr::full_join(x@coords[[i]], clusters, by='X1')

      x@cell_deconv$ESTIMATE[[i]][['purity_clusters']] <- clusters

      # p1 <-   ggplot(clusters) +
      #   geom_point(data=x@coords[[i]], aes(x=X2, y=X3), shape=1) +
      #   geom_point(aes(x=X2, y=X3, col=as.factor(cluster))) +
      #   labs(col='cluster', title = paste0('subj_', i)) +
      #   coord_fixed() +
      #   theme_classic()
      #
      # plot_list[[paste0('subj_', i, "_clust_", j)]] <- p1

#   }

  }

  # pdf("~/Desktop/modclust_purity.pdf")
  #   ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow=2)
  # dev.off()
  #

  return(x)

}






