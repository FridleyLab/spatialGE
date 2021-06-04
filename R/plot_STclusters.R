##
#' @title plot_STclusters: Plot cluster memberships of ST spots
#' @description Generates a plot of the location of spots within an spatial array,
#' and colors them according to spatially-weighted hierarchical clustering assignments.
#' @details
#' The function takes an STList with cluster memberships and plots the spots with
#' colors indicating the cluster they belong to. Optionally, the user can annotate
#' tumor/stroma compartments if ESTIMATE scores are available.
#'
#' @param x, an STList with hierarchical cluster memberships.
#' @param who, an integer indicating the spatial array(s) to be plotted.
#' @param purity, logical, whether or not annotate tumor/stroma.
#' @param color_pal, a scheme from 'khroma'. Default is 'light'.
#' @export
#
#
plot_STclusters <- function(x, plot_who=NULL, purity=F, color_pal='light'){

  require(ggplot2)

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  if(is.null(plot_who)){
    subjs <- (1:length(x@st_clusters$clust_dfs))
  } else{
    subjs <- plot_who
  }

  plot_list <- list()
  if(x@st_clusters$type == 'dtc'){
    for(i in subjs){

      df <- x@st_clusters$clust_dfs[[i]]

      if(purity){

        # Test ESTIMATE clusters are available.
        if(is.null(x@cell_deconv$ESTIMATE)){
          stop('No ESTIMATE cluster annotations available')
        }
        purity_df <- x@cell_deconv$ESTIMATE[[i]]$purity_clusters
        df <- dplyr::full_join(df, purity_df, by='X1')
        df$cluster <- as.factor(df$cluster)
        names(df)[5] <- "EstCluster"
      }

      nas <- complete.cases(df$WCluster)
      df_nonas <- df[nas, ]
      df$WCluster <- as.factor(tidyr::replace_na(as.vector(df$WCluster), "no_cluster"))
      p <- ggplot()
      if(ncol(df) == 4){
        p <- p + geom_point(data=df, aes(x=X3, y=X2, color=WCluster), size=0.7, shape=19)
        p <- p + geom_point(data=df_nonas, aes(x=X3, y=X2, color=WCluster), size=0.5)
        p <- p + scale_color_manual(values=c(as.vector(p_palette(max(as.numeric(levels(df_nonas$WCluster))))), 'gray50'))
        p <- p + labs(title=paste0("Hier. Clusters (dynamicTreeCut), subj ", i), color='HClusters')
        p <- p + guides(color=guide_legend(override.aes=list(size=2)))
      } else{
        p <- p + geom_point(data=df, aes(x=X3, y=X2, color=WCluster), size=0.5, shape=19)
        p <- p + geom_point(data=df_nonas, aes(x=X3, y=X2, color=WCluster, shape=EstCluster), size=1.2)
        p <- p + scale_color_manual(values=c(as.vector(p_palette(max(as.numeric(levels(df_nonas$WCluster))))), 'gray50'))
        p <- p + scale_shape_manual(values=c(3, 15))
        p <- p + labs(title=paste0("Hier. Clusters (dynamicTreeCut), subj ", i), color='HClusters', shape='tumor/stroma')
        p <- p + guides(shape=guide_legend(override.aes=list(size=2)), color=guide_legend(override.aes=list(size=2)))
      }

      p <- p + ylab('Y Position') + xlab('X Position')
      p <- p + scale_x_reverse() + scale_y_reverse() + coord_fixed() + theme_classic()
      p <- p + theme(plot.title=element_text(size=10), legend.text=element_text(size=10))

      plot_list[[paste0("p",i)]] <- p
    }

  } else{
    for(i in subjs){
      plot_list[[paste0("p",i)]] <- list()
      for(k in names(x@st_clusters$clust_dfs[[i]])){

        df <- x@st_clusters$clust_dfs[[i]][[k]]

        if(purity){

          # Test ESTIMATE clusters are available.
          if(is.null(x@cell_deconv$ESTIMATE)){
            stop('No ESTIMATE cluster annotations available')
          }
          purity_df <- x@cell_deconv$ESTIMATE[[i]]$purity_clusters
          df <- dplyr::full_join(df, purity_df, by='X1')
          df$cluster <- as.factor(df$cluster)
          names(df)[5] <- "EstCluster"
        }

        nas <- complete.cases(df$WCluster)
        df_nonas <- df[nas, ]
        df$WCluster <- as.factor(tidyr::replace_na(as.vector(df$WCluster), "no_cluster"))
        p <- ggplot()
        if(ncol(df) == 4){
          p <- p + geom_point(data=df, aes(x=X3, y=X2, color=WCluster), size=0.7, shape=19)
          p <- p + geom_point(data=df_nonas, aes(x=X3, y=X2, color=WCluster), size=0.5)
          p <- p + scale_color_manual(values=c(as.vector(p_palette(max(as.numeric(levels(df_nonas$WCluster))))), 'gray50'))
          p <- p + labs(title=paste0("Hier. Clusters ", k, ", subj ", i), color='HClusters')
          p <- p + guides(color=guide_legend(override.aes=list(size=2)))
        } else{
          p <- p + geom_point(data=df, aes(x=X3, y=X2, color=WCluster), size=0.5, shape=19)
          p <- p + geom_point(data=df_nonas, aes(x=X3, y=X2, color=WCluster, shape=EstCluster), size=1.2)
          p <- p + scale_color_manual(values=c(as.vector(p_palette(max(as.numeric(levels(df_nonas$WCluster))))), 'gray50'))
          p <- p + scale_shape_manual(values=c(3, 15))
          p <- p + labs(title=paste0("Hier. Clusters ", k, ", subj ", i), color='HClusters', shape='tumor/stroma')
          p <- p + guides(shape=guide_legend(override.aes=list(size=2)), color=guide_legend(override.aes=list(size=2)))
        }
        p <- p + ylab('Y Position') + xlab('X Position')
        p <- p + scale_x_reverse() + scale_y_reverse() + coord_fixed() + theme_classic()
        p <- p + theme(plot.title=element_text(size=10), legend.text=element_text(size=10))

        plot_list[[paste0("p",i)]][[k]] <- p
      }
    }
  }

  return(plot_list)

}
