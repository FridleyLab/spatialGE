##
#' @title plot_STclusters: Plot cluster memberships of ST spots
#' @description Generates a plot of the location of spots within an spatial array,
#' and colors them according to spatially-weighted hierarchical clustering assignments.
#' @details
#' The function takes an STList with cluster memberships and plots the spots with
#' colors indicating the cluster they belong to. Optionally, the user can annotate
#' tumor/stroma compartments if ESTIMATE scores are available.
#'
#' @param x an STList with hierarchical cluster memberships.
#' @param plot_who an integer indicating the spatial array(s) to be plotted.
#' @param purity logical, whether or not annotate tumor/stroma.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param visium whether or not to reverse axes for Visium slides.
#' @return a list with the requested plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # cluster_p <- plot_STclusters(melanoma, purity=T, plot_who=c(2,3), visium=F)
#'
#' @export
#
#
plot_STclusters <- function(x, plot_who=NULL, purity=F, color_pal='light', visium=T){

  require('ggplot2')
  require('magrittr')

  if(is.null(plot_who)){
    subjs <- (1:length(x@st_clusters$clust_dfs))
  } else{
    subjs <- plot_who
  }

  if(rlang::is_empty(x@st_clusters)){
    stop('No clusters found in STList.')
  }

  plot_list <- list()

  for(i in subjs){

    subj_listid <- grep(paste0("sub", i, "_"), names(x@st_clusters$clust_dfs), value=T)

    for(s in subj_listid){

      df <- x@st_clusters$clust_dfs[[s]]

      if(purity){
        # Test ESTIMATE clusters are available.
        if(is.null(x@cell_deconv$ESTIMATE)){
          stop('No ESTIMATE cluster annotations available')
        }
        purity_df <- x@cell_deconv$ESTIMATE[[i]]$purity_clusters
        df <- dplyr::full_join(df, purity_df, by='libname')
        df$cluster <- as.factor(df$cluster)
        names(df)[5] <- "EstCluster"
      }

      nas <- complete.cases(df$WCluster)
      df_nonas <- df[nas, ]
      df$WCluster <- as.factor(tidyr::replace_na(as.vector(df$WCluster), "no_cluster"))

      sweight <- as.vector(stringr::str_match(s, paste0("spw0\\.?[0-9]*"))) %>% gsub('spw', '', .)

      if(x@st_clusters$type == 'dtc'){
        title_p <- paste0("Cluster assignments (dynamicTreeCut)\n", "weight=", sweight, " - subj ", i)
      } else{
        kval <- as.vector(stringr::str_match(s, paste0("k[0-9]+"))) %>% gsub('k', '', .)
        title_p <- paste0("Cluster assignments k=", kval, "\n", "weight=", sweight, " - subj ", i)
      }

      # Determnine number of clusters generated, and cerate color palette.
      nclust <- max(as.numeric(levels(df_nonas$WCluster)))

      # Create color palette.
      clust_cols = color_parse(color_pal, n_cats=nclust)

      # Determine color scheme
      # if(length(color_pal) == 1){
      #   if(color_pal %in% khroma::info()[[1]]){
      #     # Creates color palette function.
      #     p_palette <- khroma::colour(color_pal)
      #     clust_cols <- as.vector(p_palette(nclust))
      #   }else{
      #     stop('Provide a Khroma palette name or a vector with 2 or more colors.')
      #   }
      # }else{
      #   clust_cols <- color_pal
      # }

      if(nclust > length(clust_cols)){
        stop('Provide enough colors to plot clusters.')
      }

      p <- ggplot()
      if(ncol(df) == 4){
        p <- p + geom_point(data=df, aes(x=xpos, y=ypos, color=WCluster), size=0.7, shape=19)
        p <- p + geom_point(data=df_nonas, aes(x=xpos, y=ypos, color=WCluster), size=0.5)
        p <- p + scale_color_manual(values=c(clust_cols, 'gray50'))
        p <- p + labs(title=title_p, color='Clusters')
        p <- p + guides(color=guide_legend(override.aes=list(size=2)))
      } else{
        p <- p + geom_point(data=df, aes(x=xpos, y=ypos, color=WCluster), size=0.5, shape=19)
        p <- p + geom_point(data=df_nonas, aes(x=xpos, y=ypos, color=WCluster, shape=EstCluster), size=1.2)
        p <- p + scale_color_manual(values=c(clust_cols, 'gray50'))
        p <- p + scale_shape_manual(values=c(3, 15))
        p <- p + labs(title=title_p, color='Clusters', shape='tumor/stroma')
        p <- p + guides(shape=guide_legend(override.aes=list(size=2)), color=guide_legend(override.aes=list(size=2)))
      }

      p <- p + ylab('Y Position') + xlab('X Position')

      if(visium){
        p <- p + scale_y_reverse() + scale_x_reverse() + coord_fixed(ratio=1.7)
      } else{
        p <- p + coord_fixed(ratio=1)
      }

      p <- p + theme_classic()
      p <- p + theme(plot.title=element_text(size=10), legend.text=element_text(size=10))

      plot_list[[s]] <- p

    }

  }

  return(plot_list)

}
