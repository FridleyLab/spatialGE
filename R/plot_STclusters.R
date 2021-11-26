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
#' @param plot_who an integer indicating the spatial array(s) to be plotted. Numbers
#' follow the order of `names(x@counts)`. If NULL, will plot all spatial arrays.
#' @param ks the k values to plot
#' @param weights the spatial weights to plot
#' @param purity logical, whether or not annotate tumor/stroma.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param visium whether or not to reverse axes for Visium slides.
#' @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
#' @return a list with the requested plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # cluster_p <- plot_STclusters(melanoma, purity=T, plot_who=c(2,3), visium=F)
#'
#' @export
#
#
plot_STclusters <- function(x, plot_who=NULL, ks=NULL, weights=NULL, purity=F, color_pal='light', visium=T, ptsize=NULL){

  require('ggplot2')
  require('magrittr')

  if(rlang::is_empty(x@st_clusters)){
    stop('No clusters found in STList.')
  }

  # Define which samples to plot
  if(is.null(plot_who)){
    subjs <- names(x@st_clusters)
  } else{
    subjs <- grep(paste0(names(x@tr_counts[plot_who]), '_', collapse='|'), names(x@st_clusters), value=T)
  }
  # Define which weights to plot
  if(!is.null(weights)){
    if(x@misc[['STclust_cuttype']] == 'multiK'){
      subjs = grep(paste0('spw', weights, '_', collapse = '|'), subjs, value=T)
    } else{
      subjs = grep(paste0('spw', weights, '$', collapse = '|'), subjs, value=T)
    }
  }
  # Define which k values to plot
  if(!is.null(ks) && x@misc[['STclust_cuttype']] == 'multiK'){
    subjs = grep(paste0('_k', ks, '$', collapse = '|'), subjs, value=T)
  }

  # Throw error if no cluster configurations were found
  if(length(subjs) == 0){
    stop('No cluster results found. Please check requested samples, k, or spatial weights.')
  }

  plot_list <- list()
  for(s in subjs){

    df = x@st_clusters[[s]]

    if(purity){
      # Test ESTIMATE clusters are available.
      if(is.null(x@cell_deconv[['ESTIMATE']])){
        stop('No tumor/stroma annotations available.')
      }
      sampleid = gsub("_spw[0-9\\.k_]*$", "", s)
      purity_df <- x@cell_deconv[['ESTIMATE']][[sampleid]]$purity_clusters
      df <- dplyr::full_join(df, purity_df, by='libname')
      df$cluster <- as.factor(df$cluster)
      names(df)[5] <- "EstCluster"
    }

    # Replace NAs to string
    nas <- complete.cases(df$WCluster)
    df_nonas <- df[nas, ]
    df$WCluster <- as.factor(tidyr::replace_na(as.vector(df$WCluster), "no_cluster"))

    # Prepare titles for plots
    sweight <- as.vector(stringr::str_match(s, paste0("spw0\\.?[0-9]*"))) %>% gsub('spw', '', .)
    if(x@misc[['STclust_cuttype']] == 'dtc'){
      title_p <- paste0("Cluster assignments (dynamicTreeCut)\nweight:", sweight, "\nsample: ", gsub("_spw[\\.0-9]+", '', s))
    } else{
      kval <- as.vector(stringr::str_match(s, paste0("k[0-9]+"))) %>% gsub('k', '', .)
      title_p <- paste0("Cluster assignments k=", kval, "\nweight=", sweight, "\nsample ", gsub("_spw[\\.0-9]_k[0-9]+", '', s))
    }

    # Determnine number of clusters generated, and cerate color palette.
    nclust <- max(as.numeric(levels(df_nonas$WCluster)))
    # Create color palette.
    clust_cols = color_parse(color_pal, n_cats=nclust)

    # Define size of points
    if(is.null(ptsize)){
      ptsize = 0.5
    }

    p <- ggplot()
    if(ncol(df) == 4){
      p <- p + geom_point(data=df, aes(x=xpos, y=ypos, color=WCluster), size=ptsize, shape=19)
      p <- p + geom_point(data=df_nonas, aes(x=xpos, y=ypos, color=WCluster), size=0.5)
      p <- p + scale_color_manual(values=c(clust_cols, 'gray50'))
      p <- p + labs(title=title_p, color='Clusters')
      p <- p + guides(color=guide_legend(override.aes=list(size=2)))
    } else{
      #p <- p + geom_point(data=df, aes(x=xpos, y=ypos, color=WCluster), size=(ptsize-0.2), shape=19)
      p <- p + geom_point(data=df_nonas, aes(x=xpos, y=ypos, color=WCluster, shape=EstCluster), size=ptsize)
      p <- p + scale_color_manual(values=c(clust_cols, 'gray50'))
      p <- p + scale_shape_manual(values=c(3, 15))
      p <- p + labs(title=title_p, color='Clusters', shape='tumor/stroma')
      p <- p + guides(shape=guide_legend(override.aes=list(size=2)), color=guide_legend(override.aes=list(size=2)))
    }
    p <- p + ylab('Y Position') + xlab('X Position')
    if(visium){
      p <- p + scale_y_reverse() + #scale_x_reverse() +
        coord_fixed(ratio=1.7)
    } else{
      p <- p + coord_fixed(ratio=1)
    }
    p <- p + theme_classic()
    p <- p + theme(legend.text=element_text(size=10))

    plot_list[[s]] <- p

  }

  return(plot_list)
}

