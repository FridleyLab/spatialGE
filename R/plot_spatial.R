##
#' @title plot_spatial: Plot cluster memberships of ST spots
#' @description Generates a plot of the location of spots within an spatial array,
#' and colors them according to spatially-weighted hierarchical clustering assignments.
#' @details
#' The function takes an STList with cluster memberships and plots the spots with
#' colors indicating the cluster they belong to. Optionally, the user can annotate
#' tumor/stroma compartments if ESTIMATE scores are available.
#'
#' @param x an STList with hierarchical cluster memberships.
#' @param samples an integer indicating the spatial array(s) to be plotted. Numbers
#' follow the order of `names(x@counts)`. If NULL, will plot all spatial arrays.
#' @param ks the k values to plot
#' @param ws the spatial weights to plot
#' @param plot_meta a column name in `x@spatial_meta` to plot
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param visium whether or not to reverse axes for Visium slides.
#' @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
#' @return a list with the requested plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # cluster_p <- plot_STclusters(melanoma, samples=c(2,3), visium=F)
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
plot_spatial = function(x, samples=NULL, genes=NULL, data_type='tr', ks='dtc', ws=NULL, deepSplit=NULL, plot_meta=NULL, color_pal=NULL, visium=T, ptsize=NULL, image=F){

  # Define if expression or metadata is to be plotted
  if(!is.null(genes)){
    # Set default color if NULL input
    if(is.null(color_pal)){
      color_pal = 'BuRd'
    }
    plot_list = plot_gene_expression(x=x, genes=genes, samples=samples, color_pal=color_pal, data_type=data_type, image=image, visium=visium, ptsize=ptsize)
  } else{
    # Set default color if NULL input
    if(is.null(color_pal)){
      color_pal = 'light'
    }
    plot_list = plot_spatial_meta(x=x, samples=samples, ks=ks, ws=ws, deepSplit=deepSplit, plot_meta=plot_meta, color_pal=color_pal, visium=visium, ptsize=ptsize)
  }

  return(plot_list)
}

