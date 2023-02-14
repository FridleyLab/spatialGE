##
#' @title STplot: Plots of gene expression, cluster memberships, and metadata in spatial context
#' @description Generates a plot of the location of cells/spots within an spatial sample,
#' and colors them according to gene expression levels, STclust memberships, or spatial
#' metadata
#' @details
#' The function takes an STList and plots the cells or spots in the spatial context of
#' the tissue sample. The users can color the cells/spots according to the expression of
#' selected genes, cluster memberships inferred by STclust, or any metadata included in
#' `x@spatial_meta`.
#'
#' @param x an STList
#' @param samples a vector of numbers indicating the ST samples to plot, or their
#' sample names. If vector of numbers, it follow the order of `names(x@counts)`.
#' If NULL, the function plots all samples
#' @param genes a vector of gene names or a named list of gene sets. In the latter
#' case, the averaged expression of genes within the sets is plotted
#' @param data_type one of 'tr' or 'raw', to plot transformed or raw counts
#' respectively
#' @param ks the k values to plot or 'dtc' to plot results from `dynamicTreeCut`
#' clustering solutions
#' @param ws the spatial weights used during STclust
#' @param deepSplit a logical or integer indicating the `deepSplit` parameter
#' used in STclust
#' @param plot_meta a column name in `x@spatial_meta` to plot
#' @param color_pal a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with colors with enough elements to plot categories
#' @param visium whether or not to reverse axes for Visium slides
#' @param ptsize a number specifying the size of the points. Passed to `size` aesthetic.
#' @return a list of plots
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
STplot = function(x, samples=NULL, genes=NULL, data_type='tr', ks='dtc', ws=NULL, deepSplit=NULL, plot_meta=NULL, color_pal=NULL, visium=T, ptsize=NULL, image=F){

  # Define if expression or metadata is to be plotted
  if(!is.null(genes)){
    # Set default color if NULL input
    if(is.null(color_pal)){
      color_pal = 'BuRd'
    }
    if(is.list(genes)){
      plot_list = plot_spatial_geneset(x=x, genes=genes, samples=samples, color_pal=color_pal, visium=visium, ptsize=ptsize)
    } else{
      plot_list = plot_spatial_expression(x=x, genes=genes, samples=samples, color_pal=color_pal, data_type=data_type, image=image, visium=visium, ptsize=ptsize)
    }
  } else{
    # Set default color if NULL input
    if(is.null(color_pal)){
      color_pal = 'light'
    }
    plot_list = plot_spatial_meta(x=x, samples=samples, ks=ks, ws=ws, deepSplit=deepSplit, plot_meta=plot_meta, color_pal=color_pal, visium=visium, ptsize=ptsize)
  }

  return(plot_list)
}

