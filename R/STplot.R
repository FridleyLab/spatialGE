##
#' @title STplot: Plots of gene expression, cluster memberships, and metadata in spatial context
#' @description Generates a plot of the location of spots/cells within an spatial
#' sample, and colors them according to gene expression levels or spot/cell-level
#' metadata
#' @details
#' The function takes an STlist and plots the cells or spots in their spatial context.
#' The users can color the spots/cells according to the expression of selected genes,
#' cluster memberships, or any spot/cell level metadata included in `x@spatial_meta`.
#' The function also can average expression of gene sets.
#'
#' @param x an STlist
#' @param samples a vector of numbers indicating the ST samples to plot, or their
#' sample names. If vector of numbers, it follow the order of samples in `names(x@counts)`.
#' If NULL, the function plots all samples
#' @param genes a vector of gene names or a named list of gene sets. In the latter
#' case, the averaged expression of genes within the sets is plotted
#' @param plot_meta a column name in `x@spatial_meta` to plot
#' @param data_type one of 'tr' or 'raw', to plot transformed or raw counts
#' respectively
#' @param ks the k values to plot or 'dtc' to plot results from `dynamicTreeCut`
#' clustering solutions. Requires previous analysis with `STclust`
#' @param ws the spatial weights to plot samples if `STclust` was used
#' @param deepSplit a logical or positive number indicating the `deepSplit`, if
#' samples were analyzed with `STclust`
#' @param color_pal a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with enough color names or HEX values
#' @param ptsize a number specifying the size of the points. Passed to the `size`
#' @param txsize a number controlling the size of the text in the plot title and legend title. Passed to the `element_text`
#' aesthetic.
#' @return a list of plots
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- transform_data(melanoma)
#'   STplot(melanoma, gene='MLANA', samples='ST_mel1_rep2', ptsize=1)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#
STplot = function(x, samples=NULL, genes=NULL, plot_meta=NULL, ks='dtc', ws=NULL, deepSplit=NULL, color_pal=NULL, data_type='tr', ptsize=NULL, txsize=NULL){

  # Check if data set is Visium to flip y axis (may consider remove this in the future)
  visium = F
  if(x@misc[['platform']] == 'visium'){
    visium = T
  }

  # Define if expression or metadata is to be plotted
  if(!is.null(genes)){
    # Set default color if NULL input
    if(is.null(color_pal)){
      color_pal = 'BuRd'
    }
    if(is.list(genes)){
      plot_list = plot_spatial_geneset(x=x, genes=genes, samples=samples, color_pal=color_pal, visium=visium, ptsize=ptsize)
    } else{
      plot_list = plot_spatial_expression(x=x, genes=genes, samples=samples, color_pal=color_pal, data_type=data_type, visium=visium, ptsize=ptsize)
    }
  } else{
    plot_list = plot_spatial_meta(x=x, samples=samples, ks=ks, ws=ws, deepSplit=deepSplit, plot_meta=plot_meta, color_pal=color_pal, visium=visium, ptsize=ptsize, txsize=txsize)
  }

  return(plot_list)
}

