##
#' @title spatial_stat_plot: Plots association of samples variable and spatial statistics
#' @description Plots the value of a sample/clinical variable against the spatial
#' autocorrelation statistics.
#' @details
#' This function takes the name of a gene or deconvoluted cell score and calculates
#' its spatial autocorrelation statistics. Then generates a plot showing the variation
#' of a sample-associated or clinical variable with respect to the Moran's I, Geary's C,
#' and Getis-Ord Gi values of a gene.
#'
#' @param x an STList.
#' @param var a string indicating the name of the variable in the clinical
#' data frame. If NULL, uses sample
#' @param gene the name of the gene to plot.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @return an STList containing the plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # stat_p <- spatial_stat_plot(melanoma, var='survival_months', gene='CD74')
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
#
spatial_stat_plot <- function(x=NULL, var=NULL, gene=NULL, color_pal="muted") {

  #require('magrittr')
  #require('ggplot2')

  # Spatial plots with cell scores might get removed from package
  cell=NULL

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  if(!is.null(gene) && length(gene) != 1){
    stop('At the moment, only one gene at a time is enabled.')
  }

  if(!is.null(cell) && length(cell) != 1){
    stop('At the moment, only one cell type at a time is enabled.')
  }

  if(!is.null(gene) && !is.null(cell)){
    stop('At the moment, only one gene OR one cell type at a time is enabled.')
  }

  if(is.null(gene) && is.null(cell)){
    stop('Please, enter a gene to generate plot.')
  }

  # Extract clinical data from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  if(!is.null(var)){
    clinvar_vals <- x@clinical[[var]]
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name=var)
    plot_labs <- as.character(x@clinical[[1]])
    clinvar_vals <- clinvar_vals %>%
      tibble::add_column(., samples=plot_labs)
  }else{
    var <- 'sample_id'
    clinvar_vals <- 1:length(x@counts)
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='sample_id')
    plot_labs <- as.character(1:length(x@tr_counts))
    clinvar_vals <- clinvar_vals %>%
      tibble::add_column(., samples=plot_labs)
  }

  clinvar_vals$moran <- rep(NA, length(x@tr_counts))
  clinvar_vals$geary <- rep(NA, length(x@tr_counts))
  clinvar_vals$getis <- rep(NA, length(x@tr_counts))

  if(!is.null(gene) && is.null(cell)){

    x = gene_moran_I(x, genes=gene, who=c(1:length(x@tr_counts)))
    x = gene_geary_C(x, genes=gene, who=c(1:length(x@tr_counts)))
    x = gene_getis_Gi(x, genes=gene, who=c(1:length(x@tr_counts)))

    clinvar_vals = get_gene_stats(x=x, gene=gene, clinvar_vals=clinvar_vals)

    x = spatial_stat_plot_gene(x=x, clinvar_vals=clinvar_vals, var=var, feature=gene, color_pal=color_pal)

  } else if(is.null(gene) && !is.null(cell)){
    get_cell_stats = NULL
    spatial_stat_cell = NULL
  }

  return(x)
}

# Helpers ----------------------------------------------------------------------

##
#' @title get_gene_stats
#' @param x an STList
#' @param gene a gene name
#' @param clinvar_vals a data frame to populate with spatial statistics
#

get_gene_stats = function(x=NULL, gene=NULL, clinvar_vals=NULL){

  for(i in 1:length(x@tr_counts)){
    moran_estimate = as.vector(x@gene_het[[gene]][[i]]$morans_I$estimate[[1]])
    geary_estimate = as.vector(x@gene_het[[gene]][[i]]$gearys_C$estimate[[1]])
    getis_estimate = as.vector(x@gene_het[[gene]][[i]]$getis_ord_Gi$estimate[[1]])

    if(length(moran_estimate) != 0){
      clinvar_vals$moran[i] = moran_estimate
    } else{
      clinvar_vals$moran[i] = NA
    }
    if(length(geary_estimate) != 0){
      clinvar_vals$geary[i] = geary_estimate
    } else{
      clinvar_vals$geary[i] = NA
    }
    if(length(getis_estimate) != 0){
      clinvar_vals$getis[i] = getis_estimate
    } else{
      clinvar_vals$getis[i] = NA
    }
  }
  return(clinvar_vals)
}

##
#' @title spatial_stat_plot_gene
#' @param x an STList
#' @param clinvar_vals a data frame with data frame to extract data for plot
#' @param var variable to plot from clinical data
#' @param feature the gene to plot statistics for
#' @param color_pal the name of a color palette
#

spatial_stat_plot_gene = function(x=x, clinvar_vals=NULL, var=NULL, feature=NULL, color_pal=NULL){

  # Get number of categories from selected
  n_cats = nlevels(as.factor(clinvar_vals[[var]]))
  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)

  # Assocuate colors with categories.
  names(cat_cols) <- levels(as.factor(clinvar_vals[[var]]))

  # NOTE: INSTEAD OF NUMBERS, WOULD BE GREAT TO HAVE SAMPLE ID PLOTTED
  moran_p = ggplot(clinvar_vals) +
    geom_point(aes(x=moran, y=get(var), color=samples), size=3) +
    ggtitle(paste0('Moran\'s I and ', var, '\n', feature)) +
    xlab('Moran\'s I') +
    ylab(var) +
    theme_light()

  geary_p = ggplot(clinvar_vals) +
    geom_point(aes(x=geary, y=get(var), color=samples), size=3) +
    ggtitle(paste0('Geary\'s C and ', var, '\n', feature)) +
    xlab('Geary\'s C') +
    ylab(var) +
    theme_light()

  getis_p = ggplot(clinvar_vals) +
    geom_point(aes(x=getis, y=get(var), color=samples), size=3) +
    ggtitle(paste0('Getis-Ord Gi and ', var, '\n', feature)) +
    xlab('Getis-ord Gi') +
    ylab(var) +
    theme_light()

  x@spstats_plots[[feature]] = list(moran_p, geary_p, getis_p)

  # Print plot.
  multiplist = ggpubr::ggarrange(plotlist=x@spstats_plots[[feature]], ncol=3, common.legend=T, legend='bottom')
  print(multiplist)

  return(x)
}


