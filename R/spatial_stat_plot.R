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
#' @param samplevar a string indicating the name of the variable in the clinical
#' data frame. If NULL, uses sample
#' @param gene the name of the gene to plot.
#' @param cell the name of the cell type to plot.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @return an STList containing the plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # stat_p <- spatial_stat_plot(melanoma, samplevar='survival_months', gene='CD74')
#'
#' @export
#
#
spatial_stat_plot <- function(x=NULL, samplevar=NULL, gene=NULL, cell=NULL, color_pal="muted") {

  require('magrittr')
  require('ggplot2')

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

  # Extract clinical data from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  if(!is.null(samplevar)){
    clinvar_vals <- x@clinical[[samplevar]]
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name=samplevar)
    plot_labs <- as.character(x@clinical[[1]])
    clinvar_vals <- clinvar_vals %>%
      tibble::add_column(., samples=plot_labs)
  }else{
    samplevar <- 'sample_id'
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
    x <- gene_moran_I(x, genes=gene, who=c(1:length(x@tr_counts)))
    x <- gene_geary_C(x, genes=gene, who=c(1:length(x@tr_counts)))
    x <- gene_getis_Gi(x, genes=gene, who=c(1:length(x@tr_counts)))

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
  }

  if(!is.null(gene) && is.null(cell)){
    x <- gene_moran_I(x, genes=gene, who=c(1:length(x@tr_counts)))
    x <- gene_geary_C(x, genes=gene, who=c(1:length(x@tr_counts)))
    x <- gene_getis_Gi(x, genes=gene, who=c(1:length(x@tr_counts)))

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
  }

  # Get number of categories from selected
  n_cats <- nlevels(as.factor(clinvar_vals[[samplevar]]))
  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)

  # Assocuate colors with categories.
  names(cat_cols) <- levels(as.factor(clinvar_vals[[samplevar]]))

  # Define name of cell/gene being plotted for title.
  if(!is.null(gene)){
    gene_or_cell = gene
  }

  if(!is.null(cell)){
    gene_or_cell = cell
  }

  # NOTE: INSTEAD OF NUMBERS, WOULD BE GREAT TO HAVE SAMPLE ID PLOTTED
  moran_p <-   ggplot(clinvar_vals) +
    geom_point(aes(x=moran, y=get(samplevar), color=samples), size=3) +
    ggtitle(paste0('Moran\'s I and ', samplevar, '\n', gene_or_cell)) +
    xlab('Moran\'s I') +
    ylab(samplevar) +
    theme_light()

  geary_p <-   ggplot(clinvar_vals) +
    geom_point(aes(x=geary, y=get(samplevar), color=samples), size=3) +
    ggtitle(paste0('Geary\'s C and ', samplevar, '\n', gene_or_cell)) +
    xlab('Geary\'s C') +
    ylab(samplevar) +
    theme_light()

  getis_p <-   ggplot(clinvar_vals) +
    geom_point(aes(x=getis, y=get(samplevar), color=samples), size=3) +
    ggtitle(paste0('Getis-Ord Gi and ', samplevar, '\n', gene_or_cell)) +
    xlab('Getis-ord Gi') +
    ylab(samplevar) +
    theme_light()

  plist <- list(moran_p, geary_p, getis_p)

  # Store plot in STList
  if(!is.null(gene)){
    x@pheno_plots[[gene]] = plist
  } else if(!is.null(cell)){
    x@pheno_plots[[cell]] = plist
  }

  # Print plot.
  multiplist = ggpubr::ggarrange(plotlist=plist, ncol=3, common.legend=T, legend='bottom')
  print(multiplist)

  return(x)
}
