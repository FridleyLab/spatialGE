##
#' @title phenovar_plot: Plots a phenotype variable against gene expression spatial autocorrelation
#' @description Plots the value of a phenotype/clinical variable for each spatial
#' array against the spatial autocorrelation statistics of a given gene.
#' @details
#' This function takes the name of a gene and calculates its spatial autocorrelation
#' statistics. Then generates a plot showing the variation of a phenotype/clinical
#' variable with respect to the Moran's I, Geary's C, and Getis-Ord Gi values of a
#' gene.
#'
#' @param x, an STList.
#' @param phenovar, a string indicating the name of the variable in the clinical
#' data frame. If NULL, uses sample
#' @param gene, the name of the gene to plot.
#' @param cell, the name of the cell type to plot.
#' @param color_pal, a string of a color palette from khroma.
#' @return plist, a list with three plots.
#' @export
#
#
plot_phenovar <- function(x=NULL, phenovar=NULL, gene=NULL, cell=NULL, color_pal="muted") {

  require('magrittr')
  require('ggplot2')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Extract clinical data from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  if(!is.null(phenovar)){
    clinvar_vals <- x@clinical[[phenovar]]
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name=phenovar)
    plot_labs <- as.character(x@clinical[[1]])
    clinvar_vals <- clinvar_vals %>%
      tibble::add_column(., samples=plot_labs)
  }else{
    phenovar <- 'array'
    clinvar_vals <- 1:length(x@counts)
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='array')
    plot_labs <- as.character(1:length(x@counts))
    clinvar_vals <- clinvar_vals %>%
      tibble::add_column(., samples=plot_labs)
  }

  clinvar_vals$moran <- rep(NA, length(x@counts))
  clinvar_vals$geary <- rep(NA, length(x@counts))
  clinvar_vals$getis <- rep(NA, length(x@counts))

  # Loop through count matrices in STList
  for(i in 1:length(x@counts)){

    if(!any(x@voom_counts[[i]][[1]] == gene)){
      cat(paste("\n", gene, "is not a gene in the normalized count matrix of array ", i, ".\n"))
      next
    } else{

      x <- gene_moran_I(x, genes=gene, subj=i)
      x <- gene_geary_C(x, genes=gene, subj=i)
      x <- gene_getis_Gi(x, genes=gene, subj=i)

      if(!is.null(gene)){
        clinvar_vals$moran[i] <- as.vector(x@gene_het[[gene]][[i]]$morans_I$estimate[[1]])
        clinvar_vals$geary[i] <- as.vector(x@gene_het[[gene]][[i]]$gearys_C$estimate[[1]])
        clinvar_vals$getis[i] <- as.vector(x@gene_het[[gene]][[i]]$getis_ord_Gi$estimate[[1]])
      } else if(!is.null(cell)){
        clinvar_vals$moran[i] <- as.vector(x@cell_het[[cell]][[i]]$morans_I$estimate[[1]])
        clinvar_vals$geary[i] <- as.vector(x@cell_het[[cell]][[i]]$gearys_C$estimate[[1]])
        clinvar_vals$getis[i] <- as.vector(x@cell_het[[gecellne]][[i]]$getis_ord_Gi$estimate[[1]])
      } else{
        stop('Please, specify a gene or cell type to plot.')
      }

    }
  }

  # Get color palette and number of colors needed.
  p_palette <- khroma::colour(color_pal)
  n_cats <- nlevels(as.factor(clinvar_vals[[phenovar]]))
  cat_cols <- as.vector(p_palette(n_cats))
  names(cat_cols) <- levels(as.factor(clinvar_vals[[phenovar]]))

# NOTE: INSTEAD OF NUMBERS, WOULD BE GREAT TO HAVE SAMPLE ID PLOTTED
moran_p <-   ggplot(clinvar_vals) +
  geom_point(aes(x=moran, y=get(phenovar), color=samples), size=3) +
  ggtitle(paste0('Moran\'s I and ', phenovar)) +
  xlab('Moran\'s I') +
  ylab(phenovar) +
  theme_light()

geary_p <-   ggplot(clinvar_vals) +
  geom_point(aes(x=geary, y=get(phenovar), color=samples), size=3) +
  ggtitle(paste0('Geary\'s C and ', phenovar)) +
  xlab('Geary\'s C') +
  ylab(phenovar) +
  theme_light()

getis_p <-   ggplot(clinvar_vals) +
  geom_point(aes(x=getis, y=get(phenovar), color=samples), size=3) +
  ggtitle(paste0('Getis-Ord Gi and ', phenovar)) +
  xlab('Getis-ord Gi') +
  ylab(phenovar) +
  theme_light()

plist <- list(moran_p, geary_p, getis_p)

return(plist)

}
