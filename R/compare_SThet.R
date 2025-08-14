##
#' @title compare_SThet: Compares spatial autocorrelation statistics across samples
#' @description Plots the spatial autocorrelation statistics of genes across samples and
#' colors samples acording to sample metadata.
#' @details
#' This function takes the names of genes and their Moran's I or Geary's C computed for
#' multiple samples and to provide a muti-sample comparison. Samples in the plot can
#' be colored according to sample metadata to explore potential associations between
#' spatial distribution of gene expression and sample-level data.
#'
#' @param x an STlist.
#' @param samplemeta a string indicating the name of the variable in the clinical
#' data frame. If NULL, uses sample names
#' @param genes the name(s) of the gene(s) to plot.
#' @param color_by the variable in `x@spatial_meta` used to color points in the plot.
#' If NULL, each sample is assigned a different color
#' @param categorical logical indicating whether or not to treat `color_by` as a
#' categorical variable. Default is TRUE
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param ptsize a number specifying the size of the points. Passed to the `size`
#' aesthetic.
#' @return a list of plots
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
#
compare_SThet = function(x=NULL, samplemeta=NULL, genes=NULL,  color_by=NULL,
                         categorical=TRUE, color_pal="muted", ptsize=1) {

  # To prevent NOTES in R CMD check
  . = NULL

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STList.")
  }

  if(is.null(genes)){
    stop('Please, enter one gene or more to generate plot.')
  }

  # Color by patient if no color_by
  if(is.null(color_by)){
    color_by = colnames(x@sample_meta)[1]
  }
  # If no metadata variable has been input, use sample names
  if(is.null(samplemeta)){
    samplemeta = colnames(x@sample_meta)[1]
  }

  # Extract clinical data from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  meta_df = x@sample_meta %>%
    dplyr::select(1, !!!samplemeta, !!!color_by)
  if(categorical){ # To avoid treating numbers as categories
    meta_df = meta_df %>%
      dplyr::mutate(!!color_by := as.factor(.[[!!color_by]]))
  }
  meta_df[['moran']] = NA
  meta_df[['geary']] = NA

  metadf_ls = list()
  for(gene in genes){
    metadf_ls[[gene]] = meta_df
    for(i in meta_df[[1]]){
      if(gene %in% x@gene_meta[[i]][['gene']]){
        metadf_ls[[gene]][['moran']][ metadf_ls[[gene]][[1]] == i] = x@gene_meta[[i]][['moran_i']][ x@gene_meta[[i]][['gene']] == gene ]
        metadf_ls[[gene]][['geary']][ metadf_ls[[gene]][[1]] == i] = x@gene_meta[[i]][['geary_c']][ x@gene_meta[[i]][['gene']] == gene ]
      } else{
        metadf_ls[[gene]][['moran']][ metadf_ls[[gene]][[1]] == i] = NA
        metadf_ls[[gene]][['geary']][ metadf_ls[[gene]][[1]] == i] = NA
      }
    }
    metadf_ls[[gene]][['gene']] = gene
  }
  metadf_ls = do.call(dplyr::bind_rows, metadf_ls) #%>%
    #dplyr::rename(samplename=1)

  rm(meta_df) # Clean env

  res_p = spatial_stat_plot_gene(meta_df=metadf_ls, samplemeta=samplemeta, color_by=color_by, color_pal=color_pal, ptsize=ptsize)

  return(res_p)
}

# Helpers ----------------------------------------------------------------------

##
# @title compareSThet_plot
# @param meta_df a data frame with samples and spatial stats to plot
# @param color_pal the name of a color palette
#
spatial_stat_plot_gene = function(meta_df=NULL, samplemeta=NULL, color_by=NULL, color_pal=NULL, ptsize=1){
  if(is.null(samplemeta)){
    samplemeta = names(meta_df)[1]
  }


  # Get number of categories from selected
  n_cats = nlevels(as.factor(meta_df[[color_by]]))
  # Create color palette
  cat_cols = color_parse(color_pal, n_cats=n_cats)
  # Associate colors with categories.
  names(cat_cols) = levels(as.factor(meta_df[[color_by]]))

  res_plots = list()

  if(any(!is.na(meta_df[['moran']]))){
    res_plots[['moran']] = ggplot(meta_df) +
      geom_point(aes(x=moran, y=.data[[samplemeta]], color=.data[[color_by]]), size=ptsize) +
      #ggrepel::geom_text_repel(aes(x=moran, y=.data[[samplemeta]], label=.data[[color_by]])) +
      ggtitle(paste0('Moran\'s I and ', samplemeta)) +
      xlab('Moran\'s I') +
      ylab(samplemeta)  +
      facet_wrap(~gene)

    if(is.numeric(meta_df[[color_by]])){
      res_plots[['moran']] = res_plots[['moran']] +
        scale_color_gradientn(colors=as.vector(cat_cols), # SHOULD PRBABLY CHANGE COLOR_PARSE FUNCTION TO OUPUT MIN/MID/MAX COLORS WHEN CONTINUOUS
                              guide=guide_legend(label.theme=element_text(angle=0),
                                                 override.aes=list(size=2)))
    } else{
      res_plots[['moran']] = res_plots[['moran']] +
        scale_color_manual(values=cat_cols,
                           guide=guide_legend(label.theme=element_text(angle=0),
                                              override.aes=list(size=2)))
    }

    res_plots[['moran']] = res_plots[['moran']] +
      theme_light() +
      theme(#legend.title=element_blank(),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  }

  if(any(!is.na(meta_df[['geary']]))){
    res_plots[['geary']] = ggplot(meta_df) +
      geom_point(aes(x=geary, y=.data[[samplemeta]], color=.data[[color_by]]), size=ptsize) +
      #ggrepel::geom_text_repel(aes(x=moran, y=.data[[samplemeta]], label=.data[[color_by]])) +
      ggtitle(paste0('Geary\'s C and ', samplemeta)) +
      xlab('Geary\'s C') +
      ylab(samplemeta)  +
      facet_wrap(~gene)

    if(is.numeric(meta_df[[color_by]])){
      res_plots[['geary']] = res_plots[['geary']] +
        scale_color_gradientn(colors=as.vector(cat_cols), # SHOULD PRBABLY CHANGE COLOR_PARSE FUNCTION TO OUPUT MIN/MID/MAX COLORS WHEN CONTINUOUS
                              guide=guide_legend(label.theme=element_text(angle=0),
                                                 override.aes=list(size=2)))
    } else{
      res_plots[['geary']] = res_plots[['geary']] +
        scale_color_manual(values=cat_cols,
                           guide=guide_legend(label.theme=element_text(angle=0),
                                              override.aes=list(size=2)))
    }

    res_plots[['geary']] = res_plots[['geary']] +
      theme_light() +
      theme(#legend.title=element_blank(),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  }

  # Print plot.
  res_plots = ggpubr::ggarrange(plotlist=res_plots, common.legend=T, legend='bottom', ncol=1)

  return(res_plots)
}


