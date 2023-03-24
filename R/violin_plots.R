##
#' @title violin_plots: Generates violin plots and boxplots of meta data or gene expression
#' @description Generates violin plots and boxplots of continuous variables in the spatial meta
#' data or expression of genes to make comparisons across samples
#' @details The function allows to visualize the distribution of spot/cell total counts, genes, or
#' expression of specific genes across all samples for comparative purposes.
#'
#' @param x an STlist
#' @param plot_meta vector of variables in `x@spatial_meta` to plot distributions. Default is 'total_counts'
#' @param genes vector of genes to plot expression distribution
#' @param samples samples to include in the plot. Default (NULL) includes all samples
#' @param data_type one of 'tr' or 'raw', to plot transformed or raw counts
#' @param color_pal  a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with colors
#' @param plot_type one of "violin" or "box" (violin plots or box plots)
#' @param ptsize the size of points in the plots
#' @param ptalpha the transparency of points in plots
#' @param ptjitter the dispersion of points in the plots
#'
#' @export
#'
#' @import ggplot2
#
#
violin_plots = function(x=NULL, plot_meta=NULL, genes=NULL, samples=NULL, data_type='tr',
                        color_pal='roma', plot_type='violin', ptsize=0.5, ptalpha=0.5){
  #require('ggplot2')

  # Define samples to plot if NULL or numeric
  if(is.null(samples)){
    samples = names(x@counts)
  } else if(is.numeric(samples)){
    samples = names(x@counts)[samples]
  }

  # Define if plotting meta data or genes
  if(!is.null(plot_meta) & !is.null(genes)){
    stop("Plots can be generated only for plot_meta OR genes at a time.")
  } else if(!is.null(plot_meta)){
    plist = violin_plots_spot(x=x, plot_meta=plot_meta, samples=samples,
                               plot_type=plot_type, color_pal=color_pal, ptsize=ptsize, ptalpha=ptalpha)
  } else {
    plist = violin_plots_gene(x=x, genes=genes, samples=samples, data_type=data_type,
                               plot_type=plot_type, color_pal=color_pal, ptsize=ptsize, ptalpha=ptalpha)
  }
  return(plist)
}


# Helpers ----------------------------------------------------------------------

##
# violin_plots_spot
#

violin_plots_spot = function(x=NULL, plot_meta='total_counts', samples=NULL,
                              plot_type='violin', color_pal='roma', ptsize=0.5, ptalpha=0.5){

  p_list = list()
  # Loop through meta data variables
  for(meta in plot_meta){
    # Subset samples if meta data not in all samples
    samples_tmp = samples
    df_tmp = tibble::tibble()
    for(i in samples){
      if(!(meta %in% colnames(x@spatial_meta[[i]]))){
        cat(paste0('Sample ', i, ' does not contain ', meta, '.\n'))
        samples_tmp = grep(paste0('^', i, '$'), samples_tmp, value=T, invert=T)
      } else{
        df_tmp = dplyr::bind_rows(df_tmp,
                                  x@spatial_meta[[i]] %>%
                                    dplyr::select(libname, meta) %>%
                                    tibble::add_column(samplename=i))
      }
    }

    # Define color palette
    meta_cols = color_parse(color_pal, n_cats=length(samples))

    # Create plot
    p_list[[meta]] = ggplot2::ggplot(df_tmp, aes(x=samplename, y=.data[[meta]], color=samplename))
    if(plot_type == 'box'){
      p_list[[meta]] = p_list[[meta]] +
        ggplot2::geom_boxplot(outlier.size=ptsize)
    } else if(plot_type == 'violin'){
      p_list[[meta]] = p_list[[meta]] +
        ggplot2::geom_violin() +
        ggforce::geom_sina(size=ptsize)
    }
    p_list[[meta]] = p_list[[meta]] +
      ggplot2::xlab(NULL) +
      ggplot2::scale_color_manual(values=meta_cols) +
      ggplot2::theme(axis.text.x=element_text(angle=30, hjust=1),
                     panel.border=element_rect(fill=NA, color='black'),
                     legend.title=element_blank())
  }

  return(p_list)
}


##
# violin_plots_gene
#

violin_plots_gene = function(x=NULL, genes=NULL, samples=NULL, data_type='tr',
                              plot_type='violin', color_pal='roma', ptsize=0.5, ptalpha=0.5){
  # Extract data slot
  if(data_type == 'tr'){
    expr_tmp = x@tr_counts
    p_title = 'normalized expr - '
  } else {
    expr_tmp = x@counts
    p_title = 'raw expr - '
  }

  p_list = list()
  # Loop through list of genes
  for(gene in genes){
    # Subset samples if gene not in all samples
    samples_tmp = samples
    df_tmp = tibble::tibble()
    for(i in samples){
      if(!(gene %in% rownames(expr_tmp[[i]]))){
        cat(paste0('Sample ', i, ' does not contain ', gene, '.\n'))
        samples_tmp = grep(paste0('^', i, '$'), samples_tmp, value=T, invert=T)
      } else{
        df_tmp = dplyr::bind_rows(df_tmp,
                                  tibble::tibble(libname=colnames(expr_tmp[[i]]),
                                                 geneexpr=expr_tmp[[i]][gene, ],
                                                 samplename=i))
      }
    }

    # Define color palette
    gene_cols = color_parse(color_pal, n_cats=length(samples))

    # Create plot
    p_list[[gene]] = ggplot2::ggplot(df_tmp, aes(x=samplename, y=.data[['geneexpr']], color=samplename))
    if(plot_type == 'box'){
      p_list[[gene]] = p_list[[gene]] +
        ggplot2::geom_boxplot(outlier.size=ptsize)
    } else if(plot_type == 'violin'){
      p_list[[gene]] = p_list[[gene]] +
        ggplot2::geom_violin() +
        ggforce::geom_sina(size=ptsize)
    }
    p_list[[gene]] = p_list[[gene]] +
      ggplot2::ylab(paste0(p_title, gene)) +
      ggplot2::xlab(NULL) +
      ggplot2::scale_color_manual(values=gene_cols) +
      ggplot2::theme(axis.text.x=element_text(angle=30, hjust=1),
                     panel.border=element_rect(fill=NA, color='black'),
                     legend.title=element_blank())
  }

  return(p_list)
}


