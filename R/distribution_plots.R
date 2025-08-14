##
#' @title per_unit_counts: Generates distribution plots of spot/cell meta data or gene expression
#' @description Generates violin plots, boxplots, or density plots of variables in the
#' spatial meta data or of gene expression
#' @details The function allows to visualize the distribution of spot/cell total
#' counts, total genes, or expression of specific genes across all samples for
#' comparative purposes. It also allows grouping of gene expression values by
#' categorical variables (e.g., clusters).
#'
#' @param x an STlist
#' @param plot_meta vector of variables in `x@spatial_meta` to plot distributions.
#' If 'total_counts', the function plots the counts per spot/cell. If 'total_genes',
#' the function plots the number of genes per spot/cell are plotted
#' @param genes vector of genes to plot expression distribution. If used in conjunction
#' with `plot_meta`, the expression values are grouped using that variable
#' @param samples samples to include in the plot. Default (NULL) includes all samples
#' @param data_type one of 'tr' or 'raw', to plot transformed or raw counts
#' @param color_pal  a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with colors
#' @param plot_type one of "violin", "box", or "density" (violin plots, box plots, or
#' density plots respectively). If `plot_meta` and `gene` are used together, then
#' density plots are disabled
#' @param ptsize the size of points in the plots
#' @param ptalpha the transparency of points (violin/box plot) or curves (density plots)
#' @return a list containing `ggplot2` objects
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @import ggplot2
#
#
distribution_plots = function(x=NULL, plot_meta=NULL, genes=NULL, samples=NULL, data_type='tr',
                              color_pal='okabeito', plot_type='violin', ptsize=0.5, ptalpha=0.5){

  #require('magrittr')

  # Define samples to plot if NULL or numeric
  if(is.null(samples)){
    samples = names(x@counts)
  } else if(is.numeric(samples)){
    samples = names(x@counts)[samples]
  }

  # Define if plotting meta data or genes
  if(!is.null(plot_meta) & !is.null(genes)){
    plist = cluster_gene_plots(x=x, plot_meta=plot_meta, genes=genes, samples=samples, data_type=data_type,
                               plot_type=plot_type, color_pal=color_pal, ptsize=ptsize, ptalpha=ptalpha)
  } else if(!is.null(plot_meta)){
    plist = spot_plots(x=x, plot_meta=plot_meta, samples=samples,
                       plot_type=plot_type, color_pal=color_pal, ptsize=ptsize, ptalpha=ptalpha)
  } else if(!is.null(genes)){
    plist = gene_plots(x=x, genes=genes, samples=samples, data_type=data_type,
                       plot_type=plot_type, color_pal=color_pal, ptsize=ptsize, ptalpha=ptalpha)
  } else {
    stop('Specify a spot/cell meta data variable or a gene to plot.')
  }

  return(plist)
}


# Helpers ----------------------------------------------------------------------

##
# spot_plots
#

spot_plots = function(x=NULL, plot_meta=NULL, samples=NULL,
                      plot_type='violin', color_pal='roma', ptsize=0.5, ptalpha=0.5){
  p_list = list()
  # Loop through meta data variables
  for(meta in plot_meta){
    # Subset samples if meta data not in all samples
    samples_tmp = samples
    df_tmp = tibble::tibble()
    for(i in samples){
      if(!(meta %in% colnames(x@spatial_meta[[i]]))){
        message(paste0('Sample ', i, ' does not contain ', meta, '.\n'))
        samples_tmp = grep(paste0('^', i, '$'), samples_tmp, value=T, invert=T)
      } else{
        df_tmp = dplyr::bind_rows(df_tmp,
                                  x@spatial_meta[[i]] %>%
                                    dplyr::select(libname, !!meta) %>%
                                    tibble::add_column(samplename=i))
      }
    }

    # Define color palette
    meta_cols = color_parse(color_pal, n_cats=length(samples))

    # Create plot
    if(plot_type == 'box'){
      p_list[[meta]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=samplename, y=.data[[meta]], color=samplename)) +
        ggplot2::geom_boxplot(outlier.size=ptsize)
    } else if(plot_type == 'violin'){
      p_list[[meta]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=samplename, y=.data[[meta]], color=samplename)) +
        ggplot2::geom_violin() +
        ggforce::geom_sina(size=ptsize)
    } else if(plot_type == 'density'){
      # Define bandwidth
      bandw = round(mean(unlist(lapply(unique(df_tmp[['samplename']]), function(i){
        bw = stats::bw.nrd0(df_tmp[[meta]][df_tmp[['samplename']] == i])
        return(bw)
      }))), 2)
      p_list[[meta]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=.data[[meta]], fill=samplename)) +
        ggplot2::stat_density(color='gray30', alpha=ptalpha, bw=bandw, position="identity") +
        ggplot2::ylab(paste0('Density (Average bandwidth=', bandw, ')'))
    }
    p_list[[meta]] = p_list[[meta]] +
      ggplot2::xlab(NULL) +
      ggplot2::ggtitle(paste0('Variable: ', meta)) +
      ggplot2::scale_color_manual(values=meta_cols) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1),
                     panel.border=ggplot2::element_rect(fill=NA, color='black'),
                     legend.title=ggplot2::element_blank())
  }

  return(p_list)
}


##
# gene_plots
#

gene_plots = function(x=NULL, genes=NULL, samples=NULL, data_type='tr',
                      plot_type='violin', color_pal='roma', ptsize=0.5, ptalpha=0.5){
  # Extract data slot
  if(data_type == 'tr'){
    expr_tmp = x@tr_counts
    p_title = 'Normalized expression - '
    ax_title = 'Normalized expression'
  } else {
    expr_tmp = x@counts
    p_title = 'Raw expression - '
    ax_title = 'Counts'
  }

  p_list = list()
  # Loop through list of genes
  for(gene in genes){
    # Subset samples if gene not in all samples
    samples_tmp = samples
    df_tmp = tibble::tibble()
    for(i in samples){
      if(!(gene %in% rownames(expr_tmp[[i]]))){
        message(paste0('Sample ', i, ' does not contain ', gene, '.\n'))
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
    #p_list[[gene]] =
    if(plot_type == 'box'){
      p_list[[gene]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=samplename, y=.data[['geneexpr']], color=samplename)) +
        ggplot2::geom_boxplot(outlier.size=ptsize)
    } else if(plot_type == 'violin'){
      p_list[[gene]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=samplename, y=.data[['geneexpr']], color=samplename)) +
        ggplot2::geom_violin() +
        ggforce::geom_sina(size=ptsize)
    } else if(plot_type == 'density'){
      # Define bandwidth
      bandw = round(mean(unlist(lapply(unique(df_tmp[['samplename']]), function(i){
        bw = stats::bw.nrd0(df_tmp[['geneexpr']][df_tmp[['samplename']] == i])
        return(bw)
      }))), 2)
      p_list[[gene]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=.data[['geneexpr']], fill=samplename)) +
        ggplot2::stat_density(color='gray30', alpha=ptalpha, bw=bandw, position="identity") +
        ggplot2::ylab(paste0('Density (Average bandwidth=', bandw, ')'))
    }
    p_list[[gene]] = p_list[[gene]] +
      ggplot2::ggtitle(paste0(p_title, gene)) +
      ggplot2::ylab(ax_title) +
      ggplot2::xlab(NULL) +
      ggplot2::scale_color_manual(values=gene_cols) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1),
                     panel.border=ggplot2::element_rect(fill=NA, color='black'),
                     legend.title=ggplot2::element_blank())
  }

  return(p_list)
}


##
# cluster_gene_plots
#

cluster_gene_plots = function(x=NULL, plot_meta=NULL, genes=NULL, samples=NULL, data_type='tr',
                              plot_type='violin', color_pal='roma', ptsize=0.5, ptalpha=0.5){

  if(length(plot_meta) != 1){
    message(paste0('Only one metadata variable can be plotted at a time. Plotting ', plot_meta[1], '.'))
    plot_meta = plot_meta[1]
  }

  # Check if plot_meta is available and is a discrete variable. Remove sample, otherwise
  for(i in samples){
    if(!plot_meta %in% colnames(x@spatial_meta[[i]])){
      warning(paste0('Variable ', plot_meta, ' is not present in sample ', i, '.'))
      samples = grep(paste0('^', i, '$'), samples, value=T, invert=T)
    } else if(is.double(x@spatial_meta[[i]][[plot_meta]]) | length(as.character(unique(x@spatial_meta[[i]][[plot_meta]]))) > 30){
      warning(paste0('Variable ', plot_meta, ' seems continuous or has more than 30 categories and will be skipped for sample ', i, '.'))
      samples = grep(paste0('^', i, '$'), samples, value=T, invert=T)
    }
  }
  if(length(samples) < 1){
    stop('All samples skipped.')
  }

  # Check plot type
  if(!plot_type %in% c('box', 'violin')){
    message("Only 'box' or 'violin' are valid choices. Defaulting to 'box'.")
    plot_type = 'box'
  }

  # Extract data slot
  if(data_type == 'tr'){
    expr_tmp = x@tr_counts
    p_title = 'Normalized expression - '
    ax_title = 'Normalized expression'
  } else {
    expr_tmp = x@counts
    p_title = 'Raw expression - '
    ax_title = 'Counts'
  }

  p_list = list()
  # Loop through list of genes
  for(gene in genes){
    # Subset samples if gene not in all samples
    samples_tmp = samples
    df_tmp = tibble::tibble()
    for(i in samples){
      if(!(gene %in% rownames(expr_tmp[[i]]))){
        message(paste0('Sample ', i, ' does not contain ', gene, '.\n'))
        samples_tmp = grep(paste0('^', i, '$'), samples_tmp, value=T, invert=T)
      } else{
        df_tmp = dplyr::bind_rows(df_tmp,
                                  dplyr::left_join(tibble::tibble(libname=colnames(expr_tmp[[i]]),
                                                                  geneexpr=expr_tmp[[i]][gene, ],
                                                                  samplename=i),
                                                   x@spatial_meta[[i]][, c('libname', plot_meta)], by='libname'))
      }
    }
    # Force categories as character
    df_tmp[[plot_meta]] = as.character(df_tmp[[plot_meta]])

    # Define color palette
    gene_cols = color_parse(color_pal, n_cats=length(unique(df_tmp[[plot_meta]])))
    names(gene_cols) = unique(df_tmp[[plot_meta]])

    # Create plot
    p_list[[gene]] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=samplename, y=.data[['geneexpr']], color=.data[[plot_meta]]))
    if(plot_type == 'box'){
      p_list[[gene]] = p_list[[gene]] +
        ggplot2::geom_boxplot(outlier.size=ptsize)
    } else if(plot_type == 'violin'){
      p_list[[gene]] = p_list[[gene]] +
        ggplot2::geom_violin(ggplot2::aes(fill=.data[[plot_meta]])) +
        ggforce::geom_sina(size=ptsize) +
        ggplot2::scale_fill_manual(values=gene_cols)
    }
    p_list[[gene]] = p_list[[gene]] +
      ggplot2::scale_color_manual(values=gene_cols) +
      ggplot2::ggtitle(paste0(p_title, gene, '\n', plot_meta)) +
      ggplot2::ylab(ax_title) +
      ggplot2::xlab(NULL) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30, hjust=1),
                     panel.border=ggplot2::element_rect(fill=NA, color='black')#,
                     #legend.title=ggplot2::element_blank()
      )
  }

  return(p_list)
}

