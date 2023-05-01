##
#' @title count_distribution: Generates a density plot for the distribution of counts
#' @description Generates density plots, violin plots, or boxplots for the distribution of
#' counts per spot/cell in the samples
#' @details The function allows to visualize the distribution of spot/cell counts.
#' The user can select between density plots, violin plots, or box plots as visualization
#' options. Useful for assessment of the effect of data transformations
#'
#' @param x an STlist
#' @param samples samples to include in the plot. Default (NULL) includes all samples
#' @param data_type one of 'tr' or 'raw', to plot transformed or raw counts
#' @param color_pal  a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with colors
#' @param cvalpha the transparency of the density curves
#'
#' @export
#'
#' @import ggplot2
#
#
count_distribution = function(x=NULL, samples=NULL, data_type='tr', plot_type='density',
                              color_pal='okabeito', cvalpha=0.5, distrib_subset=0.5,
                              subset_seed=12345){
  #require('ggplot2')

  # Define samples to plot if NULL or numeric
  if(is.null(samples)){
    samples = names(x@counts)
  } else if(is.numeric(samples)){
    samples = names(x@counts)[samples]
  }

  # Check that transfromed counts are available if data_type='tr'
  if(data_type == 'tr'){
    if(rlang::is_empty(x@tr_counts)){
      stop('No transformed expression data in this STlist.')
    }
  }

  # Loop through samples and make long data frame
  df_tmp = list()
  for(i in samples){
    if(data_type == 'tr'){
      df_tmp[[i]] = as.matrix(x@tr_counts[[i]]) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(cols=dplyr::everything(.), names_to='libname', values_to='expr_values') %>%
        tibble::add_column(sample_name=i)
      # Subsample data according
      set.seed(subset_seed)
      df_tmp[[i]] = df_tmp[[i]][sample(1:nrow(df_tmp[[i]]), nrow(df_tmp[[i]])*distrib_subset), ]
      # Plot and X axis title
      p_title = 'Normalized expression distribution'
      ax_title = 'Normalized expression'
    } else{
      df_tmp[[i]] = as.matrix(x@counts[[i]]) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(cols=dplyr::everything(.), names_to='libname', values_to='expr_values') %>%
        tibble::add_column(sample_name=i)
      # Subsample data according
      set.seed(subset_seed)
      df_tmp[[i]] = df_tmp[[i]][sample(1:nrow(df_tmp[[i]]), nrow(df_tmp[[i]])*distrib_subset), ]
      # Plot and X axis title
      p_title = 'Raw expression distribution'
      ax_title = 'Raw counts'
    }
  }

  # Compile all counts from all samples
  df_tmp = do.call(dplyr::bind_rows, df_tmp)

  # Define color palette
  meta_cols = color_parse(color_pal, n_cats=length(samples))

  # Create plots
  d_plot = list()
  if(any(grepl('density', plot_type))){
    d_plot[['density']] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=expr_values, fill=sample_name)) +
      ggplot2::geom_density(alpha=0.5) +
      ggplot2::xlab(ax_title) +
      ggplot2::ylab('Density') +
      ggplot2::ggtitle(p_title) +
      ggplot2::labs(fill='Sample') +
      ggplot2::scale_fill_manual(values=meta_cols) +
      ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA, color='black'))
  }
  if(any(grepl('violin', plot_type))){
    d_plot[['violin']] = ggplot2::ggplot(df_tmp, ggplot2::aes(y=expr_values, x=sample_name, fill=sample_name)) +
      ggplot2::geom_violin() +
      #ggforce::geom_sina(size=0.1) +
      ggplot2::ggtitle(p_title) +
      ggplot2::ylab(ax_title) +
      ggplot2::xlab(NULL) +
      ggplot2::labs(fill='Sample') +
      ggplot2::scale_fill_manual(values=meta_cols) +
      ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA, color='black'),
                     axis.text.x=ggplot2::element_text(angle=30, vjust=1, hjust=0.8))
  }
  if(any(grepl('box', plot_type))){
    d_plot[['boxplot']] = ggplot2::ggplot(df_tmp, ggplot2::aes(y=expr_values, x=sample_name, fill=sample_name)) +
      ggplot2::geom_boxplot() +
      ggplot2::ggtitle(p_title) +
      ggplot2::ylab(ax_title) +
      ggplot2::xlab(NULL) +
      ggplot2::labs(fill='Sample') +
      ggplot2::scale_fill_manual(values=meta_cols) +
      ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA, color='black'),
                     axis.text.x=ggplot2::element_text(angle=30, vjust=1, hjust=0.8))
  }

  return(d_plot)
}

