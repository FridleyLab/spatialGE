##
#' @title density_plots: Generates a density plot for the distribution of counts
#' @description Generates density plots for the distribution of counts per spot/cell in the samples
#' @details The function allows to visualize the distribution of spot/cell counts. Useful for
#' assessment of the effect of data transformations
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
density_plots = function(x=NULL, samples=NULL, data_type='tr', color_pal='okabeito',
                         cvalpha=0.5, dist_subset=0.5, subset_seed=12345){
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
      df_tmp[[i]] = df_tmp[[i]][sample(1:nrow(df_tmp[[i]]), nrow(df_tmp[[i]])*dist_subset), ]
      # X axis title
      p_title = 'normalized expression'
    } else{
      df_tmp[[i]] = as.matrix(x@counts[[i]]) %>%
        as.data.frame() %>%
        tidyr::pivot_longer(cols=dplyr::everything(.), names_to='libname', values_to='expr_values') %>%
        tibble::add_column(sample_name=i)
      # Subsample data according
      set.seed(subset_seed)
      df_tmp[[i]] = df_tmp[[i]][sample(1:nrow(df_tmp[[i]]), nrow(df_tmp[[i]])*dist_subset), ]
      # X axis title
      p_title = 'raw expression'
    }
  }

  # Compile all counts from all samples
  df_tmp = do.call(dplyr::bind_rows, df_tmp)

  # Define color palette
  meta_cols = color_parse(color_pal, n_cats=length(samples))

  # Create plot
  d_plot = ggplot2::ggplot(df_tmp, ggplot2::aes(x=expr_values, fill=sample_name)) +
    ggplot2::geom_density() +
    ggplot2::xlab(p_title) +
    ggplot2::labs(fill='Sample') +
    ggplot2::scale_fill_manual(values=meta_cols) +
    ggplot2::theme(panel.border=ggplot2::element_rect(fill=NA, color='black'))

  return(p_list)
}

