##
#' @title plot_counts: Generates plots for the distribution of counts
#' @description Generates density plots, violin plots, and/or boxplots for the
#' distribution of count values
#' @details
#' The function allows to visualize the distribution counts across all genes and spots
#' in the STlist. The user can select between density plots, violin plots, or box
#' plots as visualization options. Useful for assessment of the effect of filtering and
#' data transformations and to assess zero-inflation. To plot counts or genes per
#' spot/cell, the function `distribution_plots` should be used instead.
#'
#' @param x an STlist
#' @param samples samples to include in the plot. Default (NULL) includes all samples
#' @param data_type one of `tr` or `raw`, to plot transformed or raw counts
#' @param plot_type one or several of `density`, `violin`, and `box`, to generate
#' density plots, violin plots, and/or boxplots
#' @param color_pal  a string of a color palette from `khroma` or `RColorBrewer`, or a
#' vector with colors
#' @param cvalpha the transparency of the density plots
#' @param distrib_subset the proportion of spots/cells to plot. Generating these
#' plots can be time consuming due to the large amount of elements to plot. This
#' argument provides control on how many randomly values to show to speed plotting
#' @param subset_seed related to `distrib_subset`. Sets the seed number to ensure
#' the same subset of values is selected for plotting
#' @return a list of ggplot objects
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
#'   cp <- plot_counts(melanoma, data_type='raw', plot_type=c('violin', 'box'))
#'   ggpubr::ggarrange(plotlist=cp)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#'
#' @import ggplot2
#
plot_counts = function(x=NULL, samples=NULL, data_type='tr', plot_type='density',
                       color_pal='okabeito', cvalpha=0.5, distrib_subset=0.5,
                       subset_seed=12345){

  # To prevent NOTES in R CMD check
  . = NULL

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
    # Define bandwidth
    bandw = round(mean(unlist(lapply(unique(df_tmp$sample_name), function(i){
      bw = stats::bw.nrd0(df_tmp[['expr_values']][df_tmp$sample_name == i])
      return(bw)
    }))), 2)

    d_plot[['density']] = ggplot2::ggplot(df_tmp, ggplot2::aes(x=expr_values, fill=sample_name)) +
#      ggplot2::geom_density(alpha=cvalpha) +
      ggplot2::stat_density(color='gray30', alpha=cvalpha, bw=bandw, position="identity") +
      ggplot2::xlab(ax_title) +
      ggplot2::ylab(paste0('Density (Bandwidth=', bandw, ')')) +
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

