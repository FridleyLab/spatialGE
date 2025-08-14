##
#' @title STdiff_volcano: Generates volcano plots from STdiff results
#' @description Generates volcano plots of differential expression results from STdiff
#' @details The function generated volcano plots (p-value vs. log-fold change) for
#' genes tested with `STdiff`. Colors can be customized to show significance from
#' spatial and non-spatial models
#'
#' @param x the output of `STdiff`
#' @param samples samples to create plots
#' @param clusters names of the clusters to generate comparisons
#' @param pval_thr the p-value threshold to color genes with differential expression
#' @param color_pal the palette to color genes by significance
#' @return a list of ggplot objects
#'
#' @export
#'
#' @import ggplot2
#
#
STdiff_volcano = function(x=NULL, samples=NULL, clusters=NULL, pval_thr=0.05, color_pal=NULL){
  # Define samples to plot if NULL or numeric
  if(is.null(samples)){
    samples = names(x)
  } else if(is.numeric(x)){
    samples = names(x)[samples]
  }

  # Loop through samples
  plist = list()
  for(i in samples){
    df_tmp = x[[i]]
    # Find out if test was pairwise
    pairwise = F
    if(any(grepl('^cluster_2$', colnames(x[[i]])))){
      pairwise = T
    }
    # Define which clusters to plot
    if(!is.null(clusters)){
      if(pairwise){
        df_tmp = df_tmp %>% dplyr::filter(cluster_1 %in% clusters | cluster_2 %in% clusters)
      } else{
        df_tmp = df_tmp %>% dplyr::filter(cluster_1 %in% clusters)
      }
    }

    # Define combination of clusters to plot
    combo = df_tmp %>% dplyr::select(dplyr::any_of(c('cluster_1', 'cluster_2'))) %>% dplyr::distinct()
    # Create plots
    for(cb in 1:nrow(combo)){
      cl1 = as.vector(unlist(combo[cb, 'cluster_1']))
      pl_name = paste0(i, '_', cl1)
      pl_title = paste0('Sample: ', i, '\nCluster ', cl1)
      if(pairwise){
        cl2 = as.vector(unlist(combo[cb, 'cluster_2']))
        pl_name = paste0(i, '_', cl1, '_vs_', cl2)
        pl_title = paste0('Sample: ', i, '\nCluster ', cl1, ' vs. ', cl2)
        #df_plot = df_tmp %>% dplyr::filter(cluster_1 == cl1 | cluster_2 == cl2)
        df_plot = df_tmp %>% dplyr::filter( (cluster_1 == cl1 & cluster_2 == cl2) | (cluster_1 == cl2 & cluster_2 == cl1) )
      } else{
        df_plot = df_tmp %>% dplyr::filter(cluster_1 == cl1)
      }

      # Detect test performed
      pval_col = grep("mm_p_val|ttest_p_val|wilcox_p_val", colnames(df_plot), value=T)
      if(pval_col == 'mm_p_val'){
        testname = 'Mixed models'
      } else if(pval_col == 'ttest_p_val'){
        testname = 't-test'
      } else if(pval_col == 'wilcox_p_val'){
        testname = "Wilcoxon's test"
      }

      # Classify p-values
      df_plot = df_plot %>%
        dplyr::mutate(signif=dplyr::case_when(adj_p_val < pval_thr & avg_log2fc < 0 ~ 'Down',
                                              adj_p_val < pval_thr & avg_log2fc > 0 ~ 'Up',
                                              TRUE ~ 'Not DE')) %>%
        dplyr::rename(plot_pval := !!pval_col)

      # Check if results for spatial model tests are available
      if(any(colnames(df_plot) == 'exp_adj_p_val')){
        df_plot = df_plot %>%
          dplyr::mutate(signif=dplyr::case_when(exp_adj_p_val < pval_thr & avg_log2fc < 0 ~ 'Down (spatial)',
                                                exp_adj_p_val < pval_thr & avg_log2fc > 0 ~ 'Up (spatial)',
                                                TRUE ~ signif))
      }

      # Re-order values in significance legend
      df_plot = df_plot %>%
        dplyr::mutate(signif=factor(signif, levels=c('Not DE', 'Down', 'Up', 'Down (spatial)', 'Up (spatial)')))

      # Define color palette
      if(is.null(color_pal)){
        cat_cols = c('black', 'cornflowerblue', 'lightpink', 'blue', 'red')
        names(cat_cols) = c('Not DE', 'Down',  'Up', 'Down (spatial)', 'Up (spatial)')
      } else{
        if(length(color_pal) == 3){
          cat_cols = color_pal
          names(cat_cols) = c('Not DE', 'Down',  'Up')
        } else if(length(color_pal) == 5){
          cat_cols = color_pal
          names(cat_cols) = c('Not DE', 'Down',  'Up', 'Down (spatial)', 'Up (spatial)')
        } else{
          message('Custom color palettes must have three (non-spatial results) or five colors (spatial results). Setting to default palette.')
          cat_cols = c('black', 'cornflowerblue', 'lightpink', 'blue', 'red')
          names(cat_cols) = c('Not DE', 'Down',  'Up', 'Down (spatial)', 'Up (spatial)')
        }
      }

      plist[[pl_name]] = ggplot2::ggplot(df_plot, aes(x=avg_log2fc, y=-log10(plot_pval))) +
        ggplot2::geom_point(aes(color=signif)) +
        ggrepel::geom_text_repel(aes(label=gene), size=3, verbose=F, force=0, max.iter=100, max.overlaps=1, nudge_x=0.2) +
        ggplot2::xlab('Average log fold-change') +
        ggplot2::ylab(paste0('-log10(nominal p-value)\n', testname)) +
        labs(color='FDR\nsignif.') +
        ggplot2::ggtitle(pl_title) +
        ggplot2::scale_color_manual(values=cat_cols[ names(cat_cols) %in% unique(df_plot[['signif']]) ]) +
        ggplot2::theme(panel.background=element_rect(color='black', fill=NULL))

      rm(list=grep("df_plot|cl1|cl2|pl_name|pl_title|testname|pval_col", ls(), value=T)) # Clean env
    }
    rm(combo, df_tmp, pairwise) # Clean env
  }
  return(plist)
}

