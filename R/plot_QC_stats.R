##
#' @title plot_QCstats: Plots distribution of counts, gene number, or percents per sample
#' @description Generates plots with common quality control stats
#' @details
#' This function returns a list of plots with one of three QC stats:
#' \itemize{
#' \item Total number of reads per spot
#' \item Total number of genes per spot
#' \item Percent of counts of genes matching `spot_pctExpr` per spot
#'}
#'
#' @param x an STList
#' @param qc one of `total_reads`, `total_genes`, or `percent_expr.` The stat to be plotted
#' @param plot_who index from `x@counts``of the samples to be plotted
#' @param spot_pctExpr a expression to match gene names to calculate a percentage. By default '^MT-', to match mtDNA genes
#' @param ptsize the size of points in the plots
#' @param ptalpha the transparency of points in plots
#' @param ptjitter the dispersion of points in the plots
#' @param ylim a vector of two integers giving the minimum and maximum values of the y axis
#' @param saveplot File path to save a pdf with plots
#'
#' @export
#'
#' @import ggplot2
#
#
plot_QC_stats = function(x=NULL, qc='percent_expr', plot_who=NULL, spot_pctExpr='^MT-', ptsize=0.5, ptalpha=0.5, ptjitter=0.2, saveplot=NULL, ylim=NULL){

  #require('ggplot2')

  # Define samples to plot if NULL
  if(is.null(plot_who)){
    plot_who = 1:length(x@counts)
  }

  # Limit qc at one at a time
  if(length(qc) != 1){
    stop('Please, specify only one of total_reads, total_genes, or percent_expr')
  }

  # Create list of data frames with stats per spot
  # Store maximum value to scale axis
  max_val = c()
  min_val = c()
  df_list = list()
  for(i in plot_who){
    # Decompress matrix
    df = expandSparse(x@counts[[i]])

    # Get total read counts
    total_reads = colSums(df)

    if(qc == 'total_reads'){
      df_list[[names(x@counts[i])]] = dplyr::bind_cols('sample'=rep(names(x@counts[i]), ncol(df)), 'total_reads'=total_reads)
      max_val = append(max_val, max(df_list[[names(x@counts[i])]][[2]]))
      min_val = append(min_val, min(df_list[[names(x@counts[i])]][[2]]))
    }
    if(qc == 'total_genes'){
      nonzero_df = df != 0
      total_genes = colSums(nonzero_df)
      df_list[[names(x@counts[i])]] = dplyr::bind_cols('sample'=rep(names(x@counts[i]), ncol(df)), 'total_genes'=total_genes)
      max_val = append(max_val, max(df_list[[names(x@counts[i])]][[2]]))
      min_val = append(min_val, min(df_list[[names(x@counts[i])]][[2]]))
    }
    if(qc == 'percent_expr'){
      total_expr = colSums(df[grep(spot_pctExpr, rownames(df)), ])
      df_list[[names(x@counts[i])]] = dplyr::bind_cols('sample'=rep(names(x@counts[i]), ncol(df)), 'percent_expr'=total_expr/total_reads)
      max_val = append(max_val, max(df_list[[names(x@counts[i])]][[2]]))
      min_val = append(min_val, min(df_list[[names(x@counts[i])]][[2]]))
    }
  }

  # Create list with plots
  qcplot = list()
  for(i in plot_who){
    qcplot[[names(x@counts[i])]] = ggplot(df_list[[names(x@counts[i])]], aes(x=sample, y=get(qc))) +
      geom_violin() +
      geom_point(position=position_jitter(seed=1, width=ptjitter), size=ptsize, alpha=ptalpha)
    # Check if a different y limit has been requested
    if(!is.null(ylim)){
      if(length(ylim) == 2){
        qcplot[[names(x@counts[i])]] = qcplot[[names(x@counts[i])]] + ylim(ylim[1], ylim[2])
      }
    } else{
      qcplot[[names(x@counts[i])]] = qcplot[[names(x@counts[i])]] + ylim(min(min_val), max(max_val))
    }
    qcplot[[names(x@counts[i])]] = qcplot[[names(x@counts[i])]] +
      xlab('') +
      ylab(qc) +
      ggtitle(names(x@counts[i])) +
      scale_x_discrete(breaks=NULL) +
      theme(axis.title.x=element_blank())
  }

  # Generate grid of plots
  if(length(plot_who) > 3){
    mfrow_n = n2mfrow(length(plot_who))
  } else{
    mfrow_n = c(1, length(plot_who))
  }
  qc_pgrid = ggpubr::ggarrange(plotlist=qcplot, ncol=mfrow_n[2], nrow=mfrow_n[1])

  # Save plots to pdf if requested
  if(is.null(saveplot)){
    return(qc_pgrid)
  } else{
    ggsave(saveplot, plot=qc_pgrid)
  }
}

