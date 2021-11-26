##
#' @title plot_QC_quilt: Plot UMI counts at each spot from an spatial array
#' @description Plots total counts, gene counts, or counts percentages at each spot
#' within a spatially-resolved transcriptomic array.
#' @details
#' This function produces a quilt plot for total counts, gene/UMI counts, or percentage
#' of counts of a gene name matching a expression. This function is usfeul to explore
#' quality of the spots.
#'
#' @param x an STList
#' @param nreads logical, if TRUE, plots the total counts per spot
#' @param ngenes logical, if TRUE, plots the number of UMI/genes per spot
#' @param pct_Expr logical, if TRUE, plots the percentage of counts from genes with
#' names matching a expression given by `spot_pctExpr`
#' @param plot_who, a vector of numbers indicating the spatial arrays to plot.Numbers
#' follow the order of `names(x@counts)`. If NULL, will plot all spatial arrays.
#' @param color_pal, a color scheme from 'khroma' or RColorBrewer.
#' @param purity logical, whether tumor/stroma classes should be plotted.
#' @param image logical, whether to print the image stored for the spatial arrays
#' @param saveplot, a file name specifying the name of a PDF file to write plots to.
#' @param inter, whether or not quilt plot should be converted to a Plotly plot.
#' @param visium, logical, whether or not the samples are from a Visium experiment.
#' @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
#' @return a list with plots.
#'
#' @export
#
#
plot_QC_quilt = function(x = NULL, nreads=F, ngenes=F, pct_Expr=F, plot_who=NULL, spot_pctExpr='^MT-',
                         color_pal='YlOrBr', purity=F, saveplot=NULL, inter=F, visium=T, ptsize=0.5){
  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Set nreads to default
  if(!any(c(nreads, ngenes, pct_Expr))){
    nreads = T
  }

  # Find maximum value among selected spatial arrays to scale color legends
  maxvalue = c()
  minvalue = c()
  list_vals = list()
  if(nreads){
    for (i in plot_who) {
      df_current = expandSparse(x@counts[[i]])
      list_vals[[i]] = colSums(df_current)
      # Find maximum expression value for each spatial array.
      maxvalue = append(maxvalue, max(list_vals[[i]]))
      minvalue = append(minvalue, min(list_vals[[i]]))
    }
  } else if(ngenes){
    for (i in plot_who) {
      df_current = expandSparse(x@counts[[i]])
      list_vals[[i]] = df_current != 0
      list_vals[[i]] = colSums(list_vals[[i]])
      maxvalue = append(maxvalue, max(list_vals[[i]]))
      minvalue = append(minvalue, min(list_vals[[i]]))
    }
  } else if(pct_Expr){
    for (i in plot_who) {
      df_current = expandSparse(x@counts[[i]])
      total_mtct = colSums(df_current[grep(spot_pctExpr, rownames(df_current)), ])
      total_counts = colSums(df_current)
      list_vals[[i]] = total_mtct/total_counts
      maxvalue = append(maxvalue, max(list_vals[[i]]))
      minvalue = append(minvalue, min(list_vals[[i]]))
    }
  }
  maxvalue <- max(maxvalue)
  minvalue <- min(minvalue)

  # Create list of plots.
  qp_list <- list()

  # Loop through each normalized count matrix.
  for (i in plot_who) {

    df = dplyr::bind_cols(x@coords[[i]], values=list_vals[[i]])
    colnames(df) <- c('libname', 'y_pos', 'x_pos', 'values')

    # Define legend title.
    if(nreads){
      leg_name = 'reads'
    } else if(ngenes){
      leg_name = 'genes/UMIs'
    } else if(pct_Expr){
      leg_name = 'pct_Expr'
    }

    # If purity=T, get tumor/stroma classifications and add to data frame.
    # Then call the quilt plot function.
    if(purity){
      if(!(rlang::is_empty(x@cell_deconv))){
        df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
        qp <- quilt_p_purity(data_f=df, leg_name=leg_name, color_pal=color_pal,
                             title_name=names(x@counts)[[i]], ptsize=ptsize,
                             minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        qpbw <- quilt_p_purity_bw(data_f=df, visium=visium,
                                  title_name=paste0('ESTIMATE\ntumor/stroma - ', names(x@counts[i])),
                                  ptsize=ptsize)
      } else{
        stop("No tumor/stroma classification in the STList.")
      }
    }else{
      # The color palette function in khroma is created by quilt_p() function.
      qp <- quilt_p(data_f=df, leg_name=leg_name, color_pal=color_pal,
                    title_name=names(x@counts)[[i]],
                    minvalue=minvalue, maxvalue=maxvalue, visium=visium, ptsize=ptsize)
    }

    # Append plot to list.
    qp_list[[names(x@counts[i])]] = qp

    if(purity){
      qp_list[[paste0(names(x@counts[i]), '_purity')]] = qpbw
    }
  }
  row_col <- c(1, 1)

  # Test if plot should be saved to PDF.
  if(!is.null(saveplot)){
    # Print plots to PDF.
    pdf(file=saveplot)
    print(qp_list)
    dev.off()
  } else{
    # Convert ggplots to plotly plots.
    if(inter == T && purity == T){
      for(p in 1:length(qp_list)){
        qp_list[[p]] = plotly::ggplotly(qp_list[[p]], tooltip=c('Y', 'X', 'colour', 'shape'))
      }
      return(qp_list)
    }else{
      # Print plots to console.
      return(qp_list)
    }
  }
}
