##
# Plot the number of reads or number of genes (>=1 counts) per spot.
#
#
#
#
plot_QC_quilt <- function(x = NULL, nreads=F, ngenes=F, pctmt=F, plot_who=NULL, color_pal='YlOrBr',
                            purity=F, saveplot=NULL, inter=F, visium=T){

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Set nreads to default
  if(!any(c(nreads, ngenes, pctmt))){
    nreads = T
  }

  # Find maximum value among selected spatial arrays to scale color legends
  maxvalue = c()
  minvalue = c()
  list_vals = list()
  if(nreads){
    for (i in plot_who) {
      list_vals[[i]] = colSums(x@counts[[i]][, -1])
      # Find maximum expression value for each spatial array.
      maxvalue = append(maxvalue, max(list_vals[[i]]))
      minvalue = append(minvalue, min(list_vals[[i]]))
    }
  } else if(ngenes){
    for (i in plot_who) {
      list_vals[[i]] = x@counts[[i]][, -1] != 0
      list_vals[[i]] = colSums(list_vals[[i]])
      maxvalue = append(maxvalue, max(list_vals[[i]]))
      minvalue = append(minvalue, min(list_vals[[i]]))
    }
  } else if(pctmt){
    for (i in plot_who) {
      total_mtct = colSums(x@counts[[i]][grep("^MT-", x@counts[[i]][[1]]), -1])
      total_counts = colSums(x@counts[[i]][, -1])
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
    } else if(pctmt){
      leg_name = 'pct_MTreads'
    }

    # If purity=T, get tumor/stroma classifications and add to data frame.
    # Then call the quilt plot function.
    if(purity){
      if(!(rlang::is_empty(x@cell_deconv))){
        df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
        qp <- quilt_p_purity(data_f=df, leg_name=leg_name, color_pal=color_pal,
                             title_name=names(x@counts)[[i]],
                             minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        qpbw <- quilt_p_purity_bw(data_f=df, visium=visium, title_name=paste0('ESTIMATE\ntumor/stroma - ', names(x@counts)[[i]]))
      } else{
        stop("No tumor/stroma classification in the STList.")
      }
    }else{
      # The color palette function in khroma is created by quilt_p() function.
      qp <- quilt_p(data_f=df, leg_name=leg_name, color_pal=color_pal,
                    title_name=names(x@counts)[[i]],
                    minvalue=minvalue, maxvalue=maxvalue, visium=visium)
    }

    # Append plot to list.
    qp_list[[names(x@counts)[[i]]]] = qp

    if(purity){
      qp_list[[paste0(names(x@counts)[[i]], '_purity')]] = qpbw
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
