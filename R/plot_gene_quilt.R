##
#' @title plot_gene_quilt: Gene expression at each spot from an spatial array
#' @description Plot gene expression levels at each spot (quilt plot) within a
#' spatially-resolved transcriptomic array.
#' @details
#' This function produces a quilt plot for a series of HUGO gene names and spatial
#' arrays within an STList. The function can also plot tumor/stroma classifications
#' from ESTIMATE deconvolution results.
#'
#' @param x an STList with voom-norm counts.
#' @param genes a vector of one or more gene names to plot.
#' @param plot_who a vector of numbers indicating the  spatial arrays to plot
#' genes from. If NULL, will plot all spatial arrays. The numbers in the vector
#' will match the order in which the arrays were specified when creating the STList.
#' @param color_pal a color scheme from 'khroma' or RColorBrewer.
#' @param purity logical, whether tumor/stroma classes should be plotted.
#' @param saveplot logical, indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @param inter whether or not quilt plot should be converted to a Plotly plot.
#' @param visium whether or not to reverse axes for Visium slides.
#' @return a list with plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # qplots <- plot_gene_quilt(melanoma, genes='CD74', plot_who=2, purity=T, visium=F)
#'
#' @export
#
#
plot_gene_quiltV2 = function(x = NULL, trslot='log', genes=NULL, plot_who=NULL, color_pal='YlOrBr',
                            purity=F, saveplot=F, inter=F, visium=T){
  # Option to scale to 1 disabled.
  scaled=F

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Test if requested data slot is available.
  if(trslot == 'voom' & rlang::is_empty(x@voom_counts)) {
    stop("There are no voom transformed counts in this STList.")
  } else if(trslot == 'log' & rlang::is_empty(x@log_counts)){
    stop("There are no log transformed matrices in this STList.")
  } else if(trslot == 'raw' & rlang::is_empty(x@counts)){
    stop("There are no count matrices in this STList.")
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1 && genes == 'top'){
    genes = c()
    for(i in plot_who){
      genes = append(genes, x@gene_stdev[[i]]$gene[order(x@gene_stdev[[i]]$gene_stdevs, decreasing=T)][1:10])
    }
    # Get unique genes from most variable.
    genes = unique(genes)
  }

  # Select appropriate slot to take counts from
  if(trslot == 'voom'){
    counts = x@voom_counts
  }else if(trslot == 'log'){
    counts = x@log_counts
  }else if(trslot == 'raw'){
    counts = x@counts
  }

  # Store maximum expression value to standardize color legend.
  maxvalue <- c()
  minvalue <- c()
  for (i in plot_who) {
    for (gene in genes) {
      # Test if gene name exists in normalized count matrix.
      if (any(counts[[i]]$gene == gene)) {
        # Find maximum expression value for each spatial array.
        values <- unlist(counts[[i]][counts[[i]]$gene == gene,][,-1])
        maxvalue <- append(maxvalue, max(values))
        minvalue <- append(minvalue, min(values))
      }
    }
  }
  # Find maximum value among selected spatial arrays.
  maxvalue <- max(maxvalue)
  minvalue <- min(minvalue)

  # Create list of plots.
  qp_list <- list()

  # Loop through each normalized count matrix.
  for (i in plot_who) {

    # Loop though genes to plot.
    for (gene in genes) {

      # Test if gene name exists in normalized count matrix.
      if (any(counts[[i]]$gene == gene)) {
        # If scaled, standardize expression. If not, return expression values as they are.
        if(scaled){
          values <- unlist(counts[[i]][counts[[i]]$gene == gene,][,-1])/maxvalue
        } else{
          values <- unlist(counts[[i]][counts[[i]]$gene == gene,][,-1])
        }

        # Create a data frame with coordinates of spots and expression values.
        df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
        colnames(df) <- c('y_pos', 'x_pos', 'values')

        # If purity=T, get tumor/stroma classifications and add to data frame.
        # Then call the quilt plot function.
        if(purity){
          if(!(rlang::is_empty(x@cell_deconv))){
            df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
            qp <- quilt_p_purity(data_f=df, leg_name=paste0(trslot, "_expr"), color_pal=color_pal,
                                 title_name=paste0(gene, " - ", "sample ", names(counts[i])),
                                 minvalue=minvalue, maxvalue=maxvalue, visium=visium)
            qpbw <- quilt_p_purity_bw(data_f=df, visium=visium, title_name=paste0('ESTIMATE\ntumor/stroma - subj ', i))
          } else{
            stop("No tumor/stroma classification in the STList.")
          }
        }else{
          # The color palette function in khroma is created by quilt_p() function.
          qp <- quilt_p(data_f=df, leg_name=paste0(trslot, "_expr"), color_pal=color_pal,
                        title_name=paste0(gene, " - ", "sample ", names(counts[i])),
                        minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        }

      } else{
        # If gene is not in the matrix, move to next gene.
        cat(paste(gene, "is not a gene in the normalized count matrix."))
        next
      }

      # Append plot to list.
      qp_list[[paste0(gene, "_", names(counts[i]))]] <- qp
    }

    if(purity){
      qp_list[[paste0('subj', names(counts[i]))]] <- qpbw
    }
  }
  row_col <- c(1, 1)

  # Test if plot should be saved to PDF.
  if(saveplot){
    # Print plots to PDF.
    pdf(file="gene_quilt_spatarray.pdf")
    print(ggpubr::ggarrange(plotlist=qp_list, nrow=row_col[2], ncol=row_col[2], common.legend=T, legend='bottom'))
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
    #  }
  }
}
