##
#' @title plot_gene_quilt: Gene expression at each spot from an spatial array
#' @description Plot raw or transformed gene expression levels at each spot (quilt plot)
#' within a spatially-resolved transcriptomic array.
#' @details
#' This function produces a quilt plot for a series of gene names and spatial
#' arrays within an STList. The function can also plot tumor/stroma classifications
#' from ESTIMATE deconvolution results.
#'
#' @param x, an STList
#' @param genes, a vector of one or more gene names to plot.
#' @param plot_who, a vector of numbers indicating the spatial arrays to plot
#' genes from. Numbers follow the order of `names(x@counts)`. If NULL, will plot
#' all spatial arrays.
#' @param color_pal, a color scheme from 'khroma' or RColorBrewer.
#' @param data_type, one of 'tr' or 'raw', to plot transformed or raw counts
#' respectively.
#' @param purity logical, whether tumor/stroma classes should be plotted.
#' @param image logical, whether to print the image stored for the spatial arrays
#' @param saveplot, a file name specifying the name of a PDF file to write plots to.
#' @param inter, whether or not quilt plot should be converted to a Plotly plot.
#' @param visium, logical, whether or not the samples are from a Visium experiment.
#' @ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
#' @return a list with plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # qplots <- plot_gene_quilt(melanoma, genes='CD74', plot_who=2, purity=T, visium=F)
#'
#' @export
#
#
plot_gene_quilt = function(x=NULL, genes=NULL, plot_who=NULL, color_pal='YlOrBr', data_type='tr', purity=F, image=F, saveplot=NULL, inter=F, visium=T, ptsize=NULL){

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Select appropriate slot to take counts from
  if(data_type == 'tr'){
    counts = x@tr_counts
  }else if(data_type == 'raw'){
    counts = x@counts
    # Expand sparse matrices
    for(i in plot_who){
      counts[[i]] = expandSparse(counts[[i]])
      counts[[i]] = tibble::as_tibble(tibble::rownames_to_column(counts[[i]], var='gene'))
    }
  } else(
    stop('Please, select one of transformed (tr) or raw (raw) counts')
  )

  # Test if requested data slot is available.
  if(rlang::is_empty(counts)) {
    stop("Data was not found in the specified slot of this STList.")
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1 && genes == 'top'){
    genes = c()
    for(i in plot_who){
      genes = append(genes, x@gene_var[[i]]$gene[order(x@gene_var[[i]]$gene_stdevs, decreasing=T)][1:10])
    }
    # Get unique genes from most variable.
    genes = unique(genes)
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

  # Define size of points
  if(is.null(ptsize)){
    ptsize = 0.5
  }

  # Create list of plots.
  qp_list <- list()
  # Loop through each normalized count matrix.
  for (i in plot_who) {

    # Loop though genes to plot.
    for (gene in genes) {

      # Test if gene name exists in normalized count matrix.
      if (any(counts[[i]]$gene == gene)) {

        values <- unlist(counts[[i]][counts[[i]]$gene == gene,][,-1])

        # Create a data frame with coordinates of spots and expression values.
        df <- dplyr::bind_cols(x@coords[[i]][, -1], tibble::as_tibble(values))
        colnames(df) <- c('y_pos', 'x_pos', 'values')

        # If purity=T, get tumor/stroma classifications and add to data frame.
        # Then call the quilt plot function. Also get appropriate legend title
        if(data_type == 'tr'){
          qlegname = paste0(x@misc[['transform']], "_expr")
        } else{
          qlegname = 'RawNorm_expr'
        }

        if(purity){
          if(!(rlang::is_empty(x@cell_deconv))){
            df = dplyr::bind_cols(df, cluster=x@cell_deconv[['ESTIMATE']][[i]]$purity_clusters$cluster)
            qp = quilt_p_purity(data_f=df, leg_name=qlegname, color_pal=color_pal,
                                title_name=paste0(gene, "\n", "sample: ", names(counts[i])),
                                minvalue=minvalue, maxvalue=maxvalue, visium=visium, ptsize=ptsize)
            qpbw <- quilt_p_purity_bw(data_f=df, visium=visium,
                                      title_name=paste0('ESTIMATE\ntumor/stroma\nsample: ', names(counts[i])),
                                      ptsize=ptsize)
          } else{
            stop("No tumor/stroma classification in the STList.")
          }
        }else{
          # The color palette function in khroma is created by quilt_p() function.
          qp <- quilt_p(data_f=df, leg_name=qlegname, color_pal=color_pal,
                        title_name=paste0(gene, "\n", "sample: ", names(counts[i])),
                        minvalue=minvalue, maxvalue=maxvalue, visium=visium, ptsize=ptsize)
        }

      } else{
        # If gene is not in the matrix, move to next gene.
        cat(paste(gene, "is not a gene in the normalized count matrix."))
        next
      }

      # Append plot to list.
      qp_list[[paste0(gene, "_", names(counts[i]))]] <- qp
    }

    if(image && !is.null(x@misc[['sp_images']][[i]])){
      img_obj = grid::rasterGrob(x@misc[['sp_images']][[i]])
      qp_list[[paste0('image', names(counts[i]))]] = ggplot() +
        annotation_custom(img_obj)
    }

    if(purity){
      qp_list[[paste0('tumorstroma_', names(counts[i]))]] <- qpbw
    }
  }

  row_col <- c(1, 1)

  # Test if plot should be saved to PDF.
  if(!is.null(saveplot)){
    # Print plots to PDF.
    pdf(file=saveplot)
    print(ggpubr::ggarrange(plotlist=qp_list, nrow=row_col[1], ncol=row_col[2], common.legend=F, legend='bottom'))
    dev.off()
  } else{
    # Convert ggplots to plotly plots.
    if(inter == T && purity == T && image == F){
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

