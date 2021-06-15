##
#' @title plot_deconv_quilt: Deconvolution scores at each spot from an spatial array
#' @description Plot cell type deconvolution scores at each spot (quilt plot) within a
#' spatially-resolved transcriptomic array.
#' @details
#' This function produces a quilt plot for a series cell names and spatial
#' arrays within an STList. The function can also plot tumor/stroma classifications
#' if ESTIMATE deconvolution results are present in the STList.
#'
#' @param x, an STList with deconvolution matrices containing the cell names.
#' @param cells, a vector of cell names within a deconvolution matrix (one or
#' several) to plot.
#' @param plot_who, a vector of subject indexes as ordered within the STList, to
#' plot cells from. If NULL, will plot for all subjects.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.'
#' @param purity, logical, whether tumor/stroma classes should be plotted.
#' @param saveplot, logical indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @param scaled, logical, indicating if expression values should be scaled with
#' respect to the highest value among all genes to plot. WARNING: Color legends
#' are not scaled between plots, but values are.
#' @param method, the deconvolutionn method to used estimate the cell scores.
#' @param visium, whether or not to reverse axes for Visium slides.
#' @return qplist, a list with plots.
#' @export
#
#
plot_deconv_quilt <- function(x = NULL, cells=NULL, plot_who=NULL, color_pal='YlOrBr',
                            purity=F, saveplot=F, scaled=F, method='xcell', visium=T){

  #  moran_est <- round(as.vector(x@gene_het[[gene]]$morans_I$estimate[[1]]), 2)
  #  geary_est <- round(as.vector(x@gene_het[[gene]]$gearys_C$estimate[[1]]), 2)
  #  getis_est <- round(as.vector(x@gene_het[[gene]]$getis_ord_Gi$estimate[[1]]), 2)

  # Test that a gene name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cells to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Test if deconvolution data of any type is available.
  if (rlang::is_empty(x@cell_deconv)) {
    stop("There are not deconvolution matrices in this STList.")
  }

  method=tolower(method)

  # Get requested list of deconvoluted matrices.
  if(method == 'xcell'){
    deconv_list <- x@cell_deconv$xCell
  } else if(method == 'ssgsea'){
    deconv_list <- x@cell_deconv$ssGSEA
  } else{
    stop('Please, specify a deconvolution method to plot.')
  }


  # Store maximum cell score value in case 'scaled' is required.
#  if(scaled){
    maxvalue <- c()
    minvalue <- c()
    for (i in plot_who) {
      for (cell in cells) {
        # Test if cell name exists in normalized count matrix.
        #if (any(x@cell_deconv$xCell[[i]]$sqrt_scores$cell_names == cell)) {
        if (any(deconv_list[[i]]$sqrt_scores$cell_names == cell)) {
          # Find maximum score value for each spatial array.
          # values <- unlist(x@cell_deconv$xCell[[i]]$sqrt_scores[x@cell_deconv$xCell[[i]]$sqrt_scores$cell_names == cell, ][,-1])
          values <- unlist(deconv_list[[i]]$sqrt_scores[deconv_list[[i]]$sqrt_scores$cell_names == cell, ][,-1])
          maxvalue <- append(maxvalue, max(values))
          minvalue <- append(minvalue, min(values))
        }
      }
    }
    # Find maximum value among selected spatial arrays.
    maxvalue <- max(maxvalue)
    minvalue <- min(minvalue)
#  }

  # Create list of plots.
  qp_list <- list()

  # if(purity){
  #   qpbw_list <- list()
  # }

  # Loop through each normalized count matrix.
  for (i in plot_who) {

    # Create list of plots for a given subject.
#    qp_list <- list()

    # Loop though genes to plot.
    for (cell in cells) {

      # Test if gene name exists in normalized count matrix.
      # if (any(x@cell_deconv$xCell[[i]]$sqrt_scores$cell_names == cell)) {
      if (any(deconv_list[[i]]$sqrt_scores$cell_names == cell)) {
        # If scaled, standardize expression. If not, return expression values as they are.
        if(scaled){
          #values <- unlist(x@cell_deconv$xCell[[i]]$sqrt_scores[x@cell_deconv$xCell[[i]]$sqrt_scores$cell_names == cell,][,-1])/maxvalue
          values <- unlist(deconv_list[[i]]$sqrt_scores[deconv_list[[i]]$sqrt_scores$cell_names == cell,][,-1])/maxvalue
        } else{
          #values <- unlist(x@cell_deconv$xCell[[i]]$sqrt_scores[x@cell_deconv$xCell[[i]]$sqrt_scores$cell_names == cell,][,-1])
          values <- unlist(deconv_list[[i]]$sqrt_scores[deconv_list[[i]]$sqrt_scores$cell_names == cell,][,-1])
        }

        # Create data frame of gene and plot.
        df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
        colnames(df) <- c('y_pos', 'x_pos', 'values')

        # If purity=T, get tumor/stroma classifications and add to data frame.
        # Then call the quilt plot function.
        if(purity){
          if(!(rlang::is_empty(x@cell_deconv))){
            df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
            qp <- quilt_p_purity(data_f=df, leg_name="sqrt_score", color_pal=color_pal,
                                  title_name=paste0(cell, " - ", "subj ", i),
                                 minvalue=minvalue, maxvalue=maxvalue, visium=visium)
            qpbw <- quilt_p_purity_bw(data_f=df, visium=visium, title_name=paste0('ESTIMATE\ntumor/stroma - subj ', i))
          } else{
            stop("No tumor/stroma classification in the STList.")
          }
        } else{
          # The color palette function in khroma is created by quilt_p() function.
          qp <- quilt_p(data_f=df, leg_name="sqrt_score", color_pal=color_pal,
                        title_name=paste0(cell, " - ", "subj ", i),
                        minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        }

      } else{
        # If gene is not in the matrix, move to next gene.
        cat(paste(cell, "is not a cell type in the deconvoluted matrix."))
        next
      }

      # Append plot to list.
      qp_list[[paste0(cell, "_", i)]] <- qp
    }

    if(purity){
      qp_list[[paste0('subj',i)]] <- qpbw
    }
    # # Define number of columns and rows in plot and size.
    # row_col <- c(2, 2)
    # w_pdf=9
    # h_pdf=9
    # if(length(cells) == 2){
    #   row_col <- c(1, 2)
    #   w_pdf=9
    #   h_pdf=4.5
    # }else if(length(cells) == 1){
    #   row_col <- c(1, 1)
    #   w_pdf=7
    #   h_pdf=7
    # }
  }
    row_col <- c(1, 1)

    # Test if a filepath to save plots is available.
    if(saveplot){
      #dir.create(paste0(saveplot), recursive=T, showWarnings=F)
      #pdf(file=paste0("cell_quilt_spatarray_", i, ".pdf"), width=10, height=10)#,
          #width=w_pdf, height=h_pdf
      pdf(file="cell_quilt_spatarray.pdf")
      print(ggpubr::ggarrange(plotlist=qp_list,
                              nrow=row_col[2], ncol=row_col[2],
                              common.legend=T, legend='bottom'))
      dev.off()
    }else{
      # Print plots to console.
      return(qp_list)
    }
#  }
}
