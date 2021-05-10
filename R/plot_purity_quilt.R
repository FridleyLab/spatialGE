##
#' @title plot_cell_quilt
#' @description Produces a quilt plot from cell scores.
#' @details
#' This function produces a quilt plot for a series cell names and spatial
#' arrays.
#'
#' @param x, an STList with deconvolution matrices containing the cell names.
#' @param cells, a vector of cell names within a deconvolution matrix (one or
#' several) to plot.
#' @param plot_who, a vector of subject indexes as ordered within the STList, to
#' plot cells from. If NULL, will plot for all subjects.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param saveplot, logical indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @export
#
#
plot_purity_quilt <- function(x = NULL, cells=NULL, plot_who=NULL,
                            color_pal='YlOrBr', saveplot=F){

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

  # Test if voom normalized counts are available.
  if (rlang::is_empty(x@cell_deconv)) {
    stop(paste("There are not deconvolution matrices in this STList."))
  }

  # if(!is.null(saveplot)){
  #   if(dir.exists(paste0(saveplot))){
  #     unlink(paste0(saveplot), recursive=T)
  #   }
  # }

  # Loop through each normalized count matrix.
  for (i in plot_who) {

    # Create list of plots for a given subject.
    qp_list <- list()

    # Loop though genes to plot.
#    for (cell in cells) {

#     if(cell != 'stroma_score'){
        # Test if gene name exists in normalized count matrix.
#        if (any(x@cell_deconv[[i]]$transf_deconv_matrix$cell_names == cell)) {
          # Create data frame of gene and plot.

          # "StromalScore"  "ImmuneScore"   "ESTIMATEScore"

          values <- unlist(x@cell_deconv[[i]]$estimate_purity[x@cell_deconv[[i]]$estimate_purity$NAME == 'ESTIMATEScore', ][,-1])
          df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
          colnames(df) <- c('x_pos', 'y_pos', 'values')

          # The color palette function in khroma is created by quilt_p() function.
          qp <- quilt_p(data_f=df, leg_name="score", color_pal=color_pal,
                        title_name=paste0('ESTIMATEScore', " - ", "subj ", i))
#        } else{
          # If gene is not in the matrix, move to next gene.
#          cat(paste(cell, "is not a cell type in the deconvoluted matrix."))
#          next
#        }
 #     }else if(cell == 'stroma_score'){
 #       if (any(x@cell_deconv[[i]]$transf_tumorstroma$cell_names == cell)) {
          # Create data frame of gene and plot.
#          values <- unlist(x@cell_deconv[[i]]$transf_tumorstroma[x@cell_deconv[[i]]$transf_tumorstroma$cell_names == cell,][,-1])
#         df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
#          colnames(df) <- c('x_pos', 'y_pos', 'values')

          # The color palette function in khroma is created by quilt_p() function.
#          qp <- quilt_p(data_f=df, leg_name="sqrt_score", color_pal=color_pal,
#                       title_name=paste0(cell, " - ", "subj ", i))
#        }
#      }

      # Append plot to list.
      qp_list[[cell]] <- qp
#    }

    # Define number of columns and rows in plot and size.
    row_col <- c(2, 2)
    w_pdf=9
    h_pdf=9
    if(length(cells) == 2){
      row_col <- c(1, 2)
      w_pdf=9
      h_pdf=4.5
    }else if(length(cells) == 1){
      row_col <- c(1, 1)
      w_pdf=7
      h_pdf=7
    }

    # Test if a filepath to save plots is available.
    if(saveplot){
      #dir.create(paste0(saveplot), recursive=T, showWarnings=F)
      pdf(file=paste0("cell_quilt_spat_array_", i, ".pdf"),
          width=w_pdf, height=h_pdf)
      print(ggpubr::ggarrange(plotlist=qp_list,
                              nrow=row_col[1], ncol=row_col[2]))
      dev.off()
    }else{
      # Print plots to console.
      print(ggpubr::ggarrange(plotlist=qp_list,
                              nrow=row_col[1], ncol=row_col[2]))
    }
  }
}
