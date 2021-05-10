##
#' @title plot_gene_quilt
#' @description Produces a quilt plot from ST data.
#' @details
#' This function produces a quilt plot for a series of HUGO gene names and spatial
#' arrays
#'
#' @param x, an STList with voom-norm counts.
#' @param genes, a vector of gene names (one or several) to plot.
#' @param plot_who, a vector of subject indexes as ordered within the STList, to
#' plot genes from. If NULL, will plot for all subjects.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param saveplot, logical indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @export
#
#
plot_gene_quilt <- function(x = NULL, genes=NULL, plot_who=NULL,
                            color_pal='YlOrBr', purity=F, saveplot=F){

  #  moran_est <- round(as.vector(x@gene_het[[gene]]$morans_I$estimate[[1]]), 2)
  #  geary_est <- round(as.vector(x@gene_het[[gene]]$gearys_C$estimate[[1]]), 2)
  #  getis_est <- round(as.vector(x@gene_het[[gene]]$getis_ord_Gi$estimate[[1]]), 2)

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Test if voom normalized counts are available.
  if (rlang::is_empty(x@voom_counts)) {
    stop(paste("There are not normalized matrices in this STList."))
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
    for (gene in genes) {

      # Test if gene name exists in normalized count matrix.
      if (any(x@voom_counts[[i]]$gene == gene)) {
        # Create data frame of gene and plot.
        values <- unlist(x@voom_counts[[i]][x@voom_counts[[i]]$gene == gene,][,-1])
        df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
        colnames(df) <- c('x_pos', 'y_pos', 'values')

        if(purity){
          df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
          qp <- quilt_p_purity(data_f=df, leg_name="norm_expr", color_pal=color_pal,
                               title_name=paste0(gene, " - ", "subj ", i))
        }else{
        # The color palette function in khroma is created by quilt_p() function.
          qp <- quilt_p(data_f=df, leg_name="norm_expr", color_pal=color_pal,
                        title_name=paste0(gene, " - ", "subj ", i))
        }

      } else{
        # If gene is not in the matrix, move to next gene.
        cat(paste(gene, "is not a gene in the normalized count matrix."))
        next
      }

      # Append plot to list.
      qp_list[[gene]] <- qp
    }

    # Define number of columns and rows in plot and size.
    row_col <- c(2, 2)
#    w_pdf=9
#    h_pdf=9
    if(length(genes) == 2){
      row_col <- c(1, 2)
#      w_pdf=9
#      h_pdf=4.5
    }else if(length(genes) == 1){
      row_col <- c(1, 1)
#      w_pdf=7
#      h_pdf=7
    }

    # Test if a filepath to save plots is available.
    if(saveplot){
      #dir.create(paste0(saveplot), recursive=T, showWarnings=F)
      pdf(file=paste0("gene_quilt_spat_array_", i, ".pdf")#,
          #width=w_pdf, height=h_pdf
          )
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
