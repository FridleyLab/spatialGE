##
#' @title plot_gene_quilt: Gene expression at each spot from an spatial array
#' @description Plot gene expression levels at each spot (quilt plot) within a
#' spatially-resolved transcriptomic array.
#' @details
#' This function produces a quilt plot for a series of HUGO gene names and spatial
#' arrays within an STList. The function can also plot tumor/stroma classifications
#' if ESTIMATE deconvolution results are present in the STList.
#'
#' @param x, an STList with voom-norm counts.
#' @param genes, a vector of one or more gene names to plot.
#' @param plot_who, a vector of numbers indicating the  spatial arrays to plot
#' genes from. If NULL, will plot all spatial arrays. The numbers in the vector
#' will match the order in which the arrays were specified when creating the STList.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param purity, logical, whether tumor/stroma classes should be plotted.
#' @param saveplot, logical, indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @export
#
#
plot_gene_quilt <- function(x = NULL, genes=NULL, plot_who=NULL, color_pal='YlOrBr',
                            purity=F, saveplot=F, scaled=F){

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
    stop("There are no normalized matrices in this STList.")
  }

  # Store maximum expression value in case 'scaled' is required.
  if(scaled){
    maxvalue <- c()
    for (i in plot_who) {
      for (gene in genes) {
        # Test if gene name exists in normalized count matrix.
        if (any(x@voom_counts[[i]]$gene == gene)) {
          # Find maximum expression value for each spatial array.
          values <- unlist(x@voom_counts[[i]][x@voom_counts[[i]]$gene == gene,][,-1])
          maxvalue <- append(maxvalue, max(values))
        }
      }
    }
    # Find maximum value among selected spatial arrays.
    maxvalue <- max(maxvalue)
  }

  # Loop through each normalized count matrix.
  for (i in plot_who) {

    # Create list of plots for a given subject.
    qp_list <- list()

    # Loop though genes to plot.
    for (gene in genes) {

      # Test if gene name exists in normalized count matrix.
      if (any(x@voom_counts[[i]]$gene == gene)) {
        # If scaled, standardize expression. If not, return expression values as they are.
        if(scaled){
          values <- unlist(x@voom_counts[[i]][x@voom_counts[[i]]$gene == gene,][,-1])/maxvalue
        } else{
          values <- unlist(x@voom_counts[[i]][x@voom_counts[[i]]$gene == gene,][,-1])
        }

        # Create a data frame with coordinates of spots and expression values.
        df <- dplyr::bind_cols(x@coords[[i]][,-1], tibble::as_tibble(values))
        colnames(df) <- c('y_pos', 'x_pos', 'values')

        # If purity=T, get tumor/stroma classifications and add to data frame.
        # Then call the quilt plot function.
        if(purity){
          if(!(rlang::is_empty(melanoma@cell_deconv))){
          df <- dplyr::bind_cols(df, cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
          qp <- quilt_p_purity2(data_f=df, leg_name="norm_expr", color_pal=color_pal,
                                title_name=paste0(gene, " - ", "subj ", i))
          } else{
            stop("No tumor/stroma classification in the STList.")
          }
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

    #     # Define number of columns and rows in plot and size.
    #     row_col <- c(2, 2)
    # #    w_pdf=9
    # #    h_pdf=9
    #     if(length(genes) == 2){
    #       row_col <- c(1, 2)
    # #      w_pdf=9
    # #      h_pdf=4.5
    #     }else if(length(genes) == 1){
    #       row_col <- c(1, 1)
    # #      w_pdf=7
    # #      h_pdf=7
    #     }

    row_col <- c(1, 1)

    # Test if plot should be saved to PDF.
    if(saveplot){
      # Print plots to PDF.
      pdf(file=paste0("gene_quilt_spatarray_", i, ".pdf"), width=10, height=10)#,
          #width=w_pdf, height=h_pdf
      print(ggpubr::ggarrange(plotlist=qp_list,
                              nrow=row_col[1], ncol=row_col[2]))
      dev.off()
    } else{
      # Print plots to console.
      for(p in 1:length(qp_list)){
        print(qp_list[[p]])
      }
    }

  }
}
