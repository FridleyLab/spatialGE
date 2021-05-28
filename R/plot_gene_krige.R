##
#' @title plot_gene_krige: Visualize transcriptomic surfaces
#' @description Produces a transcriptomic surface from kriging interpolation of ST data.
#' @details
#' This function produces a transcriptomic surface plot for a series of HUGO gene names and
#' spatial arrays.
#'
#' @param x, an STList with kriging objects for the genes selected.
#' @param genes, a vector of gene names (one or several) to plot.
#' @param krige_type, either 'ord' (ordinary; default), or 'univ' (universal)
#' kriging. Data for the respective kriging must be generated previously with
#' cell_krige().
#' @param plot_who, a vector of subject indexes as ordered within the STList, to
#' plot genes from. If NULL, will plot all subjects.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param purity, logical, whether or not annotate tumor and stroma spots based on
#' ESTIMATE tumor purity scores.
#' @param saveplot, logical, indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @param scaled, logical, indicating if expression values should be scaled with
#' respect to the highest value among all genes to plot. WARNING: Color legends
#' are not scaled between plots, but values are.
#' @export
#
#
plot_gene_krige <- function(x=NULL, genes=NULL, krige_type='ord', plot_who=NULL,
                            color_pal='YlOrBr', purity=F, saveplot=F, scaled=F){

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Store maximum expression value in case 'scaled' is required.
  if(scaled){
    maxvalue <- c()
    for (i in plot_who) {
      for (gene in genes) {
        # Test if kriging exists for a gene and subject.
        if (rlang::has_name(x@gene_krige, gene)){
          if(length(x@gene_krige[[gene]]) >= i){
            if (rlang::has_name(x@gene_krige[[gene]][[i]], krige_type)){
              # Find maximum expression value for each spatial array.
              values <- x@gene_krige[[gene]][[i]][[krige_type]]$predict
              maxvalue <- append(maxvalue, max(values))
            }
          }
        }
        # Find maximum value among selected spatial arrays.
        maxvalue <- max(maxvalue)
      }
    }
  }

  # Loop through each of the subjects.
  for (i in plot_who) {

    # Create list of plots for a given subject.
    kp_list <- list()

    # Loop though genes to plot.
    for (gene in genes) {

      if(length(x@gene_krige[[gene]]) >= i){
        if(is.null(x@gene_krige[[gene]][[i]])){
          cat(paste0(gene, " kriging for subject ", i, " is not present in STList\n"))
          next
        }
      }else{
        cat(paste0(gene, " kriging for subject ", i, " is not present in STList\n"))
        next
      }

      # Find prediction grid.
      if(krige_type == 'univ'){
        predict_grid <- x@gene_krige[[gene]][[i]]$univ_grid
      }else if(krige_type == 'ord'){
        predict_grid <- x@gene_krige[[gene]][[i]]$ord_grid
      }

      # Create data frame with coordinates and kriging values.
      if(scaled){
        krige_vals <- (x@gene_krige[[gene]][[i]][[krige_type]]$predict)/maxvalue
      } else{
        krige_vals <- x@gene_krige[[gene]][[i]][[krige_type]]$predict
      }

      df <- dplyr::bind_cols(predict_grid, krige=krige_vals)
      names(df) <- c("x_pos", "y_pos", "krige")

      # Get coordinates of bounding box enclosing the predicted grid.
      bbox<-rbind(
        c(min(predict_grid$Var1)-1, max(predict_grid$Var2)+1),
        c(max(predict_grid$Var1)+1, max(predict_grid$Var2)+1),
        c(max(predict_grid$Var1)+1, min(predict_grid$Var2)-1),
        c(min(predict_grid$Var1)-1, min(predict_grid$Var2)-1)
      )
      bbox <- as.data.frame(bbox)
      names(bbox) <- c("V1", "V2")

      # Create SpatialPolygon object with the bounding box.
      bbox_sp <- sp::SpatialPolygons(
        list(sp::Polygons(list(sp::Polygon(bbox)), "id")))

      # Create Spatial Polygon with the inner tissue border (concave hull)
      mask_sp <- sp::SpatialPolygons(
        list(sp::Polygons(list(sp::Polygon(x@prediction_border[[i]])), "id")))

      # Substract the concave hull from the bounding box, yielding a SpatialPolygon
      # object. The convert to
      bbox_mask_diff <- raster::erase(bbox_sp, mask_sp)

      # Get spatial statistics.
      moran_est <- round(as.vector(x@gene_het[[gene]][[i]]$morans_I$estimate[[1]]), 2)
      geary_est <- round(as.vector(x@gene_het[[gene]][[i]]$gearys_C$estimate[[1]]), 2)
      getis_est <- round(as.vector(x@gene_het[[gene]][[i]]$getis_ord_Gi$estimate[[1]]), 4)

      moran_p <- as.vector(x@gene_het[[gene]][[i]]$morans_I$p.value)
      geary_p <- as.vector(x@gene_het[[gene]][[i]]$gearys_C$p.value)
      getis_p <- as.vector(x@gene_het[[gene]][[i]]$getis_ord_Gi$p.value)

      if(moran_p < 0.05){
        moran_est <- paste0(moran_est, '*')
      }
      if(geary_p < 0.05){
        geary_est <- paste0(geary_est, '*')
      }
      if(getis_p < 0.05){
        getis_est <- paste0(getis_est, '*')
      }

      # Construct title.
      if(krige_type == 'ord'){
        titlekrige <- paste0(gene, ", subj ", i, " - ordinary kriging voom-norm counts\nMorans I=",
                             moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                             getis_est)
      } else if(krige_type == 'univ'){
        titlekrige <- paste0(gene, ", subj ", i, " - universal kriging voom-norm counts\nMorans I=",
                             moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                             getis_est)
      }

      if(purity){
        tumorstroma_df <- dplyr::bind_cols(x@coords[[i]],
                                           cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
        kp <- krige_p_purity(data_f=df, mask=bbox_mask_diff, color_pal=color_pal,
                             tumorstroma=tumorstroma_df,
                             leg_name="pred_expr", title_name=titlekrige)
      } else{
        kp <- krige_p(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_expr",
                      title_name=titlekrige)
      }
      # Append plot to list.
      kp_list[[gene]] <- kp
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
      pdf(file=paste0("gene_krige_spatarray_", i, ".pdf"))#,
         #width=w_pdf, height=h_pdf
      print(ggpubr::ggarrange(plotlist=kp_list,
                              nrow=row_col[1], ncol=row_col[2]))
      dev.off()
    } else{
      # Print plots to console.
       return(kp_list)
    }
  }
}
