##
#' @title plot_deconv_krige: Visualize spatial interpolation of deconvolution scores
#' @description Produces a spatial interpolation surface from kriging of gene expression
#' deconvolution scores in ST data.
#' @details
#' This function produces a spatial interpolation surface plot for a series of
#' cell names and spatial arrays.
#'
#' @param x, an STList with kriging objects for the cells selected.
#' @param cells, a vector of cell names (one or several) within a deconvolution
#' matrix  to plot.
#' @param krige_type, either 'ord' (ordinary; default), or 'univ' (universal)
#' kriging. Data for the respective kriging must be generated previously with
#' cell_krige().
#' @param plot_who, a vector of subject indexes as ordered within the STList, to
#' plot cells from. If NULL, will plot for all subjects.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param purity, logical, whether or not annotate tumor and stroma spots based on
#' ESTIMATE tumor purity scores.
#' @param saveplot, logical indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @param pvalues, logical indicating whether or not dots where a given cell was
#' predicted to be significantly abundant (p<0.05). Default is FALSE.
#' @param scaled, logical, indicating if expression values should be scaled with
#' respect to the highest value among all genes to plot. WARNING: Color legends
#' are not scaled between plots, but values are.
#' @export
#
#
plot_deconv_krige <- function(x=NULL, cells=NULL, krige_type='ord', plot_who=NULL,
           color_pal='YlOrBr', saveplot=F, pvalues=F, purity=F, scaled=T, visium=T){

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # Store maximum expression value in case 'scaled' is required.
#  if(scaled){
    maxvalue <- c()
    minvalue <- c()
    for (i in plot_who) {
      for (cell in cells) {
        # Test if kriging exists for a cell and subject.
        if (rlang::has_name(x@cell_krige, cell)){
          if(length(x@cell_krige[[cell]]) >= i){
            if (rlang::has_name(x@cell_krige[[cell]][[i]], krige_type)){
              # Find maximum expression value for each spatial array.
              values <- x@cell_krige[[cell]][[i]][[krige_type]]$predict
              maxvalue <- append(maxvalue, max(values))
              minvalue <- append(minvalue, min(values))
            }
          }
        }
        # Find maximum value among selected spatial arrays.

      }
    }

    maxvalue <- max(maxvalue)
    minvalue <- min(minvalue)
#  }

  # Create list of plots.
  kp_list <- list()

  # Loop through each of ythe subjects.
  for (i in plot_who) {

    # Create list of plots for a given subject.
    # kp_list <- list()

    # Loop though cells to plot.
    for (cell in cells) {

      if(length(x@cell_krige[[cell]]) >= i){
        if(is.null(x@cell_krige[[cell]][[i]])){
          cat(paste0(cell, " kriging for subject ", i, " is not present in STList\n"))
          next
        }
      }else{
        cat(paste0(cell, " kriging for subject ", i, " is not present in STList\n"))
        next
      }

      # Find prediction grid.
      if(krige_type == 'univ'){
        predict_grid <- x@cell_krige[[cell]][[i]]$univ_grid
      }else if(krige_type == 'ord'){
        predict_grid <- x@cell_krige[[cell]][[i]]$ord_grid
      }

      # Create data frame with coordinates and kriging values.
      if(scaled){
        krige_vals <- (x@cell_krige[[cell]][[i]][[krige_type]]$predict)/maxvalue
      } else{
        krige_vals <- x@cell_krige[[cell]][[i]][[krige_type]]$predict
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
      moran_est <- round(as.vector(x@cell_het[[cell]][[i]]$morans_I$estimate[[1]]), 2)
      geary_est <- round(as.vector(x@cell_het[[cell]][[i]]$gearys_C$estimate[[1]]), 2)
      getis_est <- round(as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$estimate[[1]]), 4)

      moran_p <- as.vector(x@cell_het[[cell]][[i]]$morans_I$p.value)
      geary_p <- as.vector(x@cell_het[[cell]][[i]]$gearys_C$p.value)
      getis_p <- as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$p.value)

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
        titlekrige <- paste0(cell, ", subj ", i, " - ordinary kriging\nMorans I=",
                             moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                             getis_est)
      } else if(krige_type == 'univ'){
        titlekrige <- paste0(cell, ", subj ", i, " - universal kriging\nMorans I=",
                             moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                             getis_est)
      }

      if(purity | pvalues){
        if(purity){
          tumorstroma_df <- dplyr::bind_cols(x@coords[[i]],
                                             cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
          kp <- krige_p_purity(data_f=df, mask=bbox_mask_diff, color_pal=color_pal,
                               tumorstroma=tumorstroma_df,
                               leg_name="pred_score", title_name=titlekrige,
                               minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        } else if(pvalues){
          kp <- krige_p_pvals(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_score",
                              title_name=titlekrige, x=x, plot_who=i, cell=cell,
                              minvalue=minvalue, maxvalue=maxvalue, visium=visium)
        }

      } else{
        kp <- krige_p(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_score",
                      title_name=titlekrige, minvalue=minvalue, maxvalue=maxvalue,
                      visium=visium)
      }
      # Append plot to list.
      kp_list[[paste0(cell, "_", i)]] <- kp
    }

    #     # Define number of columns and rows in plot and size.
    #     row_col <- c(2, 2)
    # #    w_pdf=9
    # #    h_pdf=9
    #     if(length(cells) == 2){
    #       row_col <- c(1, 2)
    # #      w_pdf=9
    # #      h_pdf=4.5
    #     }else if(length(cells) == 1){
    #       row_col <- c(1, 1)
    # #      w_pdf=7
    # #      h_pdf=7
    #     }
}
    row_col <- c(1, 1)

    # Test if plot should be saved to PDF.
    if(saveplot){
      # Print plots to PDF.
      pdf(file="cell_krige_spatarray.pdf")#,
      print(ggpubr::ggarrange(plotlist=kp_list,
                             nrow=row_col[1], ncol=row_col[2]))
      dev.off()
    }else{
      # Print plots to console.
      return(kp_list)
    }
#  }
}
