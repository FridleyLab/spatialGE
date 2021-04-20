##
# This function produces a kriging plot for a series of cell names and subjects.
#
# @param x, an STList with kriging objects for the cells selected.
# @param genes, a vector of cell names (one or several) to plot.
# @param plot_who, a vector of subject indexes as ordered within the STList, to
# plot cells from. If NULL, will plot for all subjects.
# @color_pal, a scheme from 'khroma'.
# @param saveplot, a file path where quilt plots will be saved. If NULL, plots
# are printed to console
#
#
plot_cell_krige <- function(x=NULL, cells=NULL, krige_type='ord', plot_who=NULL,
           color_pal='YlOrBr', saveplot=F, pvalues=F){

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # if(!is.null(saveplot)){
  #   if(dir.exists(paste0(saveplot))){
  #     unlink(paste0(saveplot), recursive=T)
  #   }
  # }

  # Loop through each of ythe subjects.
  for (i in plot_who) {

    # Create list of plots for a given subject.
    kp_list <- list()

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
      df <- dplyr::bind_cols(predict_grid,
                             krige=x@cell_krige[[cell]][[i]][[krige_type]]$predict)
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
      if(cell != 'stroma_score'){
        moran_est <- round(as.vector(x@cell_het[[cell]][[i]]$morans_I$estimate[[1]]), 2)
        geary_est <- round(as.vector(x@cell_het[[cell]][[i]]$gearys_C$estimate[[1]]), 2)
        getis_est <- round(as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$estimate[[1]]), 4)
      }

      # Construct title.
      if(cell != 'stroma_score'){
        if(krige_type == 'ord'){
          titlekrige <- paste0(cell, ", subj ", i, " - ordinary kriging square root-norm scores\nMorans I=",
                               moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                               getis_est)
        }
        else if(krige_type == 'univ'){
          titlekrige <- paste0(cell, ", subj ", i, " - universal kriging square root-norm scores\nMorans I=",
                               moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                               getis_est)
        }
      # The following titles for "stroma score". In the future, purity scores will
      # replace the spatial statistics.
      }else {
        if(krige_type == 'ord'){
          titlekrige <- paste0(cell, ", subj ", i, " - ordinary kriging square root-norm scores")
        }
        else if(krige_type == 'univ'){
          titlekrige <- paste0(cell, ", subj ", i, " - universal kriging square root-norm scores")
        }
      }

      if(pvalues){
        kp <- krige_p_pvals(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_score",
                            title_name=titlekrige, x=x, plot_who=i, cell=cell)
      }else{
        kp <- krige_p(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_score",
                      title_name=titlekrige)
      }

      # Append plot to list.
      kp_list[[cell]] <- kp
    }

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
      pdf(file=paste0("cell_krige_spat_array_", i, ".pdf"),
          width=w_pdf, height=h_pdf)
      print(ggpubr::ggarrange(plotlist=kp_list,
                          nrow=row_col[1], ncol=row_col[2]))
      dev.off()
    }else{
      # Print plots to console.
      print(ggpubr::ggarrange(plotlist=kp_list,
                          nrow=row_col[1], ncol=row_col[2]))
    }
  }
}
