##
#' @title plot_deconv_krige: Visualize spatial interpolation of deconvolution scores
#' @description Produces a spatial interpolation surface from kriging of gene expression
#' deconvolution scores.
#' @details
#' This function produces a spatial interpolation surface plot for a series of
#' cell names and spatial arrays.
#'
#' @param x an STList with kriging objects for the cells selected.
#' @param cells, a vector of cell names (one or several) within a deconvolution
#' matrix  to plot. If 'top', the 10 cell types with highest standard deviation
#' from each spatial array ar plotted.
#' @param plot_who a vector of subject indexes as ordered within the STList, to
#' plot cells from. If NULL, will plot for all subjects.
#' @param color_pal a color scheme from 'khroma' or RColorBrewer.
#' @param saveplot, logical indicating whether or not save plots in a PDF file.
#' The PDFs are saved in the working directory. Default is FALSE, meaning plots
#' are printed to console.
#' @param purity, logical, whether or not to annotate tumor spots based on
#' ESTIMATE tumor purity scores.
#' @param pvalues, logical indicating whether or not dots where a given cell was
#' predicted to be significantly abundant (p<0.05). Default is FALSE.
#' @param visium, whether or not to reverse axes for Visium slides.
#' @return a list with plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # kplots <- plot_deconv_krige(melanoma,  cells=c('b_cells', 'i_dc'), plot_who=2, visium=F)
#'
#' @export
#
#
plot_deconv_krigeV2 <- function(x=NULL, cells=NULL, plot_who=NULL, color_pal='YlOrBr',
                              saveplot=F, pvalues=F, purity=F, visium=T){
  # Option to scale to 1 disabled.
  scaled=F
  # Option to use plotly disabled (not supporting geomPolypath)
  inter=F

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
  }

  # If cells='top', get names of 10 cells with the highest standard deviation.
  if(length(cells) == 1 && cells == 'top'){
    cells = c()
    for(i in plot_who){
      cells = append(cells, deconv_list[[i]]$cell_stdev$cell[order(deconv_list[[i]]$cell_stdev$cell_stdevs, decreasing=T)][1:10])
    }
    # Get unique genes from most variable.
    cells = unique(cells)
  }

  # Store maximum score value in case 'scaled' is required.
  # if(scaled){
  maxvalue <- c()
  minvalue <- c()
  for (i in plot_who) {
    for (cell in cells) {
      # Test if kriging exists for a cell and subject.
      if (rlang::has_name(x@cell_krige, cell)){
        if(length(x@cell_krige[[cell]]) >= i){
          # Find maximum expression value for each spatial array.
          values <- x@cell_krige[[cell]][[i]]
          maxvalue <- append(maxvalue, max(values))
          minvalue <- append(minvalue, min(values))
        }
      }
    }
  }
  # Find maximum value among selected spatial arrays.
  maxvalue <- max(maxvalue)
  minvalue <- min(minvalue)
  #  }

  # Create list of plots.
  kp_list <- list()

  # Loop through each of ythe subjects.
  for (i in plot_who) {
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
      predict_grid <- x@deconv_krige_data$krige_grid[[i]]

      # Create data frame with coordinates and kriging values.
      if(scaled){
        krige_vals <- (x@cell_krige[[cell]][[i]])/maxvalue
      } else{
        krige_vals <- x@cell_krige[[cell]][[i]]
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
      bbox_sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bbox)), "id")))

      # Create Spatial Polygon with the inner tissue border (concave hull)
      mask_sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(x@deconv_krige_data$krige_border[[i]])), "id")))

      # Substract the concave hull from the bounding box, yielding a SpatialPolygon object.
      bbox_mask_diff <- raster::erase(bbox_sp, mask_sp)

      # Get spatial statistics.
      # moran_est <- round(as.vector(x@cell_het[[cell]][[i]]$morans_I$estimate[[1]]), 2)
      # geary_est <- round(as.vector(x@cell_het[[cell]][[i]]$gearys_C$estimate[[1]]), 2)
      # getis_est <- round(as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$estimate[[1]]), 4)
      #
      # moran_p <- as.vector(x@cell_het[[cell]][[i]]$morans_I$p.value)
      # geary_p <- as.vector(x@cell_het[[cell]][[i]]$gearys_C$p.value)
      # getis_p <- as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$p.value)
      #
      # if(moran_p < 0.05){
      #   moran_est <- paste0(moran_est, '*')
      # }
      # if(geary_p < 0.05){
      #   geary_est <- paste0(geary_est, '*')
      # }
      # if(getis_p < 0.05){
      #   getis_est <- paste0(getis_est, '*')
      # }

      # Construct title.
      titlekrige <- paste0(cell, ", subj ", i, " - ", x@deconv_krige_data$krige_type, " kriging")
      # titlekrige <- paste0(cell, ", subj ", i, " - ", x@deconv_krige_data$krige_type, " kriging\nMorans I=",
      #                      moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=", getis_est)

      if(purity || pvalues){
        if(purity){
          tumorstroma_df <- dplyr::bind_cols(x@coords[[i]], cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
          kp <- krige_p_purity(data_f=df, mask=bbox_mask_diff, color_pal=color_pal,
                               tumorstroma=tumorstroma_df,
                               leg_name="pred_score", title_name=titlekrige, minvalue=minvalue,
                               maxvalue=maxvalue, visium=visium)
          if(!(rlang::is_empty(x@cell_deconv))){
            # Add a dummy duplicate column (expression values are added to this data frame when quilt plots)
            df_q = dplyr::bind_cols(x@coords[[i]][,-1],
                                    cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster,
                                    cluster2=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
            colnames(df_q) <- c('y_pos', 'x_pos', 'cluster', 'cluster2')
            qpbw <- quilt_p_purity_bw(data_f=df_q, visium=visium, title_name=paste0('ESTIMATE\ntumor/stroma - subj ', i))
          } else{
            stop("No tumor/stroma classification in the STList.")
          }

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
    if(purity){
      kp_list[[paste0('subj',i)]] <- qpbw
    }
  }
  row_col <- c(1, 1)

  # Test if plot should be saved to PDF.
  if(saveplot){
    # Print plots to PDF.
    pdf(file="cell_krige_spatarray.pdf")
    print(ggpubr::ggarrange(plotlist=kp_list, nrow=row_col[1], ncol=row_col[2]))
    dev.off()
  } else{
    # Convert ggplots to plotly plots.
    if(inter){
      for(p in 1:length(kp_list)){
        kp_list[[p]] = plotly::ggplotly(kp_list[[p]])
      }
      return(kp_list)
    }else{
      # Print plots to console.
      return(kp_list)
    }
    #  }
  }
}
