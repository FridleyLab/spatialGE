##
#' @title plot_gene_krige: Visualize transcriptomic surfaces
#' @description Produces a transcriptomic surface from kriging interpolation of ST data.
#' @details
#' This function produces a transcriptomic surface plot for a series of HUGO gene names and
#' spatial arrays.
#'
#' @param x an STList with kriging objects for the genes selected.
#' @param genes a vector of gene names (one or several) to plot.  If 'top', the 10
#' genes with highest standard deviation from each spatial array ar plotted.
#' @param plot_who, a vector of numbers indicating the spatial arrays to plot
#' genes from. Numbers follow the order in `names(x@counts)`. If NULL, will plot
#' all spatial arrays.
#' @param color_pal a color scheme from 'khroma' or RColorBrewer.
#' @param purity logical, whether or not annotate tumor spots based on
#' ESTIMATE tumor purity scores.
#' @param image logical, whether to print the image stored for the spatial arrays
#' @param saveplot, a file name specifying the name of a PDF file to write plots to.
#' @param visium whether or not to reverse axes for Visium slides.
#' @return a list with plots.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # kplots <- plot_gene_krige(melanoma, genes=c('CD74', 'SOX10'), plot_who=3, visium=F)
#'
#' @export
#
#
plot_gene_krige = function(x=NULL, genes=NULL, plot_who=NULL, color_pal='YlOrBr', purity=F, image=F, saveplot=NULL, visium=T, ptsize=0.5){

  # Option to use plotly disabled (not supporting geomPolypath)
  inter=F

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(plot_who)) {
    plot_who <- c(1:length(x@counts))
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

  # Store maximum and minimum expression value for plot color scaling
  maxvalue <- c()
  minvalue <- c()
  for (i in plot_who) {
    for (gene in genes) {
      # Test if kriging exists for a gene and subject.
      if (rlang::has_name(x@gene_krige, gene)){
        if(length(x@gene_krige[[gene]]) >= i){
          # Find maximum expression value for each spatial array.
          values <- x@gene_krige[[gene]][[i]]
          maxvalue <- append(maxvalue, max(values))
          minvalue <- append(minvalue, min(values))
        }
      }
    }
  }
  # Find maximum value among selected spatial arrays.
  maxvalue <- max(maxvalue)
  minvalue <- min(minvalue)

  # Create list of plots.
  kp_list <- list()

  # Loop through each of the subjects.
  for (i in plot_who) {
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
      predict_grid = x@misc[['gene_krige_grid']][[i]]

      # Create data frame with coordinates and kriging values.
      krige_vals = x@gene_krige[[gene]][[i]]

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
      mask_sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(x@misc[['krige_border']][[i]])), "id")))

      # Substract the concave hull from the bounding box, yielding a SpatialPolygon object.
      bbox_mask_diff <- raster::erase(bbox_sp, mask_sp)

      # Construct title.
      titlekrige <- paste0(gene, " (kriging)\nsample: ", names(x@tr_counts[i]))

      if(purity){
        tumorstroma_df <- dplyr::bind_cols(x@coords[[i]], cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
        kp <- krige_p_purity(data_f=df, mask=bbox_mask_diff, color_pal=color_pal,
                             tumorstroma=tumorstroma_df,
                             leg_name="pred_expr", title_name=titlekrige, minvalue=minvalue,
                             maxvalue=maxvalue, visium=visium)
        if(!(rlang::is_empty(x@cell_deconv))){
          # Add a dummy duplicate column (expression values are added to this data frame when quilt plots)
          df_q = dplyr::bind_cols(x@coords[[i]][,-1],
                                  cluster=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster,
                                  cluster2=x@cell_deconv$ESTIMATE[[i]]$purity_clusters$cluster)
          colnames(df_q) <- c('y_pos', 'x_pos', 'cluster', 'cluster2')
          qpbw <- quilt_p_purity_bw(data_f=df_q, visium=visium,
                                    title_name=paste0('ESTIMATE\ntumor/stroma\nsample: ', names(x@tr_counts[i])),
                                    ptsize=ptsize)
        } else{
          stop("No tumor/stroma classification in the STList.")
        }

      } else{
        kp <- krige_p(data_f=df, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_expr",
                      title_name=titlekrige, minvalue=minvalue, maxvalue=maxvalue, visium=visium)
      }
      # Append plot to list.
      kp_list[[paste0(gene, "_", i)]] <- kp
    }

    if(image && !is.null(x@misc[['sp_images']][[i]])){
      img_obj = grid::rasterGrob(x@misc[['sp_images']][[i]])
      kp_list[[paste0('image', names(counts[i]))]] = ggplot() +
        annotation_custom(img_obj)
    }

    if(purity){
      kp_list[[paste0('subj',i)]] <- qpbw
    }
  }

  row_col <- c(1, 1)

  # Test if plot should be saved to PDF.
  if(!is.null(saveplot)){
    # Print plots to PDF.
    pdf(file=saveplot)
    print(ggpubr::ggarrange(plotlist=kp_list, nrow=row_col[1], ncol=row_col[2], common.legend=F, legend='bottom'))
    dev.off()
  } else{
    # Convert ggplots to plotly plots.
    if(inter == T && purity == T && image == F){
      for(p in 1:length(kp_list)){
        kp_list[[p]] = plotly::ggplotly(kp_list[[p]])
      }
      return(kp_list)
    }else{
      # Print plots to console.
      return(kp_list)
    }
  }
}

