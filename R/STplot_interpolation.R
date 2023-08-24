##
#' @title STplot_interpolation: Visualize gene expression surfaces
#' @description Produces a gene expression surface from kriging interpolation of ST data.
#' @details
#' This function produces a gene expression surface plot via kriging for one or several
#' genes and spatial samples
#'
#' @param x an STlist containing results from `gene_krige` for the genes selected.
#' @param genes a vector of gene names (one or several) to plot. If 'top', the 10
#' genes with highest standard deviation from each spatial sample are plotted.
#' @param top_n an integer indicating how many top genes to perform kriging. Default is 10.
#' @param samples a vector indicating the spatial samples to plot. If vector of numbers,
#' it follows the order of `names(x@counts)`. If NULL, the function plots all samples
#' @param color_pal a color scheme from `khroma` or `RColorBrewer`.
#' @return a list of plots
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
#' count_files <- grep("counts", data_files, value=T)
#' coord_files <- grep("mapping", data_files, value=T)
#' clin_file <- grep("thrane_clinical", data_files, value=T)
#' melanoma <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#' melanoma <- transform_data(melanoma)
#' melanoma <- gene_interpolation(melanoma, genes=c('MLANA', 'COL1A1'), samples='ST_mel1_rep2')
#' kp = STplot_interpolation(melanoma, genes=c('MLANA', 'COL1A1'))
#' ggpubr::ggarrange(plotlist=kp)
#'
#' @export
#'
#' @importFrom magrittr %>%
#
#
STplot_interpolation = function(x=NULL, genes=NULL, top_n=10, samples=NULL, color_pal='BuRd'){

  top_n = as.integer(top_n)

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(samples)) {
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1 && genes == 'top'){
    genes = c()
    for(i in samples){
      # Find tops variable genes using Seurat approach. In the past, instead of Seurat, genes with the highest stdev were used
      if(any(colnames(x@gene_meta[[i]]) == 'vst.variance.standardized')){
        x@gene_meta[[i]] = x@gene_meta[[i]][, !grepl('vst.variance.standardized', colnames(x@gene_meta[[i]]))]
      }
      x@gene_meta[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::select('gene', 'vst.variance.standardized') %>%
        dplyr::left_join(x@gene_meta[[i]], ., by='gene')

      genes = append(genes, x@gene_meta[[i]][['gene']][order(x@gene_meta[[i]][['vst.variance.standardized']], decreasing=T)][1:top_n])
    }
    # Get unique genes from most variable.
    genes = unique(genes)
  }

  # Store maximum and minimum expression value for plot color scaling
  maxvalue <- c()
  minvalue <- c()
  for(i in samples) {
    for(gene in genes) {
      # Test if kriging exists for a gene and subject.
      if (rlang::has_name(x@gene_krige, gene)){
        if(x@gene_krige[[gene]][[i]][['success_or_not']] != 'error'){
        # Find maximum expression value for each spatial array.
        values <- x@gene_krige[[gene]][[i]][['krige_out']][['z']]
        maxvalue <- append(maxvalue, max(values, na.rm=T))
        minvalue <- append(minvalue, min(values, na.rm=T))
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
  for (i in samples) {
    # Loop though genes to plot.
    for (gene in genes){

      if(is.null(x@gene_krige[[gene]][[i]][['krige_out']])){
        cat(paste0("Kriging for subject ", gene, " in sample ", i, " is not present in STlist\n"))
        next
      }

      # Create data frame with coordinates and kriging values.
      krige_vals = x@gene_krige[[gene]][[i]][['krige_out']][['z']]
      colnames(krige_vals) = paste0('id_', x@gene_krige[[gene]][[i]][['krige_out']][['y']])
      rownames(krige_vals) = paste0('id_', x@gene_krige[[gene]][[i]][['krige_out']][['x']])
      suppressWarnings({krige_vals = reshape::melt.array(krige_vals)})

      # Looks like reshape::melt.array pass columns from matrix to first column, and
      # columns from matrix to second column. Values go in third column
      names(krige_vals) <- c("y_pos", "x_pos", "krige")
      krige_vals = krige_vals %>%
        dplyr::mutate(y_pos=gsub('id_', '', y_pos) %>% as.numeric()) %>%
        dplyr::mutate(x_pos=gsub('id_', '', x_pos) %>% as.numeric())

      # Get coordinates of bounding box enclosing the predicted grid.
      bbox = rbind(
        c(min(krige_vals$x_pos)-1, min(krige_vals$y_pos)-1),
        c(min(krige_vals$x_pos)-1, max(krige_vals$y_pos)+1),
        c(max(krige_vals$x_pos)+1, max(krige_vals$y_pos)+1),
        c(max(krige_vals$x_pos)+1, min(krige_vals$y_pos)-1),
        c(min(krige_vals$x_pos)-1, min(krige_vals$y_pos)-1)
      )
      bbox <- as.data.frame(bbox)
      names(bbox) <- c("x", "y")

      # Create SpatialPolygon object with the bounding box.
      bbox_sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bbox)), "id")))

      # Create Spatial Polygon with the inner tissue border (concave hull)
      mask_sp <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(x@misc[['gene_krige']][['krige_border']][[i]][, c('V1', 'V2')])), "id")))

      # Substract the concave hull from the bounding box, yielding a SpatialPolygon object.
      bbox_mask_diff <- raster::erase(bbox_sp, mask_sp)

      # Construct title.
      titlekrige <- paste0(gene, " (interpolation)\nsample: ", i)

      visium = F
      # Define if visium
      if(x@misc[['platform']] == 'visium'){
        visium=T
      }

      kp <- krige_p(data_f=krige_vals, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_expr",
                    title_name=titlekrige, minvalue=minvalue, maxvalue=maxvalue, visium=visium)

      # Append plot to list.
      kp_list[[paste0(gene, "_", i)]] <- kp
    }
  }
  return(kp_list)
}

