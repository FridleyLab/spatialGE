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
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- transform_data(melanoma)
#'   melanoma <- gene_interpolation(melanoma, genes=c('MLANA', 'COL1A1'), samples='ST_mel1_rep2')
#'   kp = STplot_interpolation(melanoma, genes=c('MLANA', 'COL1A1'), samples='ST_mel1_rep2')
#'   ggpubr::ggarrange(plotlist=kp)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%
#
STplot_interpolation = function(x=NULL, genes=NULL, top_n=10, samples=NULL, color_pal='BuRd'){

  # To prevent NOTES in R CMD check
  . = NULL

  # Force number of top variable genes as integer
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
      x@gene_meta[[i]] = Seurat_FindVariableFeatures(x@counts[[i]]) %>%
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
        if(!is.null(x@gene_krige[[gene]][[i]][['success_or_not']])){
          if(x@gene_krige[[gene]][[i]][['success_or_not']] != 'error'){
            # Find maximum expression value for each spatial array.
            values <- x@gene_krige[[gene]][[i]][['krige_out']][['krige']]
            maxvalue <- append(maxvalue, max(values, na.rm=T))
            minvalue <- append(minvalue, min(values, na.rm=T))
          }
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
        message(paste0("Kriging for subject ", gene, " in sample ", i, " is not present in STlist\n"))
        next
      }

      # Create data frame with coordinates and kriging values.
      krige_vals = x@gene_krige[[gene]][[i]][['krige_out']]
      #krige_vals[['krige']] = x@gene_krige[[gene]][[i]][['krige_out']][['var1.pred']]
#      colnames(krige_vals) = paste0('id_', x@gene_krige[[gene]][[i]][['krige_out']][['y']])
#      rownames(krige_vals) = paste0('id_', x@gene_krige[[gene]][[i]][['krige_out']][['x']])
#      suppressWarnings({krige_vals = reshape::melt.array(krige_vals)})

      # Looks like reshape::melt.array pass columns from matrix to first column, and
      # columns from matrix to second column. Values go in third column
      names(krige_vals) <- c("x_pos", "y_pos", "krige")
      # krige_vals = krige_vals %>%
      #   dplyr::mutate(y_pos=gsub('id_', '', y_pos) %>% as.numeric()) %>%
      #   dplyr::mutate(x_pos=gsub('id_', '', x_pos) %>% as.numeric())

      # Get coordinates of bounding box enclosing the predicted grid.
      bbox = rbind(
        c(min(krige_vals$x_pos)-1, min(krige_vals$y_pos)-1),
        c(min(krige_vals$x_pos)-1, max(krige_vals$y_pos)+1),
        c(max(krige_vals$x_pos)+1, max(krige_vals$y_pos)+1),
        c(max(krige_vals$x_pos)+1, min(krige_vals$y_pos)-1),
        c(min(krige_vals$x_pos)-1, min(krige_vals$y_pos)-1)
      )
      bbox <- as.data.frame(bbox)
      colnames(bbox) <- c("x", "y")

      # Extract coordinates from kriging
      mask_df = as.data.frame(x@misc[['gene_krige']][['krige_border']][[i]][, c('V1', 'V2')])

      visium = F
      # Define if visium
      if(x@misc[['platform']] == 'visium'){
        #bbox[['x']] = rev(bbox[['x']])
        bbox[['y']] = bbox[['y']] *-1
        #mask_df[['V1']] = rev(mask_df[['V1']])
        mask_df[['V2']] = mask_df[['V2']] *-1

        krige_vals[['y_pos']] = krige_vals[['y_pos']] *-1
      }

      # Create SpatialPolygon object with the bounding box.
      bbox_sp = sf::st_as_sf(bbox, coords=c('x', 'y')) %>%
        sf::st_combine() %>%
        sf::st_cast(to='POLYGON')
        # dplyr::group_by(id) %>%
        # dplyr::summarise(geometry=sf::st_combine(geometry))
      #bbox_sp =  %>% sf::st_sf()
      #sf::st_agr(bbox_sp) = 'constant'

      # Create Spatial Polygon with the inner tissue border (concave hull)
      mask_sp = sf::st_as_sf(mask_df, coords=c('V1', 'V2')) %>%
        sf::st_combine() %>%
        sf::st_cast(to='POLYGON')

      # Substract the concave hull from the bounding box, yielding a SpatialPolygon object.
      bbox_mask_diff <- sf::st_difference(bbox_sp, mask_sp)

      # Construct title.
      titlekrige <- paste0(gene, " (interpolation)\nsample: ", i)

      kp <- krige_p(data_f=krige_vals, mask=bbox_mask_diff, color_pal=color_pal, leg_name="pred_expr",
                    title_name=titlekrige, minvalue=minvalue, maxvalue=maxvalue, visium=visium)

      # Append plot to list.
      kp_list[[paste0(gene, "_", i)]] <- kp
    }
  }
  return(kp_list)
}

