##
#' @title gene_interpolation: Spatial interpolation of gene expression
#' @description Performs spatial interpolation ("kriging") of transformed gene counts
#' @details
#' This function takes an STlist and a vector of gene names and generates spatial
#' interpolation of gene expression values via "kriging". If genes='top', then
#' the 10 genes (default) with the highest standard deviation for each ST sample
#' are interpolated. The resulting interpolations can be visualized via the
#' `STplot_interpolation` function
#'
#' @param x an STlist with transformed RNA counts
#' @param genes a vector of gene names or 'top'. If 'top' (default), interpolation of
#' the 10 genes (`top_n` default) with highest standard deviation in each ST sample
#' is estimated.
#' @param top_n an integer indicating how many top genes to perform interpolation.
#' Default is 10.
#' @param samples the spatial samples for which interpolations will be performed.
#' If NULL (Default), all samples are interpolated.
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems.
#' @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
#' @return x a STlist including spatial interpolations.
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
#'   kp = STplot_interpolation(melanoma, genes=c('MLANA', 'COL1A1'))
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
gene_interpolation = function(x=NULL, genes='top', top_n=10, samples=NULL, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()
  if(verbose){
    cat('Gene interpolation started.\n')
  }

  #suppressPackageStartupMessages(require('fields'))
  #requireNamespace('fields')

  # Resolution removed after implementation of `fields` interpolation
  #res=0.5
  # Covariance structure between spots/cells
  # May later allow users to select Exponential or Mattern (both supported by fields)
  #covstr = 'Exponential'
  #covstr = 'Matern'
  #smoothness = 0.5
  # Specify number of spots in grid
  #nxy = ceiling(sqrt(ngrid))

  top_n = as.integer(top_n)
  #res = as.double(res)

  # Test that transformed counts are available
  if(rlang::is_empty(x@tr_counts)) {
    stop("There are not transformed counts in this STList.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(samples)) {
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
  }

  # Test that a gene name was entered.
  if(is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1 && genes == 'top'){
    genes = c()
    for(i in samples){
      # Find tops variable genes using Seurat approach
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

  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_krige[[gene]]) && rlang::is_empty(x@gene_krige[[gene]])){
      x@gene_krige[[gene]] = list()
      for(i in samples){
        x@gene_krige[[gene]][[names(x@tr_counts[i])]] = list()
      }
    }
  }

  # Save miscellaneous data
  x@misc[['gene_krige']] = list()
  # Store kriging type
  x@misc[['gene_krige']][['type']] = 'ordinary'

  # Specify resolution if not input by user
  # if(is.null(res)){
  #   res = 0.5
  # }
  #x@misc[['gene_krige']][['res']] = res

  # Give warning about kriging res < 0.5 for large matrices (e.g. Visium)
  # sizes = lapply(x@spatial_meta[samples], nrow)
  # if(res < 0.5 & (any(sizes > 1000))){
  #   cat('Kriging at the requested resolution is computationally expensive. Setting to res=0.5\n')
  #   res=0.5
  # }
  # rm(sizes) # Clean environment

  # Create lists to store prediction grids and borders.
  if(is.null(x@misc[['gene_krige']][['krige_border']])){
    x@misc[['gene_krige']][['krige_border']] = list()
    for(k in samples){
      x@misc[['gene_krige']][['krige_border']][[k]] = list()
    }
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in samples){
    subsetgenes_mask = genes %in% rownames(x@tr_counts[[i]])
    subsetgenes = genes[subsetgenes_mask]
    combo = dplyr::bind_rows(combo, expand.grid(i, subsetgenes))

    # Get genes not present.
    notgenes = genes[!subsetgenes_mask]

    if(!rlang::is_empty(notgenes)){
      message(paste(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i, ".\n"))
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@misc[['gene_krige']][['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@spatial_meta[[i]][, c('xpos', 'ypos')]))

    rm(subsetgenes_mask, subsetgenes, notgenes) # Clean env
  }

  # Save interpolation method
  #x@misc[['gene_krige']][['gene_krige_algorithm']] = 'fields'
  x@misc[['gene_krige']][['gene_krige_algorithm']] = 'gstat'

  # Store decompressed transformed counts list in separate object
  # tr_counts = list()
  # for(mtx in samples){
  #   tr_counts[[mtx]] = expandSparse(x@tr_counts[[mtx]])
  # }

  # Define cores available
  # Define cores available ### PARALLEL
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(unique(combo[[1]])))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Loop through combinations of samples
  kriging_list = parallel::mclapply(seq_along(1:length(as.vector(unique(combo[[1]])))), function(i_combo){
    i = as.vector(unique(combo[[1]]))[i_combo]
    genes_tmp = as.vector(unlist(combo[combo[[1]] == i, 2]))

    # Process genes
    kriging_res = list()
    for(j in genes_tmp){
      # Get transformed counts
      gene_expr = data.frame(
        position=names(unlist(x@tr_counts[[i]][rownames(x@tr_counts[[i]]) == j, ])),
        expr_val=as.vector(unlist(x@tr_counts[[i]][rownames(x@tr_counts[[i]]) == j, ])))

      # Match order of library names in counts data and coordinate data
      gene_expr = gene_expr[match(x@spatial_meta[[i]][[1]], gene_expr[[1]]), ]
      gene_geo_df = dplyr::bind_cols(x@spatial_meta[[i]][c('xpos', 'ypos')], gene_expr=as.numeric(gene_expr[[2]]))
      gene_geo_df = as.data.frame(gene_geo_df)
      sp::coordinates(gene_geo_df) = ~xpos+ypos

      rm(gene_expr) # Clean env

      # Create a grid with points in the middle of each coordinate
      ny = range(sp::coordinates(gene_geo_df)[, 'ypos'])
      ny = length(seq(ny[1], ny[2], by=100))
      nx = range(sp::coordinates(gene_geo_df)[, 'xpos'])
      nx = length(seq(nx[1], nx[2], by=100))
      grid_sf = sf::st_make_grid(gene_geo_df, n=c(nx, ny), what="centers")

      # Give hope to users...
      if(verbose > 1L){
        system(sprintf('echo "%s"', paste0("\tInterpolating ", j, ' in ', i, "....")))
      }

      # Perform kriging
      suppressMessages({
        krige_out = tryCatch({
          dat_vgm = gstat::variogram(gene_expr~1, gene_geo_df) # Calculates sample variogram values
          est_vgm = gstat::fit.variogram(dat_vgm, gstat::vgm('Exp')) # Fit variaogram
          krige_pred = gstat::krige(gene_expr~1, gene_geo_df, grid_sf, model=est_vgm) # Predict values
        },
        error=function(err){return(err)}
        )
      })

      # Compute surface if parameters available and save kriging values
      # Record if kriging was successful
      #if(any(class(cov_est) == 'simpleError') | !cov_est[['optimSuccess']]){
      if(any(class(krige_out) == 'simpleError')){
        kriging_res[[j]] = list(krige_out=NULL,
                                success_or_not='error')
      } else{
        krige_tmp = as.data.frame(sf::st_coordinates(krige_out))
        krige_tmp[['krige']] = krige_out[['var1.pred']]
        kriging_res[[j]] = list(krige_out=krige_tmp)
        rm(krige_tmp) # Clean env
        if(any(class(krige_out) == 'simpleWarning')){
          kriging_res[[j]][['success_or_not']] = 'warning'
        } else{
          kriging_res[[j]][['success_or_not']] = 'kriging_completed'
        }
      }
      rm(krige_out) # Clean env
    }

    return(kriging_res)
  }, mc.cores=cores, mc.preschedule=F)
  #names(kriging_list) = paste(combo[[1]], combo[[2]], sep='&&')
  names(kriging_list) = as.vector(unique(combo[[1]]))

  # Store kriging results in STList.
  for(i in names(kriging_list)){
    for(j in names(kriging_list[[i]])){
      x@gene_krige[[j]][[i]] = kriging_list[[i]][[j]]
    }
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('Gene interpolation completed in ', round(end_t, 2), ' min.\n'))
  }

  return(x)
}

