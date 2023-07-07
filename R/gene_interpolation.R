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
#' @param ngrid an integer indicating the number of point to predict. Default is 10000,
#' resulting in a grid of 100 x 100 points. Larger numbers provide more resolution,
#' but computing time is longer.
#' @return x a STlist including spatial interpolations.
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
#' @import fields
#' @importFrom magrittr %>%
#'
gene_interpolation = function(x=NULL, genes='top', top_n=10, samples=NULL, ngrid=10000){
  #require("magrittr")
  suppressPackageStartupMessages(require('fields'))

  # Universal kriging not implemented for now
  univ=F
  # Python kriging disabled
  python=F
  # Resolution removed after implementation of `fields` interpolation
  res=0.5
  # Covariance structure between spots/cells
  # May later allow users to select Exponential or Mattern (both supported by fields)
  covstr = 'Exponential'
  #covstr = 'Matern'
  #smoothness = 0.5
  # Specify number of spots in grid
  nxy = round(sqrt(ngrid), 0)

  top_n = as.integer(top_n)
  res = as.double(res)

  # geoR implementation of universal kriging to be implemented.
  # Probably will allow users to specify parameters from variogram
  if(python == F && univ == T){
    stop('Currently, universal kriging is only available using Python kriging (PyKrige)')
  }

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

  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_krige[[gene]]) && rlang::is_empty(x@gene_krige[[gene]])){
      x@gene_krige[[gene]] = list()
      for(i in samples){
        x@gene_krige[[gene]][[names(x@tr_counts[i])]] = list()
      }
    }
  }

  # Store kriging type.
  if(univ){
    x@misc[['gene_krige_type']] = 'universal'
  }else{
    x@misc[['gene_krige_type']] = 'ordinary'
  }

  # Specify resolution if not input by user
  if(is.null(res)){
    res = 0.5
  }
  x@misc[['gene_krige_res']] = res

  # Give warning about kriging res < 0.5 for large matrices (e.g. Visium)
  sizes = lapply(x@spatial_meta[samples], nrow)
  if(res < 0.5 & (any(sizes > 1000))){
    cat('Kriging at the requested resolution might take some time. Setting to res=0.5\n')
    res=0.5
  }
  rm(sizes) # Clean environment

  # Create lists to store prediction grids and borders.
  if(is.null(x@misc[['krige_border']])){
    x@misc[['krige_border']] = list()
    #x@misc[['gene_krige_grid']] = list()
    for(k in samples){
      x@misc[['krige_border']][[k]] = list()
      #x@misc[['gene_krige_grid']][[k]] = list()
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
      cat(paste(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i, ".\n"))
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@misc[['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@spatial_meta[[i]][, c('xpos', 'ypos')]))

    # Create grid for PyKrige or geoR.
    # if(python == T){
    #   gridx = seq(min(x@spatial_meta[[i]][['xpos']]), max(x@spatial_meta[[i]][['xpos']]), res)
    #   gridy = seq(min(x@spatial_meta[[i]][['ypos']]), max(x@spatial_meta[[i]][['ypos']]), res)
    #   gridx = gridx[-length(gridx)]
    #   gridy = gridy[-length(gridy)]
    #   gene_geo_grid = expand.grid(
    #     seq(min(gridx), max(gridx), by=res),
    #     seq(min(gridy), max(gridy), by=res)
    #   )
    #   x@misc[['gene_krige_grid']][[i]] = gene_geo_grid
    # } else if(python == F && univ == F){ #  Create grid for geoR estimation.
    # gene_geo_grid <-expand.grid(
    #   seq((min(x@spatial_meta[[i]][['xpos']])-1), (max(x@spatial_meta[[i]][['xpos']])+1), by=res),
    #   seq((min(x@spatial_meta[[i]][['ypos']])-1), (max(x@spatial_meta[[i]][['ypos']])+1), by=res)
    # )
    # gene_geo_grid = list(xpos=gene_geo_grid[[1]], ypos=gene_geo_grid[[2]])
    #   x@misc[['gene_krige_grid']][[i]] = gene_geo_grid
    # }

    rm(subsetgenes_mask, subsetgenes, notgenes) # Clean env
  }

  # Store prediction grid and kriging algorithm in STList.
  if(python){
    x@misc[['gene_krige_algorithm']] = 'pykrige'
  } else{
    x@misc[['gene_krige_algorithm']] = 'fields'
  }

  # Store decompressed transformed counts list in separate object
  tr_counts = list()
  for(mtx in samples){
    tr_counts[[mtx]] = expandSparse(x@tr_counts[[mtx]])
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  kriging_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get transformed counts.
    gene_expr = tr_counts[[i]][rownames(tr_counts[[i]]) == j, ]
    gene_expr = as.data.frame(t(gene_expr))
    gene_expr = gene_expr %>%
      tibble::rownames_to_column(., var="position")
    # Match order of library names in counts data and coordinate data.
    gene_expr = gene_expr[match(x@spatial_meta[[i]][[1]], gene_expr[[1]]), ]
    gene_geo_df = dplyr::bind_cols(x@spatial_meta[[i]][c('xpos', 'ypos')], gene_expr=as.numeric(gene_expr[[2]]))

    # Call the requested Kriging algorithm
    if(python == T){
      # Call PyKrige implementation.
      gridx = seq(min(x@spatial_meta[[i]][[3]]), max(x@spatial_meta[[i]][[3]]), res)
      gridy = seq(min(x@spatial_meta[[i]][[2]]), max(x@spatial_meta[[i]][[2]]), res)
      gridx = gridx[-length(gridx)]
      gridy = gridy[-length(gridy)]
      kriging_res = krige_py(gridx=gridx, gridy=gridy, geo_df=gene_geo_df, univ=univ)
    } else{
      # Give hope to users...
      system(sprintf('echo "%s"', crayon::yellow(paste0("\tInterpolating ", j, ' in ', i, "...."))))
      # Fit spatial process
      suppressMessages({
        cov_est = fields::spatialProcess(x=gene_geo_df[, c('ypos', 'xpos')], y=gene_geo_df[['gene_expr']],
                                         cov.args=list(Covariance=covstr#,
                                                       #smoothness=smoothness
                                                       ), verbose=F, REML=F, reltol=1e-3)
      })
      # Compute surface
      kriging_res = fields::predictSurface(cov_est, nx=nxy, ny=nxy, extrap=T#,
                                            #grid.list=gene_geo_grid
      )
    }
    return(kriging_res)
  }, mc.cores=cores, mc.preschedule=F)
  names(kriging_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(kriging_list)[i], split = '&&'))
    x@gene_krige[[combo_name[2]]][[combo_name[1]]] = kriging_list[[i]]
  }

  return(x)
}

