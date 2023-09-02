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
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems.
#' @return x a STlist including spatial interpolations.
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
#' count_files <- grep("counts", data_files, value=T)
#' coord_files <- grep("mapping", data_files, value=T)
#' clin_file <- grep("thrane_clinical", data_files, value=T)
#' melanoma <- STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples
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
gene_interpolation = function(x=NULL, genes='top', top_n=10, samples=NULL, ngrid=10000, cores=NULL){
  # Record time
  zero_t = Sys.time()
  verbose = 1L
  if(verbose > 0L){
    cat(crayon::green(paste0('Gene interpolation started.\n')))
  }

  suppressPackageStartupMessages(require('fields'))

  # Resolution removed after implementation of `fields` interpolation
  res=0.5
  # Covariance structure between spots/cells
  # May later allow users to select Exponential or Mattern (both supported by fields)
  covstr = 'Exponential'
  #covstr = 'Matern'
  #smoothness = 0.5
  # Specify number of spots in grid
  nxy = ceiling(sqrt(ngrid))

  top_n = as.integer(top_n)
  res = as.double(res)

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

  # Save miscellaneous data
  x@misc[['gene_krige']] = list()
  # Store kriging type
  x@misc[['gene_krige']][['type']] = 'ordinary'

  # Specify resolution if not input by user
  if(is.null(res)){
    res = 0.5
  }
  x@misc[['gene_krige']][['res']] = res

  # Give warning about kriging res < 0.5 for large matrices (e.g. Visium)
  sizes = lapply(x@spatial_meta[samples], nrow)
  if(res < 0.5 & (any(sizes > 1000))){
    cat('Kriging at the requested resolution is computationally expensive. Setting to res=0.5\n')
    res=0.5
  }
  rm(sizes) # Clean environment

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
      cat(paste(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i, ".\n"))
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@misc[['gene_krige']][['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@spatial_meta[[i]][, c('xpos', 'ypos')]))

    rm(subsetgenes_mask, subsetgenes, notgenes) # Clean env
  }

  # Save interpolation method
  x@misc[['gene_krige']][['gene_krige_algorithm']] = 'fields'

  # Store decompressed transformed counts list in separate object
  # tr_counts = list()
  # for(mtx in samples){
  #   tr_counts[[mtx]] = expandSparse(x@tr_counts[[mtx]])
  # }

  # Define cores available
  # Define cores available ### PARALLEL
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

      rm(gene_expr) # Clean env

      # Give hope to users...
      system(sprintf('echo "%s"', crayon::yellow(paste0("\tInterpolating ", j, ' in ', i, "...."))))
      # Fit spatial process
      suppressMessages({
        cov_est = tryCatch({
          fields::spatialProcess(x=gene_geo_df[, c('ypos', 'xpos')], y=gene_geo_df[['gene_expr']],
                                 cov.args=list(Covariance=covstr #,smoothness=smoothness
                                 ), verbose=F, REML=T#, reltol=1e-3
          )
        },
        #warning=function(warn){return(warn)},
        error=function(err){return(err)}
        )
      })

      # Compute surface if parameters available and save kriging values
      # Record if kriging was successful
      if(any(class(cov_est) == 'simpleError') | !cov_est[['optimSuccess']]){
        kriging_res[[j]] = list(krige_out=NULL,
                                success_or_not='error')
      } else{
        kriging_res[[j]] = list(krige_out=fields::predictSurface(cov_est, nx=nxy, ny=nxy, extrap=T))
        if(any(class(cov_est) == 'simpleWarning')){
          kriging_res[[j]][['success_or_not']] = 'warning'
        } else{
          kriging_res[[j]][['success_or_not']] = 'kriging_completed'
        }
      }
      rm(cov_est) # Clean env
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
  if(verbose > 0L){
    cat(crayon::green(paste0('Gene interpolation completed in ', round(end_t, 2), ' min.\n')))
  }

  return(x)
}

