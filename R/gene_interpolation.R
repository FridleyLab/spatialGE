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
#'   zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
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

  top_n = as.integer(top_n)

  # Test that transformed counts are available
  assert_transformed(x, message="There are not transformed counts in this STList.")

  # Test if no specific subject plot was requested.
  samples = resolve_samples(x, samples, slot='tr_counts', numeric_handling='plain',
                            char_transform='none', post_check='none')

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
      x@gene_meta[[i]] = Seurat_FindVariableFeatures(sample_counts(x, i)) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::select('gene', 'vst.variance.standardized') %>%
        dplyr::left_join(x@gene_meta[[i]], ., by='gene')

      genes = append(genes, x@gene_meta[[i]][['gene']][order(x@gene_meta[[i]][['vst.variance.standardized']], decreasing=T)][1:top_n])
    }
    # Get unique genes from most variable.
    genes = unique(genes)
  }

  # Drop requested genes that are absent, and fail clearly if none survive.
  # Previously an entirely absent gene set built empty kriging entries here and
  # then failed much later inside the gstat formula with "object 'xpos' not
  # found" -- an error that points at neither the gene nor this function.
  present = Reduce(union, lapply(samples, function(s) sample_genes(x, s, 'tr_counts')))
  missing_genes = setdiff(genes, present)
  if(length(missing_genes) > 0){
    message(paste0(paste(missing_genes, collapse=', '),
                    ': Not present in the transformed counts. Skipping.\n'))
  }
  genes = intersect(genes, present)
  if(length(genes) == 0){
    raise_err(err_code='error0046')
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
    subsetgenes_mask = genes %in% sample_genes(x, i, 'tr_counts')
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
  cores = resolve_cores(length(unique(combo[[1]])), cores, coerce='as_integer')
  if(is.na(cores)){
    stop('Could not recognize number of cores requested')
  }

  # Loop through combinations of samples
  kriging_list = sge_lapply(seq_along(1:length(as.vector(unique(combo[[1]])))), function(i_combo){
    i = as.vector(unique(combo[[1]]))[i_combo]
    genes_tmp = as.vector(unlist(combo[combo[[1]] == i, 2]))
    interp_one_sample(sample_tr_counts(x, i), sample_spatial_meta(x, i), genes_tmp, i, verbose)
  }, cores=cores, preschedule=F)
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

