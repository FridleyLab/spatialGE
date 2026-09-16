##
#' @title SThet: Computes global spatial autocorrelation statistics on gene expression
#' @description Computes the global spatial autocorrelation statistics Moran's I and/or
#' Geary's C for a set of genes
#' @details The function computes global spatial autocorrelation statistics (Moran's I and/or
#' Geary's C) for the requested genes and samples. Then computation uses the
#' package `spdep`. The calculated statistics are stored in the STlist, which can
#' be accessed with the `get_gene_meta` function. For visual comparative analysis,
#' the function `compare_SThet` can be used afterwards.
#'
#' @param x an STlist
#' @param genes a vector of gene names to compute statistics
#' @param samples the samples to compute statistics
#' @param method The spatial statistic(s) to estimate. It can be set to 'moran',
#' 'geary' or both. Default is 'moran'
#' @param k the number of neighbors to estimate weights. By default NULL, meaning that
#' spatial weights will be estimated from Euclidean distances. If an positive integer is
#' entered, then the faster k nearest-neighbors approach is used. Please keep in mind
#' that estimates are not as accurate as when using the default distance-based method.
#' @param overwrite logical indicating if previous statistics should be overwritten.
#' Default to FALSE (do not overwrite)
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems
#' @param verbose logical, whether to print text to console
#' @return an STlist containing spatial statistics
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
#'   melanoma <- SThet(melanoma, genes=c('MLANA', 'TP53'), method='moran')
#'   get_gene_meta(melanoma, sthet_only=TRUE)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#
SThet = function(x=NULL, genes=NULL, samples=NULL, method='moran', k=NULL, overwrite=TRUE,
                 cores=NULL, verbose=TRUE){
  # Record time
  zero_t = Sys.time()
  if(verbose){
    cat(paste0('SThet started.\n'))
  }

  # Select sample names if NULL or if number entered
  samples = resolve_samples(x, samples, slot='tr_counts', numeric_handling='plain',
                            char_transform='none', post_check='none')

  # Check that genes have been input
  if(is.null(genes)){
    stop('Please enter one or more genes to calculate statistics.')
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in samples){
    # Check if gene names are in the data set
    subsetgenes = genes[genes %in% rownames(x@tr_counts[[i]])]
    combo = dplyr::bind_rows(combo, expand.grid(i, subsetgenes))

    # Get genes not present.
    notgenes = genes[!(genes %in% rownames(x@tr_counts[[i]]))]

    if(!rlang::is_empty(notgenes)){
      message(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i), ".\n")
    }

    rm(subsetgenes, notgenes) # Clean env

    # Add columns in gene meta data if not already present
    if(!('moran_i' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['moran_i']] = NA
    }
    if(!('geary_c' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['geary_c']] = NA
    }
  }

  # Check whether or not a list of weights have been created
  if(overwrite | is.null(x@misc[['sthet']][['listws']])){
    if(verbose){
      cat(paste("\tCalculating spatial weights...\n")) ## Mostly added to make sure calculation is happening only when needed.
    }
    if(!is.null(k)){
      k = as.integer(k)
      if(!is.na(k) & k > 0){
        x@misc[['sthet']][['listws']] = create_listw_from_knn(x, ks=k)
      } else{
        stop("If using k nearest-neighbors, please input a positive integer for k.")
      }
    } else{
      x@misc[['sthet']][['listws']] = create_listw_from_dist(x, cores=cores)
    }
  }

  # Perform calculations
  if('moran' %in% method){
    x = gene_moran_i_notest(x=x, combo=combo, overwrite=overwrite, cores=cores)
  }
  if('geary' %in% method){
    x = gene_geary_c_notest(x=x, combo=combo, overwrite=overwrite, cores=cores)
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('SThet completed in ', round(end_t, 2), ' min.\n'))
  }

  return(x)
}


# Helpers ----------------------------------------------------------------------

##
# @title gene_moran_i_dist
# @description Calculates Moran's I from ST data.
#
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Moran's I
#
#
gene_moran_i_notest = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  # Define cores available
  cores = resolve_cores(length(x@spatial_meta), cores, coerce='as_integer')
  if(is.na(cores)){
    stop('Could not recognize number of cores requested')
  }

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  stat_list = sge_lapply(seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))), function(i_combo){
    i = unique(as.vector(unlist(combo[[1]])))[i_combo]
    genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
    sthet_stat_one_sample(sample_tr_counts(x, i), genes_tmp, x@misc[['sthet']][['listws']][[i]], method='moran')
  }, cores=cores, preschedule=F)

  names(stat_list) = unique(as.vector(unlist(combo[[1]])))

  # Store kriging results in STList.
  #for(i in 1:nrow(combo)){
  for(i in names(stat_list)){
    for(j in names(stat_list[[i]])){
      #combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
      combo_name = c(i, j)
      if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i']))){
        # x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
        #x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]][[j]][['estimate']][1])
        x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]][[j]][['I']])
        #print(as.vector(stat_list[[i]]$estimate[1]))
      }
    }
  }
  return(x)
}


##
# @title gene_geary_c_dist
# @description Calculates Geary's C from ST data.
#
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Geary's I
#
#
gene_geary_c_notest = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  # Define cores available
  cores = resolve_cores(length(x@spatial_meta), cores, coerce='as_integer')
  if(is.na(cores)){
    stop('Could not recognize number of cores requested')
  }

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  stat_list = sge_lapply(seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))), function(i_combo){
    i = unique(as.vector(unlist(combo[[1]])))[i_combo]
    genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
    sthet_stat_one_sample(sample_tr_counts(x, i), genes_tmp, x@misc[['sthet']][['listws']][[i]], method='geary')
  }, cores=cores, preschedule=F)

  names(stat_list) = unique(as.vector(unlist(combo[[1]])))

  # Store kriging results in STList.
  #for(i in 1:nrow(combo)){
  for(i in names(stat_list)){
    for(j in names(stat_list[[i]])){
      #combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
      combo_name = c(i, j)
      if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c']))){
        # x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
        #x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c'] = as.vector(stat_list[[i]][[j]][['estimate']][1])
        x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c'] = as.vector(stat_list[[i]][[j]][['C']])
        #print(as.vector(stat_list[[i]]$estimate[1]))
      }
    }
  }
  return(x)
}


##
# @title sthet_stat_one_sample: Compute a spatial autocorrelation statistic
# for one sample's genes
# @description The per-sample unit of work dispatched by both
# `gene_moran_i_notest()`'s and `gene_geary_c_notest()`'s `sge_lapply` calls,
# extracted so it is directly callable/benchmarkable for one sample without
# going through the whole STlist. The two callers differ only in which
# `spdep` statistic they request.
#
# @param tr_counts_mat a transformed count matrix for one sample (genes in rows).
# @param genes the gene names to compute the statistic for.
# @param listw a `spdep` listw object of spatial weights for this sample.
# @param method one of 'moran' or 'geary'.
# @return a named list (by gene) of `spdep::moran`/`spdep::geary` result objects.
#
sthet_stat_one_sample = function(tr_counts_mat, genes, listw, method=c('moran', 'geary')){
  method = match.arg(method)
  stat_list_tmp = list()
  for(j in genes){
    # Extract expression data for a given gene.
    gene_expr = tr_counts_mat[j, ]

    # Estimate statistic.
    if(method == 'moran'){
      stat_est = spdep::moran(x=gene_expr,
                              listw=listw,
                              n=length(listw$neighbours),
                              S0=spdep::Szero(listw))
    } else{
      stat_est = spdep::geary(x=gene_expr,
                              listw=listw,
                              n=length(listw$neighbours),
                              n1=length(listw$neighbours)-1,
                              S0=spdep::Szero(listw))
    }

    stat_list_tmp[[j]] = stat_est
  }
  return(stat_list_tmp)
}

