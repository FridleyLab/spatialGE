##
#' @title SThet: Computes global spatial autocorrelation statistics on gene expression
#' @description Computes the global spatial autocorrelation statistics Moran's I and/or
#' Geary's C for a set of genes
#' @details The function computes global spatial autocorrelation statistics (Moran's I and/or
#' Geary's C) for the requested genes and samples. Then computation uses the
#' package `spdep`. The calculated statistics are stored in the STlist, which can
#' be accessed with the `get_gene_meta` function.
#'
#' @param x an STlist
#' @param genes a vector of gene names to compute statistics
#' @param samples the samples to compute statistics
#' @param method The spatial statistic(s) to estimate. It can be set to 'moran',
#' 'geary' or both. Default is 'moran'
#' @param overwrite logical indicating if previous statistics should be overwritten.
#' Default to FALSE (do not overwrite)
#' @return an STlist containing spatial statistics
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
#' count_files <- grep("counts", data_files, value=T)
#' coord_files <- grep("mapping", data_files, value=T)
#' clin_file <- grep("thrane_clinical", data_files, value=T)
#' melanoma <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#' melanoma <- transform_data(melanoma, method='log')
#' melanoma <- SThet(melanoma, genes=c('MLANA', 'TP53'), method='moran')
#' get_gene_meta(melanoma, sthet_only=T)
#'
#' @export
#'
SThet = function(x=NULL, genes=NULL, samples=NULL, method='moran', overwrite=F){
  # Select sample names if NULL or if number entered
  if (is.null(samples)){
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
  }

  # Generate combination of sample x gene to for.
  combo_tmp = tibble::tibble()
  for(i in samples){
    # Check if gene names are in the data set
    subsetgenes = genes[genes %in% rownames(x@tr_counts[[i]])]
    combo_tmp = dplyr::bind_rows(combo_tmp, expand.grid(i, subsetgenes))

    # Get genes not present.
    notgenes = genes[!(genes %in% rownames(x@tr_counts[[i]]))]

    if(!rlang::is_empty(notgenes)){
      cat(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i), ".\n")
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

  # Perform calculations
  if('moran' %in% method){
    x = gene_moran_I(x=x, combo=combo_tmp, overwrite=overwrite) # Capital I instead to differentiate from new implementation
  }
  if('geary' %in% method){
    x = gene_geary_C(x=x, combo=combo_tmp, overwrite=overwrite) # Capital C instead to differentiate from new implementation
  }

  return(x)
}


# Helpers ----------------------------------------------------------------------

##
# @title gene_moran_I
# @description Calculates Moran's I from ST data.
#
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Moran's I
#
#
gene_moran_I <- function(x=NULL, combo=NULL, overwrite=overwrite) {

  genes = as.vector(unique(combo[[2]]))

  # Check whether or not a list of weights have been created
  if(is.null(x@misc[['sthet']][['listws']])){
    x@misc[['sthet']][['listws']] = create_listw(x)
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Extract expression data for a given gene.
    gene_expr = x@tr_counts[[i]][j, ]

    # Estimate statistic.
    stat_est = spdep::moran.test(x=gene_expr, listw=x@misc[['sthet']][['listws']][[i]])

    return(stat_est)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
  }

  return(x)
}


##
# @title gene_geary_C
# @description Calculates Geary's C from ST data.
#
# @param x an STlist with normalized gene counts.
# @param combo a table with combinations of samples and genes to calculate statistics
# @return x a STlist including the calculated Geary's I
#
#
gene_geary_C <- function(x=NULL, combo=NULL, overwrite=overwrite) {

  genes = as.vector(unique(combo[[2]]))

  # Check whether or not a list of weights have been created
  if(is.null(x@misc[['sthet']][['listws']])){
    x@misc[['sthet']][['listws']] = create_listw(x)
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Extract expression data for a given gene.
    gene_expr = x@tr_counts[[i]][j, ]

    # Estimate statistic.
    stat_est = spdep::geary.test(x=gene_expr, listw=x@misc[['sthet']][['listws']][[i]])

    return(stat_est)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c'] = as.vector(stat_list[[i]]$estimate[1])
  }

  return(x)
}

