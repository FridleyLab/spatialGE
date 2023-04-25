##
#' @title SThet: Compute gene-wise global spatial autocorrelation statistics
#' @description Computes Moran's I and/or Geary's C for a set of genes
#' @details The function computes global spatial autocorrelation statistics for the
#' requested genes. The spatial weights are defined by the minimum distance between
#' cells/spots. Spots within the two times the minimum distance are defined as neighbors and
#' a binary matrix is used as spatial weights for computation
#'
#' @param x an STlist
#' @param genes a vector of gene names to compute statistics
#' @param samples the samples to compute statistics
#' @param method The spatial statistic(s) to estimate. Default is 'moran',
#' @param dist_thr the number of euclidean distance units to consider a cell/spot a neighbor.
#' but a vector with any combination of 'moran' or 'geary' can be passed.
#'
#' @export
#'
#
#
SThet_2 = function(x=NULL, genes=NULL, samples=NULL, method='moran', dist_thr=20){
  # Select sample names if NULL or if number entered
  if (is.null(samples)){
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
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
      cat(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", i), ".\n")
    }

    rm(subsetgenes, notgenes) # Clean env
  }

  # Add columns in gene meta data
  for(i in samples){
    if(!('moran_i' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['moran_i']] = NA
    }
    if(!('geary_c' %in% colnames(x@gene_meta[[i]]))){
      x@gene_meta[[i]][['geary_c']] = NA
    }

    # Deactivated because this could be implemented as a plotting function instead
    # The function could provide detection of hot or coldspots
    # if(!('getisord_gi' %in% colnames(x@gene_meta[[i]]))){
    #   x@gene_meta[[i]][['getisord_gi']] = NA
    # }
  }

  # Create binary 1/0 adjacency matrix
  # Use 2 x euclidean distance to define neighbors
  dist_thr = as.numeric(dist_thr)
  adj_ls = list()
  for(i in samples){
    adj = as.matrix(dist(x@spatial_meta[[i]][, c('xpos', 'ypos')]))
    diag(adj) = NA
    colnames(adj) = x@spatial_meta[[i]][['libname']]
    rownames(adj) = x@spatial_meta[[i]][['libname']]
    min_tmp = min(adj, na.rm=T)
    adj[adj >= min_tmp * 0.5 & adj <= min_tmp * dist_thr] = 1
    adj[adj > min_tmp * 2] = 0
    diag(adj) = 0

    adj_ls[[i]] = adj

    rm(min_tmp, adj) # Clean env
  }

  if('moran' %in% method){
    x = gene_moran_i(x=x, combo=combo, adjmtx=adj_ls)
  }
  if('geary' %in% method){
    x = gene_geary_c(x=x, combo=combo, adjmtx=adj_ls)
  }
  # if('getisord' %in% method){
  #   x = gene_getisord_gi(x=x, combo=combo, adjmtx=adj_ls)
  # }

  return(x)
}


# Helpers ----------------------------------------------------------------------

##
# @title gene_moran_i
# @description Calculates Moran's I from ST data.
# @details
# This function takes an STlist and a vector with HUGO gene names and returns
# Morans' I stored in the gene_meta slot.
#
# @param x an STlist with normalized gene counts.
# @param combo data frame with combination of samples and genes to analyze
# @param adjmtx a list of adjacency matrices. Each element of the list represents
# one sample
# @return x a STlist including the values corresponding to Moran's I for each
# gene in genes.
#
#
gene_moran_i = function(x=NULL, combo=NULL, adjmtx=NULL){
  # Define cores available
  cores = count_cores(nrow(combo))

  # Use method to compute Moran's I described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get expression values
    y = x@tr_counts[[i]][j, ][match(colnames(adjmtx[[i]]), colnames(x@tr_counts[[i]]))]

    # Compute Moran's I
    moran_val = terra::autocor(x=y, w=adjmtx[[i]], method="moran")

    return(moran_val)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_meta[[ combo_name[1] ]][ x@gene_meta[[ combo_name[1] ]][['gene']] == combo_name[2], 'moran_i'] = stat_list[[i]]
  }

  return(x)
}


##
# @title gene_geary_c
# @description Calculates Geary's C from ST data.
# @details
# This function takes an STlist and a vector with HUGO gene names and returns
# Geary's C stored in the gene_meta slot.
#
# @param x an STlist with normalized gene counts.
# @param combo data frame with combination of samples and genes to analyze
# @param adjmtx a list of adjacency matrices. Each element of the list represents
# one sample
# @return x a STlist including the values corresponding to Geary's C for each
# gene in genes.
#
#
gene_geary_c = function(x=NULL, combo=NULL, adjmtx=NULL){
  # Define cores available
  cores = count_cores(nrow(combo))

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get expression values
    y = x@tr_counts[[i]][j, ][match(colnames(adjmtx[[i]]), colnames(x@tr_counts[[i]]))]

    # Compute Geary's C
    geary_val = terra::autocor(x=y, w=adjmtx[[i]], method="geary")

    return(geary_val)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_meta[[ combo_name[1] ]][ x@gene_meta[[ combo_name[1] ]][['gene']] == combo_name[2], 'geary_c'] = stat_list[[i]]
  }

  return(x)
}


##
# @title gene_getisord_gi
# @description Calculates Getis-Ord Gi from ST data.
# @details
# This function takes an STlist and a vector with HUGO gene names and returns
# Getis-Ord Gi stored in the gene_meta slot.
#
# @param x an STlist with normalized gene counts.
# @param combo data frame with combination of samples and genes to analyze
# @param adjmtx a list of adjacency matrices. Each element of the list represents
# one sample
# @return x a STlist including the values corresponding to Getis-Ord Gi for each
# gene in genes.
#
#
gene_getisord_gi = function(x=NULL, combo=NULL, adjmtx=NULL){
  # Define cores available
  cores = count_cores(nrow(combo))

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get expression values
    y = x@tr_counts[[i]][j, ][match(colnames(adjmtx[[i]]), colnames(x@tr_counts[[i]]))]

    # Compute Geary's C
    getisord_val = terra::autocor(x=y, w=adjmtx[[i]], method="Gi")

    return(getisord_val)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_meta[[ combo_name[1] ]][ x@gene_meta[[ combo_name[1] ]][['gene']] == combo_name[2], 'getisord_gi'] = stat_list[[i]]
  }

  return(x)
}

