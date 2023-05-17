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
SThet_3 = function(x=NULL, genes=NULL, samples=NULL, method='moran'){
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

  if('moran' %in% method){
    x = gene_moran_i_knn(x=x, combo=combo_tmp) # Capital I instead to differentiate from new implementation
  }
  if('geary' %in% method){
    x = gene_geary_c_knn(x=x, combo=combo_tmp) # Capital C instead to differentiate from new implementation
  }
  # if('getisord' %in% method){
  #   x = gene_getisord_gi(x=x, combo=combo_tmp, adjmtx=adj_ls)
  # }

  return(x)
}


# Helpers ----------------------------------------------------------------------

##
# @title gene_moran_I
# @description Calculates Moran's I from ST data.
# @details
# This function takes an STList and a vector with HUGO gene names and returns
# Morans' I for each element of the vector.
#
# @param x, an STList with normalized gene counts.
# @param genes, a vector with gene names in the normalized count matrix.
# @param who, the indexes of the spatial arrays for which the statistic
# will be calculated.
# @return x, a STList including the values corresponding to Moran's I for each
# gene in genes.
#
# @export
#
#
gene_moran_i_knn <- function(x=NULL, combo=NULL) {

  genes = as.vector(unique(combo[[2]]))
  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_het[[gene]]) && rlang::is_empty(x@gene_het[[gene]])){
      x@gene_het[[gene]] = list()
      for(i in 1:length(x@tr_counts)){
        x@gene_het[[gene]][[names(x@tr_counts[i])]] = list(morans_I=NULL,
                                                           gearys_C=NULL#,
                                                           #getis_ord_Gi=NULL
                                                           )
      }
    }
  }

  # Check whether or not a list of weights have been created
  if(is.null(x@misc$listws)){
    x@misc$listws = create_listw_from_knn(x, ks=4)
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Create distance matrix based on the coordinates of each sampled location.
    # subj_dists = as.matrix(dist(x@coords[[i]][2:3]))
    # subj_dists[subj_dists == 0] = 0.0001
    # subj_dists_inv = 1/subj_dists
    # diag(subj_dists_inv) = 0

    # Extract expression data for a given gene.
    gene_expr = x@tr_counts[[i]][j, ]

    # Estimate statistic.
    stat_est = spdep::moran.test(x=gene_expr, listw=x@misc$listws[[i]])

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
# @details
# This function takes an STList and a vector with HUGO gene names and returns
# Geary's C for each element of the vector.
#
# @param x, an STList with normalized gene counts.
# @param genes, a vector with gene names in the normalized count matrix.
# @param who, the index of the spatial array for which the statistic
# will be calculated.
# @return x, a STList including the values corresponding to Geary's C for each
# gene in genes.
#
# @export
#
#
gene_geary_c_knn <- function(x=NULL, combo=NULL) {

  genes = as.vector(unique(combo[[2]]))
  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_het[[gene]]) && rlang::is_empty(x@gene_het[[gene]])){
      x@gene_het[[gene]] = list()
      for(i in 1:length(x@tr_counts)){
        x@gene_het[[gene]][[names(x@tr_counts[i])]] = list(morans_I=NULL,
                                                           gearys_C=NULL#,
                                                           #getis_ord_Gi=NULL
                                                           )
      }
    }
  }

  # Check whether or not a list of weights have been created
  if(is.null(x@misc$listws)){
    x@misc$listws = create_listw_from_knn(x, ks=4)
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Create distance matrix based on the coordinates of each sampled location.
    # subj_dists = as.matrix(dist(x@coords[[i]][2:3]))
    # subj_dists[subj_dists == 0] = 0.0001
    # subj_dists_inv = 1/subj_dists
    # diag(subj_dists_inv) = 0

    # Extract expression data for a given gene.
    gene_expr = x@tr_counts[[i]][j, ]

    # Estimate statistic.
    stat_est = spdep::geary.test(x=gene_expr, listw=x@misc$listws[[i]])

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


##
# @title gene_getis_Gi
# @description Calculates Getis-Ord Gi C from ST data.
# @details
# This function takes an STList and a vector with HUGO gene names and returns
# Getis-Ord Gi for each element of the vector.
#
# @param x, an STList with normalized gene counts.
# @param genes, a vector with gene names in the normalized count matrix.
# @param who, the indexes of the spatial arrays for which the statistic
# will be calculated.
# @return x, a STList including the values corresponding to Getis-Ord Gi for each
# gene in genes.
#
# @export
#
#
# gene_getis_Gi <- function(x=NULL, genes=NULL, who=NULL) {
#   # Test if no specific subject plot was requested.
#   if (is.null(who)) {
#     who = c(1:length(x@tr_counts))
#   }
#
#   # Generate combination of sample x gene to for.
#   combo = tibble::tibble()
#   for(i in who){
#     subsetgenes_mask = genes %in% x@tr_counts[[i]]$gene
#     subsetgenes = genes[subsetgenes_mask]
#     combo = dplyr::bind_rows(combo, expand.grid(names(x@tr_counts[i]), subsetgenes))
#
#     # Get genes not present.
#     notgenes = genes[!subsetgenes_mask]
#
#     if(!rlang::is_empty(notgenes)){
#       cat(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", names(x@tr_counts[i]), ".\n"))
#     }
#   }
#
#   # Test if list with kriging exists for each gene. If not create it.
#   for(gene in genes){
#     if(is.null(x@gene_het[[gene]]) && rlang::is_empty(x@gene_het[[gene]])){
#       x@gene_het[[gene]] = list()
#       for(i in 1:length(x@tr_counts)){
#         x@gene_het[[gene]][[names(x@tr_counts[i])]] = list(morans_I=NULL,
#                                                            gearys_C=NULL,
#                                                            getis_ord_Gi=NULL)
#       }
#     }
#   }
#
#   # Check whether or not a list of weights have been created
#   if(is.null(x@misc$listws)){
#     x@misc$listws = create_listw(x)
#   }
#
#   # Define cores available
#   cores = count_cores(nrow(combo))
#   # Loop through combinations of samples x genes
#   stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
#     i = as.vector(unlist(combo[i_combo, 1]))
#     j = as.vector(unlist(combo[i_combo, 2]))
#
#     # Create distance matrix based on the coordinates of each sampled location.
#     # subj_dists = as.matrix(dist(x@coords[[i]][2:3]))
#     # subj_dists[subj_dists == 0] = 0.0001
#     # subj_dists_inv = 1/subj_dists
#     # diag(subj_dists_inv) = 0
#
#     # Extract expression data for a given gene.
#     gene_expr = x@tr_counts[[i]][j, ]
#
#     # Estimate statistic.
#     stat_est = spdep::globalG.test(x=gene_expr, listw=x@misc$listws[[i]])
#
#     return(stat_est)
#
#   }, mc.cores=cores, mc.preschedule=F)
#   names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')
#
#   # Store kriging results in STList.
#   for(i in 1:nrow(combo)){
#     combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
#     x@gene_het[[combo_name[2]]][[combo_name[1]]][['getis_ord_Gi']] = stat_list[[i]]
#   }
#
#   return(x)
# }

