##
# @title SThet: Computes global spatial autocorrelation statistics on gene expression
# @description Computes the global spatial autocorrelation statistics Moran's I and/or
# Geary's C for a set of genes
# @details The function computes global spatial autocorrelation statistics (Moran's I and/or
# Geary's C) for the requested genes and samples. Then computation uses the
# package `spdep`. The calculated statistics are stored in the STlist, which can
# be accessed with the `get_gene_meta` function. For visual comparative analysis,
# the function `compare_SThet` can be used afterwards.
#
# @param x an STlist
# @param genes a vector of gene names to compute statistics
# @param samples the samples to compute statistics
# @param method The spatial statistic(s) to estimate. It can be set to 'moran',
# 'geary' or both. Default is 'moran'
# @param k the number of neighbors to estimate weights. By default NULL, meaning that
# spatial weights will be estimated from Euclidean distances. If an positive integer is
# entered, then the faster k nearest-neighbors approach is used. Please keep in mind
# that estimates are not as accurate as when using the default distance-based method.
# @param overwrite logical indicating if previous statistics should be overwritten.
# Default to FALSE (do not overwrite)
# @param cores the number of cores to use during computations
# @return an STlist containing spatial statistics
#
# @export
#
SThet_invdist_test = function(x=NULL, genes=NULL, samples=NULL, method='moran', k=NULL, overwrite=T, cores=NULL){
  # Select sample names if NULL or if number entered
  if (is.null(samples)){
    samples = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samples = names(x@tr_counts)[samples]
    }
  }

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
    #cat("Calculating spatial weights...\n") ## Mostly added to make sure calculation is happening only when needed.
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
    x = gene_moran_i_dist(x=x, combo=combo, overwrite=overwrite, cores=cores)
  }
  if('geary' %in% method){
    x = gene_geary_c_dist(x=x, combo=combo, overwrite=overwrite, cores=cores)
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
gene_moran_i_dist = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  # Define cores available ### PARALLEL
  if(is.null(cores)){
    cores = count_cores(length(x@spatial_meta))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  #stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){ ### PARALLEL
  #stat_list = list()   #### WHEN NOT USING PARALLEL
  #for(i_combo in 1:nrow(combo)){   #### WHEN NOT USING PARALLEL
  stat_list = parallel::mclapply(seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))), function(i_combo){
    i = unique(as.vector(unlist(combo[[1]])))[i_combo]
    genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
    stat_list_tmp = list()
    for(j in genes_tmp){

      #j = as.vector(unlist(combo[i_combo, 2]))

      # Extract expression data for a given gene.
      gene_expr = x@tr_counts[[i]][j, ]

      # Estimate statistic.
      stat_est = spdep::moran.test(x=gene_expr, listw=x@misc[['sthet']][['listws']][[i]])

      #stat_list[[i_combo]] = stat_est
      stat_list_tmp[[j]] = stat_est
      #stat_list[[paste(i, j, sep='&&')]] = stat_est
      #return(stat_est) ### PARALLEL
    }
    return(stat_list_tmp)
  }, mc.cores=cores, mc.preschedule=F) ### PARALLEL

  #}
  #names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')
  names(stat_list) = unique(as.vector(unlist(combo[[1]])))

  # Store kriging results in STList.
  #for(i in 1:nrow(combo)){
  for(i in names(stat_list)){
    for(j in names(stat_list[[i]])){
      #combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
      combo_name = c(i, j)
      if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i']))){
        # x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
        x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]][[j]][['estimate']][1])
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
gene_geary_c_dist = function(x=NULL, combo=NULL, overwrite=T, cores=NULL){
  # Define cores available ### PARALLEL
  if(is.null(cores)){
    cores = count_cores(length(x@spatial_meta))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Use method to compute autocorrelation described in tutorial of https://rspatial.org/
  # Loop through combinations of samples x genes
  #stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){ ### PARALLEL
  #stat_list = list()   #### WHEN NOT USING PARALLEL
  #for(i_combo in 1:nrow(combo)){   #### WHEN NOT USING PARALLEL
  stat_list = parallel::mclapply(seq_along(1:length(unique(as.vector(unlist(combo[[1]]))))), function(i_combo){
    i = unique(as.vector(unlist(combo[[1]])))[i_combo]
    genes_tmp = unique(as.vector(unlist(combo[[2]][combo[[1]] == i])))
    stat_list_tmp = list()
    for(j in genes_tmp){

      #j = as.vector(unlist(combo[i_combo, 2]))

      # Extract expression data for a given gene.
      gene_expr = x@tr_counts[[i]][j, ]

      # Estimate statistic.
      stat_est = spdep::geary.test(x=gene_expr, listw=x@misc[['sthet']][['listws']][[i]])

      #stat_list[[i_combo]] = stat_est
      stat_list_tmp[[j]] = stat_est
      #stat_list[[paste(i, j, sep='&&')]] = stat_est
      #return(stat_est) ### PARALLEL
    }
    return(stat_list_tmp)
  }, mc.cores=cores, mc.preschedule=F) ### PARALLEL

  #}
  #names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')
  names(stat_list) = unique(as.vector(unlist(combo[[1]])))

  # Store kriging results in STList.
  #for(i in 1:nrow(combo)){
  for(i in names(stat_list)){
    for(j in names(stat_list[[i]])){
      #combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
      combo_name = c(i, j)
      if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c']))){
        # x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
        x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c'] = as.vector(stat_list[[i]][[j]][['estimate']][1])
        #print(as.vector(stat_list[[i]]$estimate[1]))
      }
    }
  }
  return(x)
}

