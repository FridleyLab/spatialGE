##
# @title SThet: Computes global spatial autocorrelation statistics on gene expression
# @description Computes the global spatial autocorrelation statistics Moran's I and/or
# Geary's C for a set of genes
# @details The function computes global spatial autocorrelation statistics (Moran's I and/or
# Geary's C) for the requested genes and samples. Then computation uses the
# package `spdep`. The calculated statistics are stored in the STlist, which can
# be accessed with the `get_gene_meta` function.
#
# @param x an STlist
# @param genes a vector of gene names to compute statistics
# @param samples the samples to compute statistics
# @param method The spatial statistic(s) to estimate. It can be set to 'moran',
# 'geary' or both. Default is 'moran'
# @param overwrite logical indicating if previous statistics should be overwritten.
# Default to TRUE (overwrite)
# @return an STlist containing spatial statistics
#
# @export
#
old_SThet = function(x=NULL, genes=NULL, samples=NULL, method='moran', overwrite=T){
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

  # Check whether or not a list of weights have been created
  if(overwrite | is.null(x@misc[['sthet']][['listws']])){
    x@misc[['sthet']][['listws']] = create_listw(x)
  }

  # Perform calculations
  if('moran' %in% method){
    x = gene_moran_i(x=x, combo=combo_tmp, overwrite=overwrite)
  }
  if('geary' %in% method){
    x = gene_geary_c(x=x, combo=combo_tmp, overwrite=overwrite)
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
gene_moran_i <- function(x=NULL, combo=NULL, overwrite=T) {

  genes = as.vector(unique(combo[[2]]))

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
    if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i']))){
      x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'moran_i'] = as.vector(stat_list[[i]]$estimate[1])
      #print(as.vector(stat_list[[i]]$estimate[1]))
    }
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
gene_geary_c <- function(x=NULL, combo=NULL, overwrite=T) {

  genes = as.vector(unique(combo[[2]]))

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
    if(overwrite | is.na(as.vector(x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c']))){
      x@gene_meta[[combo_name[1]]][x@gene_meta[[combo_name[1]]][['gene']] == combo_name[2], 'geary_c'] = as.vector(stat_list[[i]]$estimate[1])
      #print(as.vector(stat_list[[i]]$estimate[1]))
    }
  }

  return(x)
}

