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
#
gene_getis_Gi <- function(x=NULL, genes=NULL, who=NULL) {
  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who = c(1:length(x@tr_counts))
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in who){
    subsetgenes_mask = genes %in% x@tr_counts[[i]]$gene
    subsetgenes = genes[subsetgenes_mask]
    combo = dplyr::bind_rows(combo, expand.grid(names(x@tr_counts[i]), subsetgenes))

    # Get genes not present.
    notgenes = genes[!subsetgenes_mask]

    if(!rlang::is_empty(notgenes)){
      cat(paste0(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", names(x@tr_counts[i]), ".\n"))
    }
  }

  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_het[[gene]]) && rlang::is_empty(x@gene_het[[gene]])){
      x@gene_het[[gene]] = list()
      for(i in 1:length(x@tr_counts)){
        x@gene_het[[gene]][[names(x@tr_counts[i])]] = list(morans_I=NULL,
                                                           gearys_C=NULL,
                                                           getis_ord_Gi=NULL)
      }
    }
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  stat_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Create distance matrix based on the coordinates of each sampled location.
    subj_dists = as.matrix(dist(x@coords[[i]][2:3]))
    subj_dists_inv = 1/subj_dists
    diag(subj_dists_inv) = 0

    # Extract expression data (voom counts) for a given gene.
    gene_expr = unlist(x@tr_counts[[i]][x@tr_counts[[i]][[1]] == j, -1])

    # Estimate statistic.
    stat_est = spdep::globalG.test(gene_expr, spdep::mat2listw(subj_dists_inv, style='B'))

    return(stat_est)
  }, mc.cores=cores, mc.preschedule=F)
  names(stat_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(stat_list)[i], split = '&&'))
    x@gene_het[[combo_name[2]]][[combo_name[1]]][['getis_ord_Gi']] = stat_list[[i]]
  }

  return(x)
}

