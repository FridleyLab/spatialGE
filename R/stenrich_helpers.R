# Helpers for STenrich(): per-combination units of work (score/permutation
# test for one sample x gene-set pair) and the sample x gene-set score
# calculators that dispatch them.

##
# @title stenrich_score_one_combo: Compute average gene set expression for
# one sample x gene-set combination
# @description The per-combination unit of work dispatched by
# `calculate_gs_mean_exp()`'s `sge_lapply` call, extracted so it is directly
# callable/benchmarkable without going through the whole combination table.
#
# @param delayed_mat the sample's (possibly domain-subset) `DelayedArray` of
# transformed counts.
# @param pw_genes_combo the gene set's genes present in this sample's data.
# @param geneset_name the gene set's name, used as the returned row name.
# @param min_genes the minimum number of present genes required to compute
# an average (otherwise all spots/cells get NA).
# @return a one-row data frame of per-spot/cell average expression.
#
stenrich_score_one_combo = function(delayed_mat, pw_genes_combo, geneset_name, min_genes){
  # Test if genes in data set are enough
  if(length(pw_genes_combo) >= min_genes){
    # Calculate average expression of genes within pathway for each spot/cell
    expr_subset = delayed_mat[pw_genes_combo, ]
    pw_avg_exp = DelayedMatrixStats::colMeans2(expr_subset)
  } else{
    pw_avg_exp = rep(NA, ncol(delayed_mat))
    names(pw_avg_exp) = colnames(delayed_mat)
  }
  pw_avg_exp = data.frame(as.list(pw_avg_exp), check.names=F)
  rownames(pw_avg_exp) = geneset_name

  return(pw_avg_exp)
}


##
# @title stenrich_perm_one_combo: Run the spatial-enrichment permutation
# test for one sample x gene-set combination
# @description The per-combination unit of work dispatched by `STenrich()`'s
# `sge_lapply` call, extracted so it is directly callable/benchmarkable
# without going through the whole combination table.
#
# @param expr_vals the combination's per-spot/cell gene set expression values.
# @param coords_mat the sample's (possibly domain-subset) coordinate matrix.
# @param sample_name,gene_set_name the combination's sample and gene set names.
# @param size_test,size_gene_set the number of genes present vs. the gene set's total size.
# @param num_sds the number of standard deviations above the mean used as the high-expression threshold.
# @param min_units the minimum number of high-expression spots/cells required to test.
# @param reps the number of random permutations.
# @param seed the seed used before drawing the random permutations.
# @param verbose either logical or an integer (0, 1, or 2); values >= 2 print a per-combination progress message.
# @return a one-row data frame with the combination's test result.
#
stenrich_perm_one_combo = function(expr_vals, coords_mat, sample_name, gene_set_name, size_test, size_gene_set,
                                   num_sds, min_units, reps, seed, verbose){
  res_df = data.frame(sample_name=sample_name, gene_set=gene_set_name,
                      size_test=size_test,
                      size_gene_set=size_gene_set,
                      p_value=NA)

  # Perform tests
  if(!all(is.na(expr_vals))){ # Are all values for a gene set NA? (can happen if not enough genes in gene set, for example)
    if(verbose >= 2L){
      system(sprintf('echo "%s"', paste0("\tTesting sample ", sample_name, ", ", gene_set_name, "...")))
    }
    # Calculate expression or score threshold
    exp_thresh = mean(expr_vals, na.rm=T) + (num_sds*sd(expr_vals, na.rm=T))
    # Extract spots with average expression above sd threshold
    high_spots_bc = names(which(expr_vals >= exp_thresh))

    # Are there at least 'min_units' number of cells/spots?
    if(length(high_spots_bc) >= min_units){
      # Compute sum of distances among high expression spots
      coords_df_tmp = coords_mat[high_spots_bc, ]
      sum_high_distances = computeSubsampleSums(coords_df_tmp, n_subsample=nrow(coords_df_tmp), n_samples=1)
      rm(coords_df_tmp) # Clean env

      # Compute random distance permutations
      set.seed(seed)
      sum_rand_distances = computeSubsampleSums(coords_mat, n_subsample=length(high_spots_bc), n_samples=reps)

      # Compute p-value
      # Ho: sum of observed distances is larger than sum of random distances
      # Ha: sum of observed distances is smaller than sum of random distances
      count_test = sum(sum_rand_distances < sum_high_distances) # Times observed dist was higher than random dists
      p_val = count_test / reps

      res_df[['p_value']] = p_val
    }
  }

  return(res_df)
}


##
# calculate_gs_mean_exp
# The function calculates the average expression of genes in a gene set for a group of cells or spots
# The function works on an STlist and works by calculating average on a combination of samples and gene sets
# Returns a named list of data frames
#
calculate_gs_mean_exp = function(x=NULL, combo=NULL, pw_genes=NULL, min_genes=NULL, cores=NULL){
  # Loop through combinations of samples and gene sets
  result_df = sge_lapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[[1]][i])
    geneset_tmp = as.vector(combo[[2]][i])
    stenrich_score_one_combo(x[[sample_tmp]], pw_genes[[i]], geneset_tmp, min_genes)
  }, cores=cores)
  names(result_df) = paste0(combo[[1]], '&&', combo[[2]])

  # Create dataframes per sample
  result_df = lapply(1:length(unique(combo[[1]])), function(i){
    sample_tmp = unique(combo[[1]])[i]
    idx = which(combo[[1]] == sample_tmp)
    df_tmp = do.call(rbind, result_df[idx])
    rownames(df_tmp) = gsub(paste0('^', sample_tmp, '&&'), '', rownames(df_tmp))

    return(df_tmp)
  })
  names(result_df) = unique(combo[[1]])

  return(result_df)
}


##
# calculate_gs_gsva_score
# The function calculates the GSVA score of a gene set for a group of cells or spots
# The function works on an STlist and works by calculating the GSVA score on a series of samples and gene sets
# Returns a named list of data frames
#
calculate_gs_gsva_score = function(x=NULL, combo=NULL, gene_sets=NULL, pw_genes=NULL, min_genes=NULL, cores=NULL, verbose=T){
  # Loop through samples
  samples = as.vector(unique(combo[[1]]))
  result_df = lapply(1:length(samples), function(i){
    sample_tmp = as.vector(unique(samples))[i]

    if(!is.null(pw_genes)){
      pw_genes_tmp = pw_genes[combo[[1]] == sample_tmp]
      names(pw_genes_tmp) = as.vector(combo[[2]][combo[[1]] == sample_tmp])
    }
    gene_sets_tmp = gene_sets[ names(pw_genes_tmp)[unlist(lapply(pw_genes_tmp, length)) >= min_genes] ]

    # Calculate GSVA scores for each spot or cell
    if(length(gene_sets_tmp) > 0){
      gsvapar = GSVA::gsvaParam(as.array(x[[sample_tmp]]), geneSets=gene_sets_tmp)
      pw_gsva_exp = GSVA::gsva(gsvapar, BPPARAM=sge_bpparam(cores), verbose=verbose)
      pw_gsva_exp = as.data.frame(pw_gsva_exp)
    } else{
      pw_gsva_exp = data.frame(matrix(nrow=length(pw_genes_tmp), ncol=ncol(x[[sample_tmp]])))
      rownames(pw_gsva_exp) = names(pw_genes_tmp)
      colnames(pw_gsva_exp) = colnames(x[[sample_tmp]])
    }

    # Add rows with NA to mimic output from "average expression" approach
    if(length(gene_sets_tmp) != length(gene_sets)){
      pw_gsva_exp_list = lapply(names(gene_sets), function(j){
        if(j %in% rownames(pw_gsva_exp)){
          return(pw_gsva_exp[j, , drop=F])
        } else{
          df_tmp = as.data.frame(as.list(setNames(rep(NA, ncol(pw_gsva_exp)), colnames(pw_gsva_exp))), check.names=F)
          rownames(df_tmp) = j
          return(df_tmp)
        }
      })
      pw_gsva_exp = do.call(rbind, pw_gsva_exp_list)
    }

    return(pw_gsva_exp)
  })
  names(result_df) = unique(samples)

  return(result_df)
}
