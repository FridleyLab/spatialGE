##
#' @title STenrich
#' @description Test for spatial enrichment of gene expression sets in ST data sets
#' @details The function performs a randomization test to assess if the sum of
#' distances between cells/spots with high expression of a gene set is lower than
#' the sum of distances among randomly selected cells/spots. The cells/spots are
#' considered as having high gene set expression if the average expression of genes in a
#' set is higher than the average expression plus `num_sds` times the standard deviation.
#' Control over the size of regions with high expression is provided by setting the
#' minimum number of cells/spots (`min_units`). This method is a modification of
#' the method devised by Hunter et al. 2021 (zebrafish melanoma study).
#'
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test. The names of the list should
#' identify the gene sets to be tested
#' @param score_type Controls how gene set expression is calculated. The options are
#' the average expression among genes in a set ('avg'), or a GSEA score ('gsva'). The
#' default is 'avg'
#' @param reps the number of random samples to be extracted. Default is 1000 replicates
#' @param annot name of the annotation within `x@spatial_meta` containing the spot/cell
#' categories. Needs to be used in conjunction with `domain`
#' @param domain the domain to restrict the analysis. Must exist within the spot/cell
#' categories included in the selected annotation (i.e., `annot`)
#' @param num_sds the number of standard deviations to set the minimum gene set
#' expression threshold. Default is one (1) standard deviation
#' @param min_units Minimum number of spots with high expression of a pathway for
#' that gene set to be considered in the analysis. Defaults to 20 spots or cells
#' @param min_genes the minimum number of genes of a gene set present in the data set
#' for that gene set to be included. Default is 5 genes
#' @param pval_adj_method the method for multiple comparison adjustment of p-values.
#' Options are the same as that of `p.adjust`. Default is 'BH'
#' @param seed the seed number for the selection of random samples. Default is 12345
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @param verbose either logical or an integer (0, 1, or 2) to increase verbosity
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats p.adjust
#
#
STenrich = function(x=NULL, samples=NULL, gene_sets=NULL, score_type='avg', reps=1000,
                    annot=NULL, domain=NULL, num_sds=1, min_units=20, min_genes=5,
                    pval_adj_method='BH', seed=12345, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()

  if(verbose){
    cat("Running STenrich...\n")
  }

  reps = as.integer(reps)
  num_sds = as.double(num_sds)
  min_units = as.integer(min_units)
  min_genes = as.integer(min_genes)

  # Define samples using names (convert indexes to names if necessary)
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = as.vector(na.omit(names(x@spatial_meta)[samples]))
    } else{
      samples = samples[samples %in% names(x@spatial_meta)]
    }
    # Verify that sample names exist
    if(length(samples) == 0 | !any(samples %in% names(x@spatial_meta))){
      raise_err(err_code="error0020")
    }
  }

  # Check if user request analysis within a specific domain(s)
  # Select spots/cells from specific domain(s)
  tissue_spots = NULL
  sample_rm = c() # Remove samples if annotation or domain(s) are not present
  if(!is.null(annot) & !is.null(domain)){
    if(verbose){
      cat('\tUsing domain-restricted mode.\n')
    }
    tissue_spots = list()
    for(i in samples){
      if(any(colnames(x@spatial_meta[[i]]) == annot)){
        if(any(x@spatial_meta[[i]][[annot]] %in% domain)){
          tissue_spots_tmp = x@spatial_meta[[i]][[1]][ x@spatial_meta[[i]][[annot]] %in% domain ]
          if(length(tissue_spots_tmp) >= min_units){
            tissue_spots[[i]] = tissue_spots_tmp
          } else{
            message(paste0('\tSample ', i, ' has less than two spots/cells assigned to domain. Skipping.\n'))
            sample_rm = append(sample_rm, i)
          }
          rm(tissue_spots_tmp) # Clean env
        } else{
          message(paste0('\tSample ', i, ' does not contain the specified domains. Skipping.\n'))
          sample_rm = append(sample_rm, i)
        }
      } else{
        message(paste0('\tSample ', i, ' does not contain the specified annotation. Skipping.\n'))
        sample_rm = append(sample_rm, i)
      }
    }
  }
  samples = samples[!(samples %in% sample_rm)]

  # Check if no samples were left in data set after checking for annotations
  if(length(samples) < 1){
    raise_err(err_code="error0021")
  }

  # Define number of cores for parallelization of tests
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(samples))
    user_cores = F
  } else{
    cores = ceiling(cores)
    user_cores = T
  }

  # Extract spots to be used in analysis (in case analysis restricted to domain was requested)
    # Extract spot coordinates and match order of spots
  coords_df = lapply(samples, function(i){
    df_tmp = x@spatial_meta[[i]][, c('libname', 'xpos', 'ypos')]
    if(length(tissue_spots) > 0){
      df_tmp = df_tmp[df_tmp[['libname']] %in% tissue_spots[[i]], ]
    }
    df_tmp = df_tmp %>% tibble::column_to_rownames(var='libname')
    df_tmp = as.matrix(df_tmp)

    return(df_tmp)
  })
  names(coords_df) = samples

  # Make sure temp directory is empty (needed for out-of-memory operations)
  #  unlink(dir(tempdir(), full.names=T), force=T, recursive=T)
  #library('DelayedArray')
  #library('HDF5Array')
  # Loop through samples in STlist and calculate distances
  # dists_mtx = parallel::mclapply(samples, function(i){
  #   system(sprintf('echo "%s"', paste0("\tCalculating distances in sample: ", i, "...")))
  #
  #
  #
  #
  #
  #   # Calculate distances for the sample
  #   distances_spots = calculate_euclidean_distances(coords_df)
  #
  #   return(distances_spots)
  # }, mc.cores=cores)
  # names(dists_mtx) = samples

  # Create data frame with combinations of samples and gene sets
  combo = expand.grid(sample_name=samples, gene_set=names(gene_sets))

  # Find genes from each gene set that are present in data
  pw_genes = lapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    pw_genes_tmp = unique(gene_sets[[gs_tmp]])
    pw_genes_tmp = intersect(pw_genes_tmp, rownames(x@tr_counts[[sample_tmp]]))

    return(pw_genes_tmp)
  })

  # Make sure valid score type was input
  score_type = tolower(score_type)
  if(!(score_type %in% c('avg', 'gsva'))){
    warning('Only `avg` or `gsva` are allowed in `score_type`. Using `avg`.')
    score_type = 'avg'
  }

  # Convert expression data to delayed matrix
  delayed_x = lapply(samples, function(i){
    mtx_tmp = x@tr_counts[[i]]
    if(length(tissue_spots) > 0){
      mtx_tmp = mtx_tmp[, tissue_spots[[i]]]
    }
    mtx_tmp = DelayedArray::DelayedArray(mtx_tmp)
    return(mtx_tmp)
  })
  names(delayed_x) = samples

  # Clear memory by offloading STlist
  x = NULL

  # Check requested score type (mean or GSVA) and calculate average expresion or scores
  if(score_type == 'avg'){
    if(verbose){
      cat("\tCalculating average gene set expression...\n")
    }
    # Define if user input number of cores, otherwise determine
    if(!user_cores){
      cores = count_cores(nrow(combo))
    }
    # Calculate average gene set expression
    result_df = calculate_gs_mean_exp(delayed_x, combo, pw_genes, min_genes, cores)

  } else if(score_type == 'gsva'){
    if(verbose){
      cat("\tCalculating GSVA score...\n")
    }
    # Define if user input number of cores, otherwise determine
    if(!user_cores){
      cores = count_cores(length(gene_sets))
    }
    # Calculate GSVA scores
    result_df = calculate_gs_gsva_score(delayed_x, combo, gene_sets=gene_sets, pw_genes=pw_genes, min_genes=min_genes, cores=cores)
  }

  # Find spots that highly express this pathway (mean + stdev_t*std in this case), calculate observed sum distances, and perform permutations
  # Define if user input number of cores, otherwise determine
  if(!user_cores){
    cores = count_cores(nrow(combo))
  }
  pval_res = parallel::mclapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    expr_vals = unlist(result_df[[sample_tmp]][gs_tmp, ])

    res_df = data.frame(sample_name=sample_tmp, gene_set=gs_tmp,
                        size_test=length(pw_genes[[i]]),
                        size_gene_set=length(gene_sets[[gs_tmp]]),
                        p_value=NA)

    # Perform tests
    if(!all(is.na(expr_vals))){ # Are all values for a gene set NA? (can happen if not enough genes in gene set, for example)
      if(verbose >= 2L){
        system(sprintf('echo "%s"', paste0("\tTesting sample ", sample_tmp, ", ", gs_tmp, "...")))
      }
      # Calculate expression or score threshold
      exp_thresh = mean(expr_vals, na.rm=T) + (num_sds*sd(expr_vals, na.rm=T))
      # Extract spots with average expression above sd threshold
      high_spots_bc = names(which(expr_vals >= exp_thresh))

      # Are there at least 'min_units' number of cells/spots?
      if(length(high_spots_bc) >= min_units){
        # Compute sum of distances among high expression spots
        #sum_high_distances = calculate_sum_distances(dists_mtx=dists_mtx[[sample_tmp]], select_spots=high_spots_bc)
        coords_df_tmp = coords_df[[sample_tmp]][high_spots_bc, ]
        sum_high_distances = computeSubsampleSums(coords_df_tmp, n_subsample=nrow(coords_df_tmp), n_samples=1)
        rm(coords_df_tmp) # Clean env

        # Compute random distance permutations
        set.seed(seed)
        # sum_rand_distances = lapply(1:reps, function(rep){
        #   rand_idx = sample(x=colnames(dists_mtx[[sample_tmp]]), size=length(high_spots_bc))
        #   rand_dists = calculate_sum_distances(dists_mtx=dists_mtx[[sample_tmp]], select_spots=rand_idx)
        #   return(rand_dists)
        # })
        # sum_rand_distances = unlist(sum_rand_distances)
        sum_rand_distances = computeSubsampleSums(coords_df[[sample_tmp]], n_subsample=length(high_spots_bc), n_samples=reps)

        # Compute p-value
        # Ho: sum of observed distances is larger than sum of random distances
        # Ha: sum of observed distances is smaller than sum of random distances
        count_test = sum(sum_rand_distances < sum_high_distances) # Times observed dist was higher than random dists
        p_val = count_test / reps

        res_df[['p_value']] = p_val #data.frame(sample_name=sample_tmp, gene_set=gs_tmp, size_gene_set=length(gene_sets[[gs_tmp]]), p_value=p_val)

        rm(high_spots_bc, sum_high_distances, sum_rand_distances, count_test, p_val) # Clean env
      }
    }

    return(res_df)
  }, mc.cores=cores)

  # Compile results in a single data frame
  pval_res = do.call(rbind, pval_res)

  # Calculate gene set "coverage"
  pval_res = pval_res %>%
    tibble::add_column(prop_size_test=.[['size_test']]/.[['size_gene_set']], .before='p_value') %>%
    dplyr::mutate(prop_size_test=round(prop_size_test, 3))

  # Split result among samples (for compatibility with old STenrich implementation)
  # Also, adjust p-values
  sample_names_tmp = unique(pval_res[['sample_name']])
  pval_res = lapply(sample_names_tmp, function(i){
    df_tmp = pval_res[pval_res[['sample_name']] == i, ]
    df_tmp[['adj_p_value']] = p.adjust(df_tmp[['p_value']], method=pval_adj_method)
    df_tmp = df_tmp[order(df_tmp[['adj_p_value']]), ]

    return(df_tmp)
  })
  names(pval_res) = sample_names_tmp

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STenrich completed in ', round(end_t, 2), ' min.\n'))
  }

  return(pval_res)
}


# Helpers ----------------------------------------------------------------------

##
# calculate_euclidean_distances
# The function returns a matrix of euclidean distances realized in HDF5Array out-of-memory format
# Returns an HDF5Array
#
calculate_euclidean_distances = function(coords_df=NULL){
  distances_spots = DelayedArray::DelayedArray(coords_df)
  distances_spots = DelayedArray::DelayedArray(as.matrix(stats::dist(distances_spots[, c('xpos', 'ypos')], method='euclidean')))
  #distances_spots = DelayedArray::realize(distances_spots, BACKEND='HDF5Array')
  distances_spots = DelayedArray::realize(distances_spots, BACKEND='RleArray')

  return(distances_spots)
}


##
# calculate_gs_mean_exp
# The function calculates the average expression of genes in a gene set for a group of cells or spots
# The function works on an STlist and works by calculating average on a combination of samples and gene sets
# Returns a named list of data frames
#
calculate_gs_mean_exp = function(x=NULL, combo=NULL, pw_genes=NULL, min_genes=NULL, cores=NULL){
  # Loop through combinations of samples and gene sets
  result_df = parallel::mclapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[[1]][i])
    geneset_tmp = as.vector(combo[[2]][i])
    pw_genes_tmp = pw_genes[[i]]
    # Test if genes in data set are enough
    if(length(pw_genes_tmp) >= min_genes){
      # Calculate average expression of genes within pathway for each spot/cell
      expr_subset = x[[sample_tmp]][pw_genes_tmp, ]
      pw_avg_exp = DelayedMatrixStats::colMeans2(expr_subset)
    } else{
      pw_avg_exp = rep(NA, ncol(x[[sample_tmp]]))
      names(pw_avg_exp) = colnames(x[[sample_tmp]])
    }
    pw_avg_exp = data.frame(as.list(pw_avg_exp), check.names=F)
    rownames(pw_avg_exp) = geneset_tmp

    return(pw_avg_exp)
  }, mc.cores=cores)
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
      pw_gsva_exp = GSVA::gsva(gsvapar, BPPARAM=BiocParallel::MulticoreParam(workers=cores), verbose=verbose)
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


##
# calculate_sum_distances
# The function calculates the sum of all pairwise distances in a distance matrix
# Returns a numeric value, the sum of distances
#
calculate_sum_distances = function(dists_mtx=NULL, select_spots=NULL){
  distances_high_spots = DelayedArray::DelayedArray(dists_mtx[select_spots, select_spots])
  upper_mtx_mask = DelayedArray::DelayedArray(upper.tri(distances_high_spots)) # Need to be DelayedArray to subset DelayedArray
  #upper_mtx_mask = DelayedArray::realize(upper.tri(distances_high_spots), BACKEND='RleArray') # Need to be DelayedArray to subset DelayedArray
  distances_high_spots[upper_mtx_mask] = 0 # Make upper half zero to avoid sum of distances twice
  #sum_high_distances = sum(distances_high_spots)
  sum_high_distances = DelayedMatrixStats::colSums2(distances_high_spots)
  sum_high_distances = sum(sum_high_distances)

  return(sum_high_distances)
}

