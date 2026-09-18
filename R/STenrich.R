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
  samples = resolve_samples(x, samples, slot='spatial_meta', numeric_handling='na_omit',
                            char_transform='exact', post_check='any_in',
                            err_code='error0020')

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
  user_cores = !is.null(cores)
  cores = resolve_cores(length(samples), cores, coerce='ceiling')

  # Extract spots to be used in analysis (in case analysis restricted to domain was requested)
    # Extract spot coordinates and match order of spots
  coords_df = lapply(samples, function(i){
    df_tmp = x@spatial_meta[[i]][, c('libname', 'xpos', 'ypos')]
    if(length(tissue_spots) > 0){
      df_tmp = df_tmp[df_tmp[['libname']] %in% tissue_spots[[i]], ]
    }
    # The row subset above leaves non-contiguous row names behind (e.g. 1, 5, 9),
    # which column_to_rownames() rejects. Drop them first so the domain-restricted
    # path behaves like the unrestricted one.
    df_tmp = df_tmp %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var='libname')
    df_tmp = as.matrix(df_tmp)

    return(df_tmp)
  })
  names(coords_df) = samples

  # Create data frame with combinations of samples and gene sets
  combo = expand.grid(sample_name=samples, gene_set=names(gene_sets))

  # Find genes from each gene set that are present in data
  pw_genes = lapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    pw_genes_tmp = unique(gene_sets[[gs_tmp]])
    pw_genes_tmp = intersect(pw_genes_tmp, sample_genes(x, sample_tmp, 'tr_counts'))

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
      cores = resolve_cores(nrow(combo), NULL, coerce='ceiling')
    }
    # Calculate average gene set expression
    result_df = calculate_gs_mean_exp(delayed_x, combo, pw_genes, min_genes, cores)

  } else if(score_type == 'gsva'){
    if(verbose){
      cat("\tCalculating GSVA score...\n")
    }
    # Define if user input number of cores, otherwise determine
    if(!user_cores){
      cores = resolve_cores(length(gene_sets), NULL, coerce='ceiling')
    }
    # Calculate GSVA scores
    result_df = calculate_gs_gsva_score(delayed_x, combo, gene_sets=gene_sets, pw_genes=pw_genes, min_genes=min_genes, cores=cores)
  }

  # Find spots that highly express this pathway (mean + stdev_t*std in this case), calculate observed sum distances, and perform permutations
  # Define if user input number of cores, otherwise determine
  if(!user_cores){
    cores = resolve_cores(nrow(combo), NULL, coerce='ceiling')
  }
  pval_res = sge_lapply(1:nrow(combo), function(i){
    sample_tmp = as.vector(combo[i, 1])
    gs_tmp = as.vector(combo[i, 2])
    expr_vals = unlist(result_df[[sample_tmp]][gs_tmp, ])
    stenrich_perm_one_combo(expr_vals, coords_df[[sample_tmp]], sample_tmp, gs_tmp,
                            length(pw_genes[[i]]), length(gene_sets[[gs_tmp]]),
                            num_sds, min_units, reps, seed, verbose)
  }, cores=cores)

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

