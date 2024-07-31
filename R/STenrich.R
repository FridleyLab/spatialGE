##
#' @title STenrich
#' @description Test for spatial enrichment of gene expression sets in ST data sets
#' @details The function performs a randomization test to assess if the sum of
#' distances between cells/spots with high expression of a gene set is lower than
#' the sum of distances of randomly selected cells/spots. The cells/spots are
#' considered as having high gene set expression if the average expression of genes in a
#' set is higher than the average expression plus a `num_sds` times the standard deviation.
#' Control over the size of regions with high expression is provided by setting the
#' minimum number of cells or spots (`min_units`). This method is a modification of
#' the method devised by Hunter et al. 2021 (zebrafish melanoma study)
#'
#' @param x an STlist with transformed gene expression
#' @param samples a vector with sample names or indexes to run analysis
#' @param gene_sets a named list of gene sets to test. The names of the list should
#' identify the gene sets to be tested
#' @param reps the number of random samples to be extracted. Default is 1000 replicates
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
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#
#
STenrich = function(x=NULL, samples=NULL, gene_sets=NULL, reps=1000, num_sds=1, min_units=20, min_genes=5, pval_adj_method='BH', seed=12345, cores=NULL){
  # Record time
  zero_t = Sys.time()

  cat(crayon::yellow("Running STenrich...\n"))

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
      stop('None of the requested samples are present in the STlist.')
    }
  }

  # Define number of cores for parallelization of tests
  if(is.null(cores)){
    cores = count_cores(length(samples))
  } else{
    cores = ceiling(cores)
  }

  # Loop through samples in STlist
  result_dfs = parallel::mclapply(samples, function(i){
    system(sprintf('echo "%s"', crayon::yellow(paste0("\tSample: ", i, "..."))))
    # Extract spots to be used in analysis
    # This selection implemented proactively as analysis might later be applied to tissue niches within samples
    tissue_spots = x@spatial_meta[[i]][['libname']]

    # Extract gene expression
    exp = expandSparse(x@tr_counts[[i]])
    exp = exp[, tissue_spots]

    # Extract spot coordinates and match order of spots
    coords_df = x@spatial_meta[[i]][, c('libname', 'xpos', 'ypos')]
    coords_df = coords_df[match(colnames(exp), coords_df[['libname']]), ]
    coords_df = coords_df %>% tibble::column_to_rownames(var='libname')

    # Calculate distances for the sample
    distances_spots = as.matrix(stats::dist(coords_df[, c('xpos', 'ypos')], method='euclidean'))
    rm(coords_df) # Clean env

    # Perform tests in parallel
    pval_df = tibble::tibble()
    for(pw in 1:length(gene_sets)){
      pw_genes = gene_sets[[pw]]
      pw_genes = pw_genes[pw_genes %in% rownames(exp)] # Filter genes that aren't in the expression matrix.

      # Test if genes in data set are enough
      if(length(pw_genes) >= min_genes){
        # Average expression of genes within pathway for each spot
        pw_avg_exp = apply(exp[pw_genes, ], 2, mean)

        # Find spots that highly express this pathway (mean + stdev_t*std in this case)
        avg = mean(pw_avg_exp)
        std = sd(pw_avg_exp)
        exp_thresh = avg + (num_sds*std)

        # Extract spots with average expression above sd threshold
        high_spots_bc = names(which(pw_avg_exp >= exp_thresh))

        # Are there at least 'min_units' number of cells/spots?
        if(length(high_spots_bc) >= min_units){

          # Compute the distance between these spots' XY coordinates with high expression
          distances_high_spots = distances_spots[high_spots_bc, high_spots_bc]
          distances_high_spots[upper.tri(distances_high_spots)] = 0 # Make upper half zero to avoid sum of distances twice
          sum_high_distances = sum(distances_high_spots)

          rm(distances_high_spots) # Clean env

          # Compute random distance permutations
          set.seed(seed)
          sum_rand_distances = c()
          for(rep in 1:reps){
            rand_idx = sample(x=colnames(distances_spots), size=length(high_spots_bc))
            rand_dist = distances_spots[rand_idx, rand_idx]
            rand_dist[upper.tri(rand_dist)] = 0 # Make upper half zero to avoid sum of distances twice
            sum_rand_distances = append(sum_rand_distances, sum(rand_dist))

            rm(rand_idx, rand_dist) # Clean env
          }

          # Compute p-value
          # Ho: sum of observed distances is larger than sum of random distances
          # Ha: sum of observed distances is smaller than sum of random distances
          count_test = sum(sum_rand_distances < sum_high_distances) # Times observed dist was higher than random dists
          p_val = count_test / reps

          rm(high_spots_bc) # Clean env
        } else{
          p_val=NA
        } # Close test minimum spots
      } else{
        p_val=NA
      } # Close test minimum genes

      # Get results in data frame form
      pval_tmp = tibble::tibble(gene_set=names(gene_sets)[pw],
                                size_test=length(pw_genes),
                                size_gene_set=length(gene_sets[[pw]]),
                                p_value=p_val)

      # Compile results in a single data frame
      pval_df = dplyr::bind_rows(pval_df,
                                 pval_tmp %>%
                                   tibble::add_column(sample_name=i, .before=1))
    }

    # Adjust p-values
    pval_df[['adj_p_value']] = p.adjust(pval_df[['p_value']], method=pval_adj_method)
    pval_df = pval_df[order(pval_df[['adj_p_value']]), ]

    return(pval_df)
  }, mc.cores=cores)
  names(result_dfs) = samples

  # Print time
  verbose = 1L
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose > 0L){
    cat(crayon::green(paste0('STenrich completed in ', round(end_t, 2), ' min.\n')))
  }

  return(result_dfs)
}

