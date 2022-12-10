##
#' @title STenrich
#' @description Test for spatial enrichment of gene expression sets in ST data sets
#' @details The function performs a randomization test to assess if the sum of
#' distances between cells/spots with high expression of a gene set is lower than
#' the sum of distances of randomly selected cells/spots. The cells/spots are
#' considered as having high gene set expression if the average expression of genes in a
#' set is higher than the average expression plus a number of standard deviations
#' defined by `min_sd`. Control over the size of regions with high expression is
#' provided by setting the minimum number of cells or spots (`min_units`). This method
#' is a modification of the method devised by Hunter et al. 2021 (zebrafish melanoma study)
#'
#' @param x an STlist with transformed gene expression
#' @param gene_sets a named list of gene sets to test
#' @param reps the numbe rof random samples to extract
#' @param min_sd the number of standard deviations to set the minimum gene set expression threshold
#' @param min_units Minimum number of spots with high expression of a pathway for
#' that gene set to be considered in the analysis
#' @param min_genes the minimum number of genes of a gene set present in the data set
#' for that gene set to be included
#' @param pval_adj_method the method for multiple comparison adjustment of p-values. Options
#' passed to `p.adjust`
#' @param seed the seed number for the selection of random samples
#' @param cores the number of cores used during parallelization
#' @return a list of data frames with the results of the test
#'
#' @importFrom magrittr %>%
#
STenrich = function(x=NULL, gene_sets=NULL, reps=1000, min_sd=1, min_units=20, min_genes=5, pval_adj_method='BH', seed=12345, cores=NULL){
  cat(crayon::yellow("Running STenrich...\n"))
  set.seed(seed)
  # Loop through samples in STlist
  result_dfs = list()
  for(i in names(x@tr_counts)){
    cat(crayon::yellow(paste0('\tSample: ', i, '...\n')))
    # Extract spots to be used in analysis
    # This selection implemented proactively as analysis might later be applied to tissue nuches within samples
    tissue_spots = x@spatial_meta[[i]][['libname']]

    # Extract gene expression
    exp = expandSparse(x@tr_counts[[i]])
    exp = exp[, tissue_spots]

    # Extract spot coordinates and match order of spots
    coords_df = x@spatial_meta[[i]][, c('libname', 'xpos', 'ypos')]
    coords_df = coords_df[match(colnames(exp), coords_df[['libname']]), ]

    # Define number of cores to use
    if(is.null(cores)){
      cores = count_cores(length(names(gene_sets)))
    }
    # Perform tests in parallel
    pval_ls = parallel::mclapply(seq(names(gene_sets)), function(pw){
      pw_genes = gene_sets[[pw]]
      pw_genes = pw_genes[pw_genes %in% rownames(exp)] # Filter genes that aren't in the expression matrix.

      # Test if genes in data set are enough
      if(length(pw_genes) >= min_genes){
        # Average expression of genes within pathway for each spot
        pw_avg_exp = apply(exp[pw_genes, ], 2, mean)

        # Find spots that highly express this pathway (mean + stdev_t*std in this case)
        avg = mean(pw_avg_exp)
        std = sd(pw_avg_exp)
        exp_thresh = avg + (min_sd*std)

        # Extract spots with average expression above sd threshold
        high_spots_bc = names(which(pw_avg_exp >= exp_thresh))

        # Are there at least X number of cells/spots?
        if(length(high_spots_bc) >= min_units){

          # Compute the distance between these spots' XY coordinates with high expression
          coords_high_spots = coords_df[coords_df[['libname']] %in% high_spots_bc, ]
          distances_high_spots = as.matrix(rdist::rdist(coords_high_spots[, c('xpos', 'ypos')], metric='euclidean')) # Distance computation
          distances_high_spots[upper.tri(distances_high_spots)] = 0 # Make upper half zero to avoid sum of distances twice
          sum_high_distances = sum(distances_high_spots)

          # Compute random distance permutations
          sum_rand_distances = c()
          for(rep in 1:reps){
            rand_idx = sample(x=1:nrow(coords_df), size=length(high_spots_bc))
            rand_coord = coords_df[rand_idx, ]
            rand_dist = as.matrix(rdist::rdist(rand_coord[, c('xpos', 'ypos')], metric='euclidean')) # Distance computation for random spots
            rand_dist[upper.tri(rand_dist)] = 0 # Make upper half zero to avoid sum of distances twice
            sum_rand_distances = append(sum_rand_distances, sum(rand_dist))
          }

          # Compute p-value
          # Ho: sum of observed distances is larger than sum of random distances
          # Ha: sum of observed distances is smaller than sum of random distances
          count_test = sum(sum_high_distances > sum_rand_distances) # Times observed dist was higher than random dists
          p_val = count_test / reps

          # Get results in data frame form
          pval_tmp = tibble::tibble(gene_set=names(gene_sets)[pw],
                                    size_test=length(pw_genes),
                                    size_gene_set=length(gene_sets[[pw]]),
                                    p_value=p_val)

        } # Close test minimum spots
      } else{
        pval_tmp = tibble::tibble(gene_set=names(gene_sets)[pw],
                                  size_test=length(pw_genes),
                                  size_gene_set=length(gene_sets[[pw]]),
                                  p_value=NA)
      } # Close test minimum genes

      return(pval_tmp)
    }, mc.cores=cores) # Close mclappy

    # Compile results in a single data frame
    pval_df = dplyr::bind_rows(pval_ls) %>%
      tibble::add_column(sample_name=i, .before=1)

    # Adjust p-values
    pval_df[['adj_p_value']] = p.adjust(pval_df[['p_value']], method=pval_adj_method)
    pval_df = pval_df[order(pval_df[['adj_p_value']]), ]
    result_dfs[[i]] = pval_df
  }
  cat(crayon::green("STenrich completed."))
  return(result_dfs)
}

