##
#' @title STgradient: Tests of gene expression spatial gradients
#' @description Calculates Spearman's coefficients to detect genes showing expression spatial gradients
#' @details
#' The `STgradient` function fits linear models and calculates Spearman coefficients
#' between the expression of a gene and the minimum or average distance of spots or
#' cells to a reference tissue domain. In other wordsm the `STgradient` function
#' can be used to investigate if a gene is expressed higher in spots/cells closer to
#' a specific reference tissue domain, compared to spots/cells farther from the
#' reference domain (or viceversa as indicated by the Spearman's cofficient).
#'
#' @param x an STlist with transformed gene expression
#' @param samples the samples on which the test should be executed
#' @param topgenes the number of high-variance genes to be tested. These genes are
#' selected in descending order of variance as caclulated using Seurat's vst method
#' @param annot the name of a column in `@spatial_meta` containing the tissue domain
#' assignments for each spot or cell. These assignments can be generated using the
#' `STclust` function
#' @param ref one of the tissue domains in the column specified in `annot`,
#' corresponding to the "reference" cluster or domain. Spearman's correlations will
#' be calculated using spots assigned to domains other than this reference domain
#' (or domains specified in `exclude`).
#' @param exclude optional, a cluster/domain to exclude from the analysis
#' @param out_rm logical (optional), remove gene expression outliers defined by
#' the interquartile method. This option is only valid when `robust=F`
#' @param limit limite the analysis to spots/cells with distances to `ref` shorther
#' than the value specified here. Useful when gradients might occur at smaller scales
#' or when the domain in `ref` is scattered through the tissue. Caution must be used
#' due to difficult interpretation of imposed limits. It is suggested to run analysis
#' without restricted distances in addition for comparison.
#' @param distsumm the distance summary metric to use in correlations. One of `min` or `avg`
#' @param min_nb the minimum number of immediate neighbors a spot or cell has to
#' have in order to be included in the analysis. This parameter seeks to reduce the
#' effect of isolated `ref` spots on the correlation
#' @param robust logical, whether to use robust regression (`MASS` and `sfsmisc` packages)
#' @param nb_dist_thr a numeric vector of length two indicating the tolerance interval to assign
#' spots/cells to neighborhoods. The wider the range of the interval, the more likely
#' distinct neighbors to be considered. If NULL, `c(0.75, 1.25)` and `c(0.25, 3)` is assigned
#' for Visium and CosMx respectively.
#' @param log_dist logical, whether to apply the natural logarithm to the spot/cell
#' distances. It applies to all distances a constant (1e-200) to avoid log(0)
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats IQR lm cor.test p.adjust quantile
#
#
STgradient = function(x=NULL, samples=NULL, topgenes=2000, annot=NULL, ref=NULL, exclude=NULL,
                      out_rm=F, limit=NULL, distsumm='min', min_nb=3, robust=T, nb_dist_thr=NULL, log_dist=F, cores=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()

  # Make sure the reference cluster is character
  ref = as.character(ref)

  # Define samples using names (convert indexes to names if necessary)
  if(is.null(samples)){
    samplenames = names(x@tr_counts)
  } else{
    if(is.numeric(samples)){
      samplenames = as.vector(na.omit(names(x@tr_counts)[samples]))
    } else{
      samplenames = samples[samples %in% names(x@tr_counts)]
    }
    # Verify that sample names exist
    if(length(samplenames) == 0 | !any(samplenames %in% names(x@tr_counts))){
      raise_err(err_code="error0041")
    }
  }

  # Remove samples for which the requested annotation is not present
  sample_rm = c()
  for(i in samplenames){
    if( !(annot %in% colnames(x@spatial_meta[[i]])) ){
      sample_rm = append(sample_rm, i)
      cat(paste0('The annotation specified in `annot` is not present in ', i, '\n'))
    }
  }
  samplenames = samplenames[ !(samplenames %in% sample_rm) ]
  rm(sample_rm) # Clean env

  # Remove samples that do not have the requested reference cluster
  sample_rm = c()
  for(i in samplenames){
    if(sum(x@spatial_meta[[i]][[annot]] == ref) < 1){
      sample_rm = append(sample_rm, i)
    }
  }
  samplenames = samplenames[ !(samplenames %in% sample_rm) ]
  rm(sample_rm) # Clean env

  # Define number of cores to use
  if(is.null(cores)){
    cores = count_cores(length(samplenames))
  }

  # Define neighborhood tolerance
  if(is.null(nb_dist_thr) | !is.numeric(nb_dist_thr) | length(nb_dist_thr) != 2){
    nb_dist_thr = c(0.75, 1.25)
    if(x@misc[['platform']] != 'visium'){
      nb_dist_thr = c(0.25, 3)
    }
  }

  results_ls = parallel::mclapply(seq(samplenames), function(i){
    # Calculate euclidean distances
    coords_tmp = x@spatial_meta[[samplenames[i]]] %>%
      dplyr::select(c('libname', 'ypos', 'xpos')) %>% tibble::column_to_rownames('libname')
    dist_tmp = as.matrix(stats::dist(coords_tmp, method='euclidean'))
    rm(coords_tmp) # Clean env

    # Save spots in the different categories (ref, nonref, excl)
    ref_tmp = x@spatial_meta[[samplenames[i]]][['libname']][x@spatial_meta[[samplenames[i]]][[annot]] == ref]
    #excl_tmp = x@spatial_meta[[samplenames[i]]][['libname']][x@spatial_meta[[samplenames[i]]][[annot]] == exclude]
    nonref_tmp = x@spatial_meta[[samplenames[i]]][['libname']][!(x@spatial_meta[[samplenames[i]]][[annot]] %in% c(ref, exclude))]

    # Identify spots to be removed from reference if not enough neighbors
    # Get minimum distance among all spots within a sample (for Visium would be approximately the same for any sample)
    min_sample = min(as.data.frame(dist_tmp[lower.tri(dist_tmp)]))
    # Get distances among reference spots
    dists_ref_tmp = dist_tmp[ref_tmp, ref_tmp, drop=F]

    # Get number of neighbors within minimum distance
    # NOTE: When dealing with other technologies like SMI, will need to be more flexible with
    # minimum distances as not an array of equally distant spots. In this case, allowed a "buffer"
    # of a quarter of the minimum distance
    #if(x@misc[['platform']] == 'cosmx'){
      nbs = colSums(dists_ref_tmp >= min_sample * nb_dist_thr[1] & dists_ref_tmp <= min_sample * nb_dist_thr[2])
    #} else{
    #  nbs = colSums(dists_ref_tmp >= min_sample * 0.75 & dists_ref_tmp <= min_sample * 1.25)
    #}
    if(sum(nbs >= min_nb) < 1){ # At least 1 cluster of spots to continue with analysis
      nbs_keep = c()
    } else{
      nbs_keep = names(nbs)[nbs >= min_nb] # Save spots to be kept (enough neighbors)
    }
    rm(nbs, dists_ref_tmp) # Clean environment

    # Get summarized distances from the reference for each spot in the non-reference
    # Select spots in analysis (non reference in rows, reference in columns)
    dists_nonref_tmp = as.data.frame(dist_tmp[nonref_tmp, ref_tmp, drop=F])
    # Remove columns corresponding to spots without enough neighbors
    dists_nonref_tmp = dists_nonref_tmp[, colnames(dists_nonref_tmp) %in% nbs_keep, drop=F]

    # Check that distances are available for the comparison
    # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
    if(nrow(dists_nonref_tmp) > 1 & ncol(dists_nonref_tmp) > 0){
      if(distsumm == 'avg'){ # Summarize non-reference spots using minimun or mean distance
        dists_summ_tmp = tibble::tibble(barcode=rownames(dists_nonref_tmp),
                                        dist2ref=apply(dists_nonref_tmp, 1, mean))
      } else{
        dists_summ_tmp = tibble::tibble(barcode=rownames(dists_nonref_tmp),
                                        dist2ref=apply(dists_nonref_tmp, 1, min))
      }
    } else{
      dists_summ_tmp = tibble::tibble()
    }
    rm(dists_nonref_tmp) # Clean environment

    # Remove distances if outside user-specified limit
    if(!is.null(limit) & nrow(dists_summ_tmp) > 1){
      # Get lower and upper distance limits
      # If lower limit is higher than user limit, then set lower limit as upper limit
      if(!all(is.na(dists_summ_tmp[['dist2ref']]))){
        dist2reflower = min(dists_summ_tmp[['dist2ref']], na.rm=T)
        if(dist2reflower > limit){
          dist2refupper = dist2reflower
        } else{
          dist2refupper = limit
        }
        # Make NA the distances outside range
        dists_summ_tmp = dists_summ_tmp %>%
          dplyr::mutate(dist2ref=dplyr::case_when(dist2ref <= dist2refupper ~ as.numeric(dist2ref)))

        rm(dist2reflower, dist2refupper) # Clean environment
      }
    }

    # Get expression from variable genes
    # Genes are identified within the range limit
    # Extract expression data (non-transformed counts to be passed to FindVariableFeatures)
    dist_cor = tibble::tibble() # Initialize result data frame in case no computations can be done (e.g., sample without non-ref spots/cells)
    if(nrow(dists_summ_tmp) > 1){
      #raw_cts = expandSparse(x@counts[[samplenames[i]]])
      raw_cts = x@counts[[samplenames[i]]]
      # Get spots that have at least 1 distance value
      # However, if only one spot, then sample will be removed from analysis as cannot detect variable genes from single spot
      raw_cts = raw_cts[, dists_summ_tmp[['barcode']][ !is.na(dists_summ_tmp[['dist2ref']]) ], drop=F]

      # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
      # Variable genes in minimum distance range
      if(ncol(raw_cts) > 1){
        #vargenes = Seurat::FindVariableFeatures(raw_cts, verbose=F) %>%
        vargenes = calculate_vst(x=raw_cts) %>%
          dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
          dplyr::select(gene) %>%
          unlist() %>%
          as.vector()
        vargenes = vargenes[1:topgenes] # Get number of genes defined by user

        # Get transformed gene expression data (will be used for the correlations with distance)
        # Matrices will contain only the non-reference spots (as defined by non-NA value in distance)
        if(length(vargenes) > 0){
          vargenes_expr = expandSparse(x@tr_counts[[samplenames[i]]])
          vargenes_expr = vargenes_expr[(rownames(vargenes_expr) %in% vargenes), , drop=F] %>%
            t() %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var='barcode') %>%
            dplyr::left_join(dists_summ_tmp, ., by='barcode') %>%
            dplyr::filter(barcode %in% colnames(vargenes_expr)) %>%
            dplyr::right_join(x@spatial_meta[[samplenames[i]]] %>%
                                tibble::rownames_to_column(var='barcode') %>%
                                dplyr::select(barcode=libname, ypos, xpos), ., by='barcode') %>%
            tibble::column_to_rownames('barcode')

          rm(vargenes) # Clean environment
        } else{
          vargenes_expr = tibble::tibble()
        }

        # Detect gene expression outlier spots for each sample and gene
        if(out_rm & !robust){
          outs_dist2ref = list()
          dfdist2ref = vargenes_expr[!is.na(vargenes_expr[['dist2ref']]), ] %>%
            dplyr::select(-c('ypos', 'xpos', 'dist2ref'))

          for(gene in colnames(dfdist2ref)){
            # Calculate gene expression quartiles
            quarts = stats::quantile(dfdist2ref[[gene]], probs=c(0.25, 0.75))
            # Calculate inter-quartile range
            iqr_dist2ref = stats::IQR(dfdist2ref[[gene]])
            # Calculate distribution lower and upper limits
            low_up_limits = c((quarts[1]-1.5*iqr_dist2ref),
                              (quarts[2]+1.5*iqr_dist2ref))

            # Save outliers (barcodes)
            outs_dist2ref[[gene]] = rownames(dfdist2ref)[ dfdist2ref[[gene]] < low_up_limits[1] | dfdist2ref[[gene]] > low_up_limits[2] ]
          }
          rm(list=grep("iqr|quarts|low_up|dfdist2ref", ls(), value=T)) # Clean environment
        }

        # Calculate Spearman's correlations
        # Initialize data frame to store results
        dist_cor = tibble::tibble(sample_name=character(),
                                  gene=character(),
                                  lm_coef=numeric(),
                                  lm_pval=numeric(),
                                  spearman_r=numeric(),
                                  spearman_r_pval=numeric(),
                                  pval_comment=character())

        # CORRELATIONS DISTANCE TO REFERENCE CLUSTER
        genes_sample = colnames(vargenes_expr %>% dplyr::select(-c('ypos', 'xpos', 'dist2ref')))
        for(gene in genes_sample){
          tibble_tmp = tibble::tibble(sample_name=character(),
                                      gene=character(),
                                      lm_coef=numeric(),
                                      lm_pval=numeric(),
                                      spearman_r=numeric(),
                                      spearman_r_pval=numeric(),
                                      pval_comment=character())

          df_gene = vargenes_expr %>% dplyr::select(dist2ref, !!!gene)

          lm_res = list(estimate=NA, estimate_p=NA)
          cor_res = list(estimate=NA, p.value=NA)
          if(out_rm & !robust){ # Regular linear models after removal of outliers
            # Remove outliers
            if(length(outs_dist2ref[[gene]]) > 0){
              df_gene_outrm = df_gene %>%
                dplyr::filter( !(rownames(.) %in% outs_dist2ref[[gene]]) )
            } else{
              df_gene_outrm = df_gene
            }
            if(nrow(df_gene_outrm) > 1){
              # log-transform distances if selected by user
              if(log_dist){
                df_gene_outrm[['dist2ref']] = log(df_gene_outrm[['dist2ref']] + 1e-200)
              }

              # Run linear model and get summary
              lm_tmp = lm(df_gene_outrm[[gene]] ~ df_gene_outrm[['dist2ref']])
              lm_summ_tmp = summary(lm_tmp)[['coefficients']]
              if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
                lm_res = list(estimate=lm_summ_tmp[2,1],
                              estimate_p=lm_summ_tmp[2,4])
              }
              # Calculate Spearman correlation
              cor_res = tryCatch({cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman')}, warning=function(w){return(w)})
              pval_warn = NA_character_
              if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
                if(grepl('standard deviation is zero', cor_res$message)){
                  pval_warn = 'zero_st_deviation'
                } #else if(grepl('Cannot compute exact p-value with ties', cor_res$message)){ ## WARNING REMOVED AS MOST GENES WILL HAVE TIES = NON EXACT P-VAL
                  #pval_warn = 'non_exact_pvalue'
                #}
                cor_res = cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman', exact=F)
              }
            }
          } else {
            if(robust){ # Robust linear models?
              df_gene_range = df_gene
              if(nrow(df_gene_range) > 1){
                pval_warn = NA_character_

                # log-transform distances if selected by user
                if(log_dist){
                  df_gene_range[['dist2ref']] = log(df_gene_range[['dist2ref']] + 1e-200)
                }

                # Run robust linear model and get summary
                lm_tmp = MASS::rlm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']], maxit=100)
                if(lm_tmp[['converged']] & lm_tmp[['coefficients']][2] != 0){ # Check the model converged and an effect was estimated
                  # Run Wald test (MASS::rlm does not provide a p-value)
                  lm_test_tmp = sfsmisc::f.robftest(lm_tmp)
                  lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                                estimate_p=lm_test_tmp[['p.value']])
                  # Calculate Spearman correlation
                  cor_res = tryCatch({cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')}, warning=function(w){return(w)})
                  if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
                    if(grepl('standard deviation is zero', cor_res$message)){
                      pval_warn = 'zero_st_deviation'
                    } #else if(grepl('Cannot compute exact p-value with ties', cor_res$message)){ ## WARNING REMOVED AS MOST GENES WILL HAVE TIES = NON EXACT P-VAL
                      #pval_warn = 'non_exact_pvalue'
                    #}
                    cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman', exact=F)
                  }
                } else{
                  pval_warn = 'rob_regr_no_convergence'
                }
              }
            } else{ # Regular linear models without outlier removal
              df_gene_range = df_gene
              if(nrow(df_gene_range) > 1){

                # log-transform distances if selected by user
                if(log_dist){
                  df_gene_range[['dist2ref']] = log(df_gene_range[['dist2ref']] + 1e-200)
                }

                lm_tmp = lm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']])
                lm_summ_tmp = summary(lm_tmp)[['coefficients']]
                if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
                  lm_res = list(estimate=lm_summ_tmp[2,1],
                                estimate_p=lm_summ_tmp[2,4])
                }
                # Calculate Spearman correlation
                cor_res = tryCatch({cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')}, warning=function(w){return(w)})
                pval_warn = NA_character_
                if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
                  if(grepl('standard deviation is zero', cor_res$message)){
                    pval_warn = 'zero_st_deviation'
                  } #else if(grepl('Cannot compute exact p-value with ties', cor_res$message)){ ## WARNING REMOVED AS MOST GENES WILL HAVE TIES = NON EXACT P-VAL
                    #pval_warn = 'non_exact_pvalue'
                  #}
                  cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman', exact=F)
                }
              }
            }
          }

          # Create row with results
          tibble_tmp = tibble::tibble(sample_name=samplenames[i],
                                      gene=gene,
                                      lm_coef=lm_res[['estimate']],
                                      lm_pval=lm_res[['estimate_p']],
                                      spearman_r=as.vector(cor_res[['estimate']]),
                                      spearman_r_pval=cor_res[['p.value']],
                                      pval_comment=pval_warn)

          rm(list=grep("lm_|_res|_test|cor_|df_gene|exact_p", ls(), value=T)) # Clean environment

          # Add row to result table if there is one row with results
          if(nrow(tibble_tmp) == 1){
            dist_cor = dplyr::bind_rows(dist_cor, tibble_tmp)
            rm(tibble_tmp) # Clean environment
          }
        }
        rm(genes_sample) # Clean environment
      }
    }

    if(nrow(dist_cor) > 0){
      # Adjust p-values for multiple comparison
      dist_cor[['spearman_r_pval_adj']] = p.adjust(dist_cor[['spearman_r_pval']], method='BH')
      dist_cor = dist_cor %>%
        dplyr::relocate(spearman_r_pval_adj, .after=spearman_r_pval) %>%
        dplyr::arrange(spearman_r_pval_adj)

      # Rename columns
      colnames(dist_cor) = c('sample_name', 'gene', paste0(distsumm, '_', colnames(dist_cor[, -c(1,2)])))
    }

    return(dist_cor)
  }, mc.cores=cores)

  names(results_ls) = samplenames

  sample_rm = c()
  for(i in names(results_ls)){
    if(nrow(results_ls[[i]]) == 0){
      sample_rm = append(sample_rm, i)
    }
  }
  if(length(sample_rm) > 0){
    results_ls = results_ls[ !(names(results_ls) %in% sample_rm) ]
  }

  # Print time
  verbose = 1L
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose > 0L){
    cat(paste0('STgradient completed in ', round(end_t, 2), ' min.\n'))
  }

  return(results_ls)
}

