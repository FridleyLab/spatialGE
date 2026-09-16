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
#' @param verbose logical, whether to print text to console
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom stats IQR lm cor.test p.adjust quantile
#
#
STgradient = function(x=NULL, samples=NULL, topgenes=2000, annot=NULL, ref=NULL, exclude=NULL,
                      out_rm=FALSE, limit=NULL, distsumm='min', min_nb=3, robust=TRUE,
                      nb_dist_thr=NULL, log_dist=FALSE, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()

  # Make sure the reference cluster is character
  ref = as.character(ref)

  # Define samples using names (convert indexes to names if necessary)
  samplenames = resolve_samples(x, samples, slot='tr_counts', numeric_handling='na_omit',
                                char_transform='exact', post_check='any_in',
                                err_code='error0041')

  # Remove samples for which the requested annotation is not present
  sample_rm = c()
  for(i in samplenames){
    if( !(annot %in% colnames(x@spatial_meta[[i]])) ){
      sample_rm = append(sample_rm, i)
      message(paste0('The annotation specified in `annot` is not present in ', i, '\n'))
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
  cores = resolve_cores(length(samplenames), cores, coerce='none')

  # Define neighborhood tolerance
  if(is.null(nb_dist_thr) | !is.numeric(nb_dist_thr) | length(nb_dist_thr) != 2){
    nb_dist_thr = c(0.75, 1.25)
    if(x@misc[['platform']] != 'visium'){
      nb_dist_thr = c(0.25, 3)
    }
  }

  results_ls = sge_lapply(seq(samplenames), function(i){
    stgradient_one_sample(sample_counts(x, samplenames[i]), sample_tr_counts(x, samplenames[i]),
                          sample_spatial_meta(x, samplenames[i]), samplenames[i], annot, ref, exclude,
                          distsumm, nb_dist_thr, min_nb, limit, topgenes, out_rm, robust, log_dist)
  }, cores=cores)

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
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STgradient completed in ', round(end_t, 2), ' min.\n'))
  }

  return(results_ls)
}


# Helpers ----------------------------------------------------------------------

##
# @title stgradient_one_sample: Run the distance-to-reference gradient test
# for one sample
# @description The per-sample unit of work dispatched by `STgradient()`'s
# `sge_lapply` call, extracted so it is directly callable/benchmarkable for
# one sample without going through the whole STlist. Delegates to
# `stgradient_dists_one_sample()`, `stgradient_expr_one_sample()`, and
# `stgradient_fit_one_gene()` for the distance summarization, variable-gene
# expression extraction, and per-gene model fitting steps respectively.
#
# @param counts_mat the sample's raw count matrix.
# @param tr_counts_mat the sample's transformed count matrix.
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param sample_name the sample's name, stored in the result rows.
# @param annot,ref,exclude,distsumm,nb_dist_thr,min_nb,limit,topgenes,out_rm,robust,log_dist
# same as the corresponding `STgradient()` parameters.
# @return a tibble with one row per gene tested (possibly zero rows).
#
stgradient_one_sample = function(counts_mat, tr_counts_mat, spatial_meta_df, sample_name, annot, ref, exclude,
                                 distsumm, nb_dist_thr, min_nb, limit, topgenes, out_rm, robust, log_dist){

  dists_summ_tmp = stgradient_dists_one_sample(spatial_meta_df, annot, ref, exclude, distsumm, nb_dist_thr, min_nb, limit)

  # Get expression from variable genes
  # Genes are identified within the range limit
  # Extract expression data (non-transformed counts to be passed to FindVariableFeatures)
  dist_cor = tibble::tibble() # Initialize result data frame in case no computations can be done (e.g., sample without non-ref spots/cells)
  if(nrow(dists_summ_tmp) > 1){
    expr_res = stgradient_expr_one_sample(counts_mat, tr_counts_mat, spatial_meta_df, dists_summ_tmp, topgenes)
    vargenes_expr = expr_res[['vargenes_expr']]

    # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
    if(expr_res[['have_expr']]){

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
        rm(dfdist2ref, quarts, iqr_dist2ref, low_up_limits) # Clean environment
      } else{
        outs_dist2ref = list()
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
        df_gene = vargenes_expr %>% dplyr::select(dist2ref, !!!gene)

        fit_res = stgradient_fit_one_gene(df_gene, gene, out_rm, robust, log_dist, outs_dist2ref[[gene]])

        # Create row with results
        tibble_tmp = tibble::tibble(sample_name=sample_name,
                                    gene=gene,
                                    lm_coef=fit_res[['lm_res']][['estimate']],
                                    lm_pval=fit_res[['lm_res']][['estimate_p']],
                                    spearman_r=as.vector(fit_res[['cor_res']][['estimate']]),
                                    spearman_r_pval=fit_res[['cor_res']][['p.value']],
                                    pval_comment=fit_res[['pval_warn']])

        # Add row to result table if there is one row with results
        if(nrow(tibble_tmp) == 1){
          dist_cor = dplyr::bind_rows(dist_cor, tibble_tmp)
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
}


##
# @title stgradient_dists_one_sample: Summarize each non-reference spot's
# distance to the reference domain, for one sample
# @description The distance-summarization step of `stgradient_one_sample()`,
# extracted for readability and standalone benchmarking.
#
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param annot the column in `spatial_meta_df` with the tissue domain assignments.
# @param ref the reference tissue domain.
# @param exclude optional, a domain to exclude from the analysis.
# @param distsumm the distance summary metric, 'min' or 'avg'.
# @param nb_dist_thr a numeric vector of length two, the neighborhood tolerance interval.
# @param min_nb the minimum number of immediate neighbors required to keep a reference spot.
# @param limit optional, restrict distances to spots within this limit of the reference.
# @return a tibble with columns `barcode` and `dist2ref` (possibly zero rows).
#
stgradient_dists_one_sample = function(spatial_meta_df, annot, ref, exclude, distsumm, nb_dist_thr, min_nb, limit){

  # Calculate euclidean distances
  coords_tmp = spatial_meta_df %>%
    dplyr::select(c('libname', 'ypos', 'xpos')) %>% tibble::column_to_rownames('libname')
  dist_tmp = as.matrix(stats::dist(coords_tmp, method='euclidean'))
  rm(coords_tmp) # Clean env

  # Save spots in the different categories (ref, nonref, excl)
  ref_tmp = spatial_meta_df[['libname']][spatial_meta_df[[annot]] == ref]
  nonref_tmp = spatial_meta_df[['libname']][!(spatial_meta_df[[annot]] %in% c(ref, exclude))]

  # Identify spots to be removed from reference if not enough neighbors
  # Get minimum distance among all spots within a sample (for Visium would be approximately the same for any sample)
  min_sample = min(as.data.frame(dist_tmp[lower.tri(dist_tmp)]))
  # Get distances among reference spots
  dists_ref_tmp = dist_tmp[ref_tmp, ref_tmp, drop=F]

  # Get number of neighbors within minimum distance
  # NOTE: When dealing with other technologies like SMI, will need to be more flexible with
  # minimum distances as not an array of equally distant spots. In this case, allowed a "buffer"
  # of a quarter of the minimum distance
  nbs = colSums(dists_ref_tmp >= min_sample * nb_dist_thr[1] & dists_ref_tmp <= min_sample * nb_dist_thr[2])
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

  return(dists_summ_tmp)
}


##
# @title stgradient_expr_one_sample: Extract top variable-gene expression
# for the spots with a valid distance-to-reference, for one sample
# @description The variable-gene expression extraction step of
# `stgradient_one_sample()`, extracted for readability and standalone
# benchmarking.
#
# @param counts_mat the sample's raw count matrix.
# @param tr_counts_mat the sample's transformed count matrix.
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param dists_summ_tmp the tibble returned by `stgradient_dists_one_sample()`.
# @param topgenes the number of top variable genes to keep.
# @return a list with `vargenes_expr` (a tibble of gene expression joined with
# distance and coordinates, possibly empty) and `have_expr` (logical, whether
# more than one spot had usable counts).
#
stgradient_expr_one_sample = function(counts_mat, tr_counts_mat, spatial_meta_df, dists_summ_tmp, topgenes){

  # Get spots that have at least 1 distance value
  # However, if only one spot, then sample will be removed from analysis as cannot detect variable genes from single spot
  raw_cts = counts_mat
  raw_cts = raw_cts[, dists_summ_tmp[['barcode']][ !is.na(dists_summ_tmp[['dist2ref']]) ], drop=F]

  # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
  # Variable genes in minimum distance range
  vargenes_expr = tibble::tibble()
  if(ncol(raw_cts) > 1){
    vargenes = calculate_vst(x=raw_cts) %>%
      dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
      dplyr::select(gene) %>%
      unlist() %>%
      as.vector()
    vargenes = vargenes[1:topgenes] # Get number of genes defined by user

    # Get transformed gene expression data (will be used for the correlations with distance)
    # Matrices will contain only the non-reference spots (as defined by non-NA value in distance)
    if(length(vargenes) > 0){
      vargenes_mat = expandSparse(tr_counts_mat)
      vargenes_expr = vargenes_mat[(rownames(vargenes_mat) %in% vargenes), , drop=F] %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='barcode') %>%
        dplyr::left_join(dists_summ_tmp, ., by='barcode') %>%
        dplyr::filter(barcode %in% colnames(vargenes_mat)) %>%
        dplyr::right_join(spatial_meta_df %>%
                            tibble::rownames_to_column(var='barcode') %>%
                            dplyr::select(barcode=libname, ypos, xpos), ., by='barcode') %>%
        tibble::column_to_rownames('barcode')

      rm(vargenes) # Clean environment
    }
  }

  return(list(vargenes_expr=vargenes_expr, have_expr=(ncol(raw_cts) > 1)))
}


##
# @title stgradient_fit_one_gene: Fit the distance-vs-expression model for
# one gene
# @description The per-gene model-fitting step of `stgradient_one_sample()`,
# extracted so it is directly callable/benchmarkable for a single gene.
#
# @param df_gene a two-column tibble/data frame with `dist2ref` and the gene's
# expression column (named `gene`).
# @param gene the gene name (must match a column name in `df_gene`).
# @param out_rm logical, whether to remove expression outliers before fitting
# (only applied when `robust=FALSE`).
# @param robust logical, whether to use robust regression (`MASS::rlm`).
# @param log_dist logical, whether to log-transform distances before fitting.
# @param outliers a vector of barcodes (rownames of `df_gene`) to exclude when
# `out_rm & !robust`; ignored otherwise.
# @return a list with `lm_res` (list with `estimate`/`estimate_p`), `cor_res`
# (a `cor.test`/`simpleWarning` result), and `pval_warn`.
#
stgradient_fit_one_gene = function(df_gene, gene, out_rm, robust, log_dist, outliers=NULL){

  lm_res = list(estimate=NA, estimate_p=NA)
  cor_res = list(estimate=NA, p.value=NA)
  pval_warn = NA_character_

  if(out_rm & !robust){ # Regular linear models after removal of outliers
    # Remove outliers
    if(length(outliers) > 0){
      df_gene_outrm = df_gene %>%
        dplyr::filter( !(rownames(.) %in% outliers) )
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
      if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
        if(grepl('standard deviation is zero', cor_res$message)){
          pval_warn = 'zero_st_deviation'
        }
        cor_res = cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman', exact=F)
      }
    }
  } else {
    if(robust){ # Robust linear models?
      df_gene_range = df_gene
      if(nrow(df_gene_range) > 1){

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
            }
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
        if(any(class(cor_res) == 'simpleWarning')){ # Let known user if p-value could not be exactly calculated
          if(grepl('standard deviation is zero', cor_res$message)){
            pval_warn = 'zero_st_deviation'
          }
          cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman', exact=F)
        }
      }
    }
  }

  return(list(lm_res=lm_res, cor_res=cor_res, pval_warn=pval_warn))
}

