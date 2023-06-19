##
#' @title STgradient
#' @description 
#' @details XXXXXXX
#' XXXXXXX
#'
#' @param x an STlist with transformed gene expression
#' @param samples 
#' @param topgenes 
#' @param annot 
#' @param ref 
#' @param exclude 
#' @param out_rm 
#' @param limit
#' @param distsumm
#' @param min_nb
#' @param robust
#' @param cores the number of cores used during parallelization. If NULL (default),
#' the number of cores is defined automatically
#' @return a list of data frames with the results of the test
#'
#' @export
#'
#' @importFrom magrittr %>%
#
#
STgradient = function(x=NULL, samples=NULL, topgenes=2000, annot=NULL, ref=NULL, exclude=NULL, 
                      out_rm=F, limit=NULL, distsumm='min', min_nb=3, robust=T, cores=NULL){
  
  # x=normalized_stlist
  # samples=NULL
  # topgenes=100
  # annot="stclust_spw0_dspl1.5"
  # ref=2
  # exclude=NULL
  # out_rm=F
  # limit=NULL
  # distsumm='min'
  # min_nb=3
  # robust=T
  # cores=3
  
  # Make sure the reference cluster is character
  ref = as.character(ref)
  
  # Get sample names in case no specific sample was selected
  if(is.null(samples)){
    samplenames = names(x@tr_counts)
  }
  
  # Make sure that the requested annotation is present in all requested samples
  sample_rm = c()
  for(i in samplenames){
    if( !(annot %in% colnames(x@spatial_meta[[i]])) ){
      sample_rm = append(sample_rm, i)
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
    cores = count_cores(length(names(gene_sets)))
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
    nbs = colSums(dists_ref_tmp >= min_sample * 0.75 & dists_ref_tmp <= min_sample * 1.25)
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
          mutate(dist2ref=case_when(dist2ref <= dist2refupper ~ as.numeric(dist2ref)))
        
        rm(dist2reflower, dist2refupper) # Clean environment
      }
    }
    
    # Get expression from variable genes
    # Genes are identified within the range limit
    # Extract expression data (non-transformed counts to be passed to FindVariableFeatures)
    if(nrow(dists_summ_tmp) > 1){
      raw_cts = expandSparse(x@counts[[samplenames[i]]])
      
      # Get spots that have at least 1 distance value
      # However, if only one spot, then sample will be removed from analysis as cannot detect variable genes from single spot
      raw_cts = raw_cts[, dists_summ_tmp[['barcode']][ !is.na(dists_summ_tmp[['dist2ref']]) ], drop=F]
      
      # Number of rows larger than 1, because cannot compute variable genes with a single non-reference spot
      # Variable genes in minimum distance range
      if(ncol(raw_cts) > 1){
        vargenes = Seurat::FindVariableFeatures(raw_cts, verbose=F) %>%
          tibble::rownames_to_column(var='gene') %>%
          dplyr::arrange(desc(vst.variance.standardized)) %>%
          dplyr::select(gene) %>%
          unlist() %>%
          as.vector()
        vargenes = vargenes[1:topgenes] # Get number of genes defined by user
      }
    }
    
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
        quarts = quantile(dfdist2ref[[gene]], probs=c(0.25, 0.75))
        # Calculate inter-quartile range
        iqr_dist2ref = IQR(dfdist2ref[[gene]])
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
    dist_cor = tibble::tibble(gene=character(),
                              lm_coef=numeric(),
                              lm_pval=numeric(),
                              spearman_r=numeric(),
                              spearman_r_pval=numeric())
    
    # CORRELATIONS DISTANCE TO REFERENCE CLUSTER
    genes_sample = colnames(vargenes_expr %>% dplyr::select(-c('ypos', 'xpos', 'dist2ref')))
    for(gene in genes_sample){
      tibble_tmp = tibble::tibble(gene=character(),
                                  lm_coef=numeric(),
                                  lm_pval=numeric(),
                                  spearman_r=numeric(),
                                  spearman_r_pval=numeric())
      
      df_gene = vargenes_expr %>% dplyr::select(dist2ref, !!!gene)
      
      lm_res = list(estimate=NA, estimate_p=NA) 
      cor_res = list(estimate=NA, p.value=NA) 
      if(out_rm & !robust){ # Regular linear models after removal of outliers
        # Remove outliers
        df_gene_outrm = df_gene %>% 
          filter(!(rownames(df_gene) %in% outs_dist2ref[[gene]]))
        if(nrow(df_gene_outrm) > 1){
          # Run linear model and get summary
          lm_tmp = lm(df_gene_outrm[[gene]] ~ df_gene_outrm[['dist2ref']])
          lm_summ_tmp = summary(lm_tmp)[['coefficients']]
          if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
            lm_res = list(estimate=lm_summ_tmp[2,1],
                          estimate_p=lm_summ_tmp[2,4])
          }
          # Calculate Spearman correlation
          cor_res = cor.test(df_gene_outrm[['dist2ref']], df_gene_outrm[[gene]], method='spearman')
        }
      } else {
        if(robust){ # Robust linear models?
          df_gene_range = df_gene
          if(nrow(df_gene_range) > 1){
            # Run robust linear model and get summary
            lm_tmp = MASS::rlm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']], maxit=100)
            if(lm_tmp[['converged']] & lm_tmp[['coefficients']][2] != 0){ # Check the model converged and an effect was estimated
              # Run Wald test (MASS::rlm does not provide a p-value)
              lm_test_tmp = sfsmisc::f.robftest(lm_tmp)
              lm_res = list(estimate=summary(lm_tmp)[['coefficients']][2,1],
                            estimate_p=lm_test_tmp[['p.value']])
              # Calculate Spearman correlation
              cor_res= cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')
            }
          }
        } else{ # Regular linear models without outlier removal
          df_gene_range = df_gene
          if(nrow(df_gene_range) > 1){
            lm_tmp = lm(df_gene_range[[gene]] ~ df_gene_range[['dist2ref']])
            lm_summ_tmp = summary(lm_tmp)[['coefficients']]
            if(nrow(lm_summ_tmp) > 1){ # Test a linear model could be run
              lm_res = list(estimate=lm_summ_tmp[2,1],
                            estimate_p=lm_summ_tmp[2,4])
            }
            cor_res = cor.test(df_gene_range[['dist2ref']], df_gene_range[[gene]], method='spearman')
          }
        }
      }
      
      # Create row with results
      tibble_tmp = tibble::tibble(gene=gene, 
                                  lm_coef=lm_res[['estimate']],
                                  lm_pval=lm_res[['estimate_p']],
                                  spearman_r=cor_res[['estimate']],
                                  spearman_r_pval=cor_res[['p.value']])
      
      rm(list=grep("lm_|_res|_test|cor_|df_gene", ls(), value=T)) # Clean environment
      
      # Add row to result table if there is one row with results
      if(nrow(tibble_tmp) == 1){
        dist_cor = dplyr::bind_rows(dist_cor, tibble_tmp)
        rm(tibble_tmp) # Clean environment
      }
    }
    rm(genes_sample) # Clean environment
    
    # Adjust p-values for multiple comparison
    dist_cor[['adjpval_dist2ref']] = p.adjust(dist_cor[['dist2ref_spearman_p']], method='BH')
    dist_cor = dist_cor %>% 
      dplyr::relocate(adjpval_dist2ref, .after=dist2ref_spearman_p) %>%
      dplyr::arrange(adjpval_dist2ref)
    
    # Rename columns
    colnames(dist_cor) = c('gene', paste0(distsumm, '_', colnames(dist_cor[, -1])))
    
  }, mc.cores=cores)
  
  names(results_ls) = samplenames
  
  return(results_ls)
}

