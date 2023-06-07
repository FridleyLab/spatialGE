##
#' @title STdiff: Differential gene expression analysis for spatial transcriptomics data
#' @description Tests for differentially expressed genes using linear models with or
#' without spatial covariance structures
#' @details
#' The method tests for differentially expressed genes between groups of spots/cells
#' (e.g., clusters) in a spatial transcriptomics sample. Specifically, the function
#' tests for genes with significantly higher or lower gene expression in one group of
#' spots/cells with respect to the rest of spots/cells in the sample. The method first
#' runs non-spatial linear models on the genes to detect differentially expressed genes.
#' Then spatial linear models with exponential covariance structure are fit on a
#' subset of genes detected as differentially expressed by the non-linear models (`sp_topgenes`).
#' If running on clusters detected via STclust, the user can specify the assignments
#' using the same parameters (`ws`, `ks`, `deepSplit`). Otherwise, the assignments are
#' specified by indicating the column name in `x@spatial_meta`. The function uses `spaMM::fitme` and
#' is computationally expensive even on HPC environments. To run the STdiff using the
#' non-spatial approach (faster), set `sp_topgenes=0`.
#'
#' @param x an STlist
#' @param samples an integer indicating the spatial samples to be included in the DE tests.
#' Numbers follow the order in `names(x@counts)`. Sample names are also allowed.
#' If NULL, performs tests on all samples
#' @param annot a column name in `x@spatial_meta` containing the groups/clusters to
#' be tested
#' @param ws the spatial weight used in STclust
#' @param ks the k value used in STclust, or `dtc` (default) for dynamicTreeCut clusters
#' @param deepSplit the deepSplit value if used in STclust
#' @param topgenes an integer indicating the top variable genes to select from each sample
#' based on variance (default=5000). If NULL, all genes are selected.
#' @param pval_adj Method to adjust p-values. Defaults to `FDR`. Other options as
#' available from `p.adjust`
#' @param pval_thr cut-off of adjusted p-values to define differentially expressed genes from
#' non-spatial linear models. A proportion of genes (`sp_topgenes`) under this cut-off
#' will be applied the spatial models. Default=0.05
#' @param sp_topgenes Proportion of differentially expressed genes from non-spatial
#' linear models (and controlled by `pval_thr`) to use in differential gene expression
#' analysis with spatial linear models. If 0 (zero), no spatial models are fit. Default=0.2
#' @param clusters cluster name(s) to test DE genes, as opposed to all clusters.
#' @param pairwise whether or not to carry tests on a pairwise manner. The default is
#' `pairwise=F`, meaning that DE genes are tested by comparing each cluster to the
#' rest of the pooled cell/spots.
#' @param verbose output progress indicators. If `verbose=0`, no text is shown in console.
#' Other values are 1 and 2 indicating increasing level of verbosity. Default is
#' `verbose=1`
#' @param cores Number of cores to use in parallelization. If `NULL`, the number of
#' cores to use is detected automatically
#' @return a data frame with results of differential gene expression analysis
#'
#' @export
#'
#' @importFrom magrittr %>%
#
#
STdiff = function(x=NULL, samples=NULL, annot=NULL, ws=NULL, ks='dtc', deepSplit=NULL,
                topgenes=5000, pval_thr=0.05, pval_adj='fdr', test_type='mm', sp_topgenes=0.2,
                clusters=NULL, pairwise=F, verbose=1L, cores=NULL){
  # Record time
  zero_t = Sys.time()

  topgenes = as.integer(topgenes)
  pval_thr = as.double(pval_thr)
  sp_topgenes = as.double(sp_topgenes)

  verbose = as.integer(verbose)
  if(!is.integer(verbose) | !(verbose %in% c(0L, 1L, 2L))){
    verbose = 1L
  }

  # Stop function if topgenes and sp_topgenes are not valid
  if((round(topgenes, 0) <= 0) | (sp_topgenes < 0 | sp_topgenes > 1)){
    stop('topgenes or sp_topgenes contain invalid values.')
  }

  # Define samples
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = names(x@spatial_meta)[samples]
    }
    if(length(grep(paste0(samples, collapse='|'), names(x@spatial_meta))) == 0){
      stop('The requested samples are not present in the STlist spatial metadata.')
    }
  }

  # Define column to test
  if(is.null(annot)){
    annot = grep('^stclust_spw', colnames(x@spatial_meta[[1]]), value=T)
    if(!is.null(ws)){
      if(any(ws == 0)){ # To avoid zero ('0') matching other weights
        ws_tmp = ws[ws != 0]
        annot_tmp = grep('stclust_spw0_|stclust_spw0$', annot, value=T)
        if(length(ws_tmp) > 0){ # In case there are other ws in addition to '0'
          annot = c(annot_tmp,
                    grep(paste0('stclust_spw', ws_tmp, collapse='|'), annot, value=T))
        } else{
          annot = annot_tmp
        }
        rm(ws_tmp, annot_tmp)
      } else{
        annot = grep(paste0('stclust_spw', ws, collapse='|'), annot, value=T)
      }
    }

    if(ks[1] != 'dtc'){
      annot = grep(paste0('_k', ks,'$', collapse='|'), annot, value=T)
    } else if(ks[1] == 'dtc'){
      if(!is.null(deepSplit)){
        annot = grep(paste0('_dspl', deepSplit, '$', collapse='|'), annot, value=T)
      } else{
        annot = grep('_dspl', annot, value=T)
      }
    } else{
      stop('Specify a k value to test, or use ks="dtc" (default).')
    }
  }

  # Check that the meta data column exists (if not in specific sample, remove sample)
  samples_tmp = samples
  for(i in samples){
    if(length(grep(paste0('^', annot, '$', collapse='|'),  colnames(x@spatial_meta[[i]]))) == 0){
      samples_tmp = samples_tmp[!(samples_tmp %in% i)]
    } else if(length(grep(paste0('^', annot, '$', collapse='|'),  colnames(x@spatial_meta[[i]]))) > 1){
      stop('Only one clustering solution or metadata column can be tested at a time.')
    }

    if(length(samples_tmp) == 0){
      stop('No samples left to test. Are the requested annotations/clusters present in at least one sample?')
    }
  }
  samples = samples_tmp
  rm(samples_tmp) # Clean env

  # Calculate standardized variance for all samples (regardless if sample is not in requested test)
  for(i in 1:length(x@tr_counts)){
    # Check that vst is not already calculated
    if(length(grep('vst.variance.standardized', colnames(x@gene_meta[[i]]))) == 0){
      x@gene_meta[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::select('gene', 'vst.variance.standardized') %>%
        dplyr::left_join(x@gene_meta[[i]], ., by='gene')
    }
  }

  # Print metadata assignment to be tested
  if(grepl('stclust_', annot[1])){
    spw_print = stringr::str_extract(annot[1], '(?<=spw)0\\.[0-9]+')
    if(grepl('_dspl', annot[1])){
      cut_print = stringr::str_extract(annot[1], '(?<=_dspl)[0-9]') %>% paste0('dtc deepSplit=', ., ')...\n')
    } else{
      cut_print = stringr::str_extract(annot[1], '(?<=_k)[0-9]') %>% paste0('k=', ., ')...\n')
    }
    test_print = paste0('Testing STclust assignment (w=', spw_print, ', ', cut_print)
    rm(spw_print, cut_print) # Clean environment
  } else{
    test_print = paste0('Testing metadata: ', annot[1], '...\n')
  }

  if(verbose != 0L){
    cat(crayon::green(test_print))
  }
  rm(test_print) # Clean environment

  # Extract sample annotations for each sample
  metas = tibble::tibble()
  for(sample_name in samples){
    metas = dplyr::bind_rows(metas,
                             tibble::tibble(samplename=sample_name,
                                            meta_orig=as.character(unique(x@spatial_meta[[sample_name]][[annot]]))))
  }

  # Find variable genes for each sample
  genes = tibble::tibble()
  for(sample_name in samples){
    # Get top variable genes
    genes_tmp = x@gene_meta[[sample_name]] %>%
      dplyr::arrange(desc(vst.variance.standardized))
    if(!is.null(topgenes)){
      genes_tmp = genes_tmp %>% dplyr::slice_head(n=topgenes)
    }
    genes_tmp = genes_tmp %>%
      dplyr::select(gene) %>% unlist() %>% as.vector()
    genes = dplyr::bind_rows(genes, tibble::tibble(samplename=sample_name, gene=genes_tmp))

    rm(genes_tmp) # Clean env
  }

  # Merge sample names, genes, and annotations in the same data frame
  gene_and_meta = dplyr::full_join(genes, metas, by='samplename', multiple='all')

  rm(genes, metas) # Clean env

  # Create "dictionary" with coded annotations (to avoid potentially problematic characters)
  meta_dict = tibble::tibble(orig_annot=unique(gene_and_meta[['meta_orig']]),
                             coded_annot=paste0('c', 1:length(unique(gene_and_meta[['meta_orig']]))))

  # Recode annotations gene_and_meta
  for(spotrow in 1:nrow(gene_and_meta)){
    gene_and_meta[['meta']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == gene_and_meta[['meta_orig']][spotrow] ]
  }
  rm(spotrow) # Clean env

  # Create a table with unique combinations of genes and annotations to test using parallelization
  combo_df = prepare_stdiff_combo(to_expand=gene_and_meta, pairwise=pairwise)
  rm(gene_and_meta) # Clean env

  ######## ######## ######## BEGIN NON-SPATIAL TESTS ######## ######## ########
  start_t = Sys.time()
  if(test_type == 'mm'){
    test_print = 'non-spatial mixed models'
  } else if(test_type == 't_test'){
    test_print = 'T-tests'
  } else{
    test_print = "Wilcoxon's tests"
  }
  if(verbose != 0L){
    cat(crayon::green(paste0('\tRunning ', test_print, '...\n')))
  }

  # Get annotations to run in parellel
  clusters = unique(combo_df[['meta1']])
  if(pairwise){
    clusters = unique(append(clusters, unique(combo_df[['meta2']])))
  }

  #Paralellize spaMM models
  if(is.null(cores)){
    cores = count_cores(nrow(combo_df))
  }
  non_sp_models = pbmcapply::pbmclapply(1:length(clusters), function(i){
    # Subset combo_df to those of a given cluster
    combo_clust = combo_df %>%
      dplyr::filter(meta1 == clusters[i])

    # Loop through samples
    res = list()
    for(sample_name in unique(combo_clust[['samplename']])){
      # Subset combinations per sample
      combo_tmp = combo_clust %>% dplyr::filter(samplename == sample_name)

      # Create data frame with expression, coordinate, and cluster data
      # Add group dummy column and select relevant columns
      expr_df = expandSparse(x@tr_counts[[sample_name]]) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::filter(gene %in% unique(combo_tmp[['gene']])) %>%
        tibble::column_to_rownames(var='gene') %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='libname') %>%
        dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                            tibble::add_column(group=1, .after='libname') %>%
                            dplyr::select(c('libname', 'group', 'ypos', 'xpos'), meta_orig:=!!annot[1]),. , by='libname') %>%
        tibble::column_to_rownames(var='libname')

      # Recode annotations gene_and_meta
      for(spotrow in 1:nrow(expr_df)){
        expr_df[['meta']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == expr_df[['meta_orig']][spotrow] ]
      }
      rm(spotrow) # Clean env

      # Remove original annotation and keep 'meta' as the coded annotations column
      expr_df = expr_df %>% dplyr::select(-c('meta_orig')) %>% dplyr::relocate(meta, .before=1)

      # Run non-spatial models
      if(verbose == 1L){
        system(sprintf('echo "%s"', crayon::green(paste0("\n\t\t",
                                                         sample_name, ", ",
                                                         meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == clusters[i] ],
                                                         " (", nrow(combo_clust), " tests)..."))))
      }
      if(test_type == 'mm'){
        res_tmp = non_spatial_de(expr_data=expr_df, combo=combo_tmp, pairwise=pairwise)
      } else if(test_type == 't_test' | test_type == 'wilcoxon'){
        res_tmp = stdiff_mean_test(expr_data=expr_df, combo=combo_tmp, pairwise=pairwise, test_type=test_type)
      }
      res = append(res, res_tmp)

      rm(res_tmp, expr_df, combo_tmp) # Clean env
    }

    return(res)
  }, mc.cores=cores)
  names(non_sp_models) = clusters

  # Summarize DE analyses
  result_de = tibble::tibble()
  for(i in names(non_sp_models)){
    for(mod in names(non_sp_models[[i]])){
      sample_tmp = non_sp_models[[i]][[mod]][['samplename']]
      gene_tmp = non_sp_models[[i]][[mod]][['gene']]
      avglfc_tmp = non_sp_models[[i]][[mod]][['avglfc']]
      # Compile results in a row
      df_tmp = tibble::tibble(sample=sample_tmp,
                              gene=gene_tmp,
                              avg_log2fc=avglfc_tmp,
                              cluster_1=meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == i ])
      if(pairwise){
        if(test_type == 'mm'){
          meta2_tmp = unique(non_sp_models[[i]][[mod]][['nonspmod']][['data']][['meta']]) %>% grep(i, ., value=T, invert=T)
        } else if(test_type == 't_test' | test_type == 'wilcoxon'){
          meta2_tmp = non_sp_models[[i]][[mod]][['meta2']]
        }
        df_tmp[['cluster_2']] = meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta2_tmp ]
        rm(meta2_tmp) # Clean env
      }

      if(test_type == 'mm'){
        df_tmp[['mm_p_val']] = spaMM::summary.HLfit(non_sp_models[[i]][[mod]][['nonspmod']],
                                                 details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]
      } else if(test_type == 't_test'){
        df_tmp[['ttest_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      } else if(test_type == 'wilcoxon'){
        df_tmp[['wilcox_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      }

      # Add results to list of data frames
      result_de = dplyr::bind_rows(result_de, df_tmp)

      rm(sample_tmp, gene_tmp, avglfc_tmp, df_tmp) # Clean env
    }
  }
  # Separate data frame into samples
  non_sp_de_tmp = list()
  for(i in unique(result_de[['sample']])){
    non_sp_de_tmp[[i]] = result_de %>% dplyr::filter(sample == i)
    # Adjust p-values
    if(test_type == 'mm'){
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['mm_p_val']], method=pval_adj)
    } else if(test_type == 't_test'){
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['ttest_p_val']], method=pval_adj)
    }else if(test_type == 'wilcoxon'){
      non_sp_de_tmp[[i]][['adj_p_val']] = p.adjust(non_sp_de_tmp[[i]][['wilcox_p_val']], method=pval_adj)
    }
    # Arrange results
    if(pairwise){
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, cluster_2, adj_p_val, desc(avg_log2fc))
    } else{
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, adj_p_val, desc(avg_log2fc))
    }
  }
  result_de = non_sp_de_tmp
  rm(non_sp_de_tmp) # Clean env

  end_t = difftime(Sys.time(), start_t, units='min')
  cat(crayon::green(paste0('\tCompleted ', test_print, ' (', round(end_t, 2), ' min).\n')))

  ######## ######## ######## BEGIN SPATIAL TESTS ######## ######## ########
  if(test_type == 'mm' & sp_topgenes > 0){
    # Subset list of genes based on adjusted p-values
    models_keep = non_sp_de[[sample_name]] %>%
      dplyr::filter(adj_p_val < pval_thr) %>%
      dplyr::arrange(adj_p_val) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_head(prop=sp_topgenes) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, cluster) %>%
      tibble::add_column(genexcluster=paste0(.$gene, '_m', .$cluster)) %>%
      dplyr::select(genexcluster) %>%
      unlist() %>%
      as.vector()
    non_sp_models = non_sp_models[models_keep]
    rm(models_keep) # Clean environment


    start_t = Sys.time()
    cat(crayon::yellow(paste0('\t\tRunning ', length(non_sp_models), ' spatial tests...')))

    # Run models
    sp_models = spatial_de(non_sp_mods=non_sp_models, sp_topgenes=sp_topgenes, verbose=verbose, cores=cores)

    rm(expr_df) # Clean environment
    end_t = difftime(Sys.time(), start_t, units='min')
    cat(crayon::green(paste0('Completed (', round(end_t, 2), ' min).\n')))

    # Compile results of spatial models and merge with non-spatial results
    sp_de = tibble::tibble()
    for(i in names(sp_models)){
      # Check that model results exists and extract values, gene names, and clusters tested
      if(class(sp_models[[i]][['sph']]) == 'lme'){
        sph_summ = summary(sp_models[[i]][['sph']])
        meta_sph_tmp = names(sp_models[[i]][['sph']]$coefficients$fixed) %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_sph_tmp = gsub(paste0('_m', meta_sph_tmp, '$'), '', i)
      } else{
        sph_summ = list(AIC=-9999, BIC=-9999, tTable=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
        meta_sph_tmp = sp_models[[i]][['sph']][[2]] %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_sph_tmp = gsub(paste0('_m', meta_sph_tmp, '$'), '', i)
      }
      if(class(sp_models[[i]][['exp']]) == 'lme'){
        exp_summ = summary(sp_models[[i]][['exp']])
        meta_exp_tmp = names(sp_models[[i]][['exp']]$coefficients$fixed) %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_exp_tmp = gsub(paste0('_m', meta_exp_tmp, '$'), '', i)
      } else{
        exp_summ = list(AIC=-9999, BIC=-9999, tTable=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
        meta_exp_tmp = sp_models[[i]][['exp']][[2]] %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_exp_tmp = gsub(paste0('_m', meta_exp_tmp, '$'), '', i)
      }

      # Make sure genes and clusters from both exponential and spherical models match
      if((meta_exp_tmp == meta_sph_tmp) & (gene_exp_tmp == gene_sph_tmp)){
        sp_de = dplyr::bind_rows(sp_de,
                                 tibble::tibble(sample=sample_name,
                                                gene=gene_exp_tmp,
                                                cluster=meta_exp_tmp,
                                                sph_aic=sph_summ[['AIC']],
                                                sph_bic=sph_summ[['BIC']],
                                                sph_p_val=sph_summ[['tTable']][2, 5],
                                                exp_aic=exp_summ[['AIC']],
                                                exp_bic=exp_summ[['BIC']],
                                                exp_p_val=exp_summ[['tTable']][2, 5]
                                 ))
      } else{
        stop('Gene and cluster do not match between exponential and spherical models...')
      }
    }

    # Adjust spatial p-values
    sph_adj_df = sp_de %>%
      dplyr::filter(!is.na(sph_p_val) & sph_p_val != -9999) %>%
      dplyr::select(sample, gene, cluster, sph_p_val) %>%
      tibble::add_column(sph_adj_p_val=p.adjust(.[['sph_p_val']], method=pval_adj)) %>%
      dplyr::select(-sph_p_val)

    exp_adj_df = sp_de %>%
      dplyr::filter(!is.na(exp_p_val) & exp_p_val != -9999) %>%
      dplyr::select(sample, gene, cluster, exp_p_val) %>%
      tibble::add_column(exp_adj_p_val=p.adjust(.[['exp_p_val']], method=pval_adj)) %>%
      dplyr::select(-exp_p_val)

    # Merge spatial and non-spatial DE results
    non_sp_de[[sample_name]] = dplyr::full_join(non_sp_de[[sample_name]], sp_de, by=c('sample', 'gene', 'cluster')) %>%
      dplyr::full_join(., sph_adj_df, by=c('sample', 'gene', 'cluster')) %>%
      dplyr::full_join(., exp_adj_df, by=c('sample', 'gene', 'cluster')) %>%
      dplyr::relocate(sph_adj_p_val, .after=sph_p_val) %>%
      dplyr::relocate(exp_adj_p_val, .after=exp_p_val) %>%
      dplyr::mutate(sph_adj_p_val=dplyr::case_when(sph_p_val == -9999 ~ -9999, TRUE ~ sph_adj_p_val)) %>%
      dplyr::mutate(exp_adj_p_val=dplyr::case_when(exp_p_val == -9999 ~ -9999, TRUE ~ exp_adj_p_val)) %>%
      dplyr::mutate(spatialmod=dplyr::case_when(is.na(sph_p_val) ~ '2', TRUE ~ '1')) %>% # To put on top of table the genes with spatial models
      dplyr::arrange(spatialmod, cluster, adj_p_val) %>%
      dplyr::select(-spatialmod)

    rm(sp_de) # Clean environment
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  cat(crayon::green(paste0('STdiff completed in ', round(end_t, 2), ' min.\n')))

  return(result_de)
}


# Helpers ----------------------------------------------------------------------

##
# non_spatial_de
# Run non-spatial linear models for a set of gene x cluster combinations
# @param expr_df a data frame (rows=spots/cells) with normalized expression data
# as well as coordinates (x,y), and a group variable for random effects.
# @param combo a data frame with two columns (meta and gene) containing all
# combinations of gene x cluster to test.
# @param cores number of cores to use in parallelization. If `NULL`, the number of
# cores is to use is detected automatically.
# @return a list containing spaMM models
#
non_spatial_de = function(expr_data=NULL, combo=NULL, pairwise=NULL){
  models_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Create non-spatial model
    res_mod = spaMM::fitme(formula=as.formula('exprval~meta'), data=expr_tmp, method="REML")

    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['nonspmod']] = res_mod
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp

  }
  return(models_ls)
}


##
# spatial_de
# Run spatial linear models with exponential and spherical covariance structures
# @param non_sp_mods a list of non-spatial models to update with covariance structures
# @param sp_topgenes an integer indicating the number of top DE genes to apply spatial models
# @param cores number of cores to use in parallelization. If `NULL`, the number of
# cores is to use is detected automatically.
# @return a list of spatial models
#
spatial_de = function(non_sp_mods=NULL, sp_topgenes=NULL, verbose=verbose, cores=NULL){
  # Paralellize lme models
  if(is.null(cores)){
    cores = count_cores(length(non_sp_mods))
  }
  sp_mods = parallel::mclapply(names(non_sp_mods), function(i){
    res_ls = list() # Initiates list to storage results
    # Extract cluster identifiers
    meta_tmp = rownames(summary(non_sp_mods[[i]])[['tTable']])[2] %>% gsub('meta', '', .)
    # Extract data from the same model object
    expr_tmp = non_sp_mods[[i]][['data']]
    # Extract gene name from data stored in lme object (should have workaround to avoid dashes in gene names as nlme does not handle them well)
    genenodash_tmp = names(expr_tmp)[5]

    # Rewrite the lme command call so that it works on the new environment
    non_sp_mods[[i]][['call']] = str2lang('nlme::lme.formula(fixed=as.formula(paste0("`", genenodash_tmp, "`", "~meta")), data=expr_tmp, random=~1|group, method="REML")')

    if(verbose == 2){
      system(sprintf('echo "%s"', crayon::green(i)))
    }

    # Run spherical model and catch error if no convergence
    sph_out = tryCatch({
      update(non_sp_mods[[i]], correlation=nlme::corSpher(form=~xpos+ypos, nugget=T, metric='euclidean'), method='REML', control=nlme::lmeControl(returnObject=F, msVerbose=F, gradHess=T))
    }, error=function(err){return(err)})

    if(any(class(sph_out) == 'simpleError')){
      res_ls[['sph']][[1]] = 'no_conv'
      res_ls[['sph']][[2]] = names(non_sp_mods[[i]]$coefficients$fixed)
    } else{
      res_ls[['sph']] = sph_out
    }

    # Run exponential model and catch error if no convergence
    exp_out = tryCatch({
      update(non_sp_mods[[i]], correlation=nlme::corExp(form=~xpos+ypos, nugget=T, metric='euclidean'), method='REML', control=nlme::lmeControl(returnObject=F, msVerbose=F, gradHess=T))
    }, error=function(err){return(err)})

    if(any(class(exp_out) == 'simpleError')){
      res_ls[['exp']][[1]] = 'no_conv'
      res_ls[['exp']][[2]] = names(non_sp_mods[[i]]$coefficients$fixed)
    } else{
      res_ls[['exp']] = exp_out
    }
    return(res_ls)
  }, mc.cores=cores, mc.preschedule=F)

  if(length(non_sp_mods) == length(sp_mods)){
    names(sp_mods) = names(non_sp_mods) # Transfer names of non-spatial list to  spatial list
  } else{
    stop('Looks like the number of non-spatial models do not match the number of spatial models')
  }

  return(sp_mods)
}


##
# create_stdiff_combo
# @return data frame with combinations
#
prepare_stdiff_combo = function(to_expand=NULL, pairwise=NULL){
  combo_df = tibble::tibble()
  for(sample_name in unique(to_expand[['samplename']])){
    # Extract annotations for a given sample
    annots_tmp = to_expand %>% dplyr::filter(samplename == sample_name) %>%
      dplyr::select('meta') %>% unlist() %>% unique()
    if(pairwise){
      # Get unique combinations of each two annotation, without repeating (e.g., c1 vs c2 and c2 vs c1)
      # Get combinations for each cluster
      combo_meta = tibble::tibble()
      for(cl_tmp in annots_tmp){
        combo_meta = dplyr::bind_rows(combo_meta,
                                      expand.grid(cl_tmp,
                                                  grep(cl_tmp, annots_tmp, value=T, invert=T),
                                                  stringsAsFactors=F))
        annots_tmp = grep(cl_tmp, annots_tmp, value=T, invert=T)
      }
      names(combo_meta) = c('meta1', 'meta2')
      rm(cl_tmp) #Clean env

      # Add genes to each combination
      combo_meta_gene = tibble::tibble()
      for(comb in 1:nrow(combo_meta)){
        combo_meta_gene = dplyr::bind_rows(combo_meta_gene,
                                           tibble::tibble(combo_meta[comb, 1],
                                                          combo_meta[comb, 2],
                                                          gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                                            dplyr::select(gene) %>% unlist() %>% unique()))
      }
      combo_meta = combo_meta_gene
      rm(combo_meta_gene, comb) # Clean env
    } else{
      combo_meta = expand.grid(meta1=annots_tmp, meta2='other',
                               gene=to_expand %>% dplyr::filter(samplename == sample_name) %>%
                                 dplyr::select(gene) %>% unlist() %>% unique(),
                               stringsAsFactors=F)
    }
    rm(annots_tmp) # Clean env
    combo_df = dplyr::bind_rows(combo_df,
                                combo_meta %>%
                                  tibble::add_column(samplename=sample_name, .before=1))
  }
  combo_df = as.data.frame(combo_df) %>%
    dplyr::arrange(samplename, meta1, meta2, gene)

  return(combo_df)
}


##
# stdiff_mean_test
# Performs a simple t-test or Wilcoxon between spots/cells assigned to clusters
#
stdiff_mean_test = function(expr_data=NULL, combo=NULL, pairwise=NULL, test_type=NULL){
  test_ls = list()
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]

    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      expr_tmp = expr_data %>%
        dplyr::filter(meta %in% c(meta1_tmp, meta2_tmp))
    } else{
      expr_tmp = expr_data %>%
        dplyr::mutate(meta=factor(dplyr::case_when(meta != meta1_tmp ~ 'other', TRUE ~ meta),
                                  levels=c('other', meta1_tmp)))
    }

    # Get relevant gene column
    expr_tmp = expr_tmp %>%
      dplyr::select(c('group', 'ypos', 'xpos', 'meta'), exprval:=!!gene_tmp)

    # Calculate log-fold change
    avgexpr_cl1 = mean(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select(exprval) %>% unlist())
    avgexpr_cl2 = mean(expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select(exprval) %>% unlist())
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)

    # Perform test for diferences in mean
    if(test_type == 't_test'){
      res_test = t.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                        expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())
    } else if(test_type == 'wilcoxon'){
      # Create non-spatial model
      res_test = wilcox.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                             expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())
    }

    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['mean_test']] = res_test
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    test_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
  }
  return(test_ls)
}

