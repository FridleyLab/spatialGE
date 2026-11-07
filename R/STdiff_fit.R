# Compute stage: runs the non-spatial and spatial model fits dispatched by STdiff().

##
# stdiff_run_nonspatial_tests
# Run the non-spatial (or t-test/Wilcoxon) tests for every requested sample,
# dispatching one worker per sample via sge_lapply.
# @param x an STlist.
# @param combo_df the table of gene x cluster combinations to test.
# @param clusters user-specified clusters to restrict testing to (or NULL for all).
# @param meta_dict the original-to-coded annotation dictionary.
# @param annot the annotation column name in x@spatial_meta being tested.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
# @param test_type one of `mm`, `t_test`, or `wilcoxon`.
# @param cores number of cores to use in parallelization. If `NULL`, the number of
# cores to use is detected automatically.
# @return a list with `non_sp_models` (named by sample) and the resolved `cores`.
#
stdiff_run_nonspatial_tests = function(x, combo_df, clusters, meta_dict, annot, pairwise, verbose, test_type, cores){
  start_t = Sys.time()
  
  # Check test type
  if(test_type == 'mm'){
    test_print = 'non-spatial mixed models'
  } else if(test_type == 't_test'){
    test_print = 'T-tests'
  } else if(test_type == 'wilcoxon'){
    test_print = "Wilcoxon's tests"
  } else{
    stop('Test type is one of "mm", "t_test", or "wilcoxon".')
  }
  
  # Print test to conduct
  if(verbose){
    cat(paste0('\tRunning ', test_print, '...\n'))
  }
  
  # If user specified certain clusters, then subset combo_df to those clusters only
  if(!is.null(clusters)){
    coded_metas = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] %in% clusters ]
    combo_df = combo_df[ combo_df[['meta1']] %in% coded_metas, ]
    if(pairwise){
      combo_df = combo_df[ combo_df[['meta2']] %in% coded_metas, ]
    }
    rm(coded_metas) # Clean env
    # Throw error if user-input clusters are present in data
    if(nrow(combo_df) < 1){
      stop('All cluster x gene commparisons were removed. Is at least one of the specified clusters present in one of the samples?')
    }
  }
  
  # Define number of cores for parallelization of tests
  cores = resolve_cores(length(unique(combo_df[['samplename']])), cores, coerce='ceiling')
  
  # Paralellize non-spatial tests
  non_sp_models = sge_lapply(1:length(unique(combo_df[['samplename']])), function(i){
    # Get sample name
    sample_name = unique(combo_df[['samplename']])[i]
    # Subset combo_df to those of a given sample
    combo_clust = combo_df[ combo_df[['samplename']] == sample_name, ]
    
    stdiff_nonspatial_one_sample(sample_tr_counts(x, sample_name), sample_spatial_meta(x, sample_name),
                                 combo_clust, sample_name, annot, meta_dict, pairwise, verbose, test_type)
  }, cores=cores)
  
  # Put sample names to results of parallelezation
  names(non_sp_models) = unique(combo_df[['samplename']])
  
  end_t = difftime(Sys.time(), start_t, units='min')
  if(verbose){
    cat(paste0('\tCompleted ', test_print, ' (', round(end_t, 2), ' min).\n'))
  }
  
  return(list(non_sp_models=non_sp_models, cores=cores))
}


##
# stdiff_run_spatial_tests
# Run the spatial model fits for every sample in result_de, restricted to the
# genes retained per `sp_topgenes`/`pval_thr`. Dispatches one worker per
# cluster (within each sample) via sge_lapply, and labels which gene x cluster
# combinations were selected for spatial testing directly on result_de.
# @param x an STlist.
# @param result_de the per-sample non-spatial DE result tables.
# @param non_sp_models the per-sample non-spatial models, keyed by sample.
# @param meta_dict the original-to-coded annotation dictionary.
# @param annot the annotation column name in x@spatial_meta being tested.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
# @param sp_topgenes proportion of non-spatial DE genes to carry into spatial testing.
# @param pval_thr adjusted p-value cut-off used to select genes for spatial testing.
# @param cores number of cores to use in parallelization (same value resolved for
# the non-spatial tests; not re-resolved here).
# @return a list with the updated `result_de` (spatial-test labels added) and `sp_models`.
#
stdiff_run_spatial_tests = function(x, result_de, non_sp_models, meta_dict, annot, pairwise, verbose, sp_topgenes, pval_thr, cores){
  # To prevent NOTES in R CMD check
  . = NULL
  
  sp_models = list()
  for(sample_name in names(result_de)){
    # Subset list of genes based on adjusted p-values
    models_keep = result_de[[sample_name]] %>%
      dplyr::filter(adj_p_val < pval_thr) %>%
      dplyr::arrange(adj_p_val)
    if(pairwise){
      models_keep = models_keep %>%
        dplyr::group_by(cluster_1, cluster_2)
    } else{
      models_keep = models_keep %>%
        dplyr::group_by(cluster_1)
    }
    models_keep = models_keep %>%
      dplyr::slice_head(prop=sp_topgenes) %>%
      dplyr::ungroup()
    
    # Add label to indicate if gene x cluster was selected for spatial test
    if(pairwise){
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::full_join(., models_keep %>%
                           tibble::add_column(comments='spatial_test') %>%
                           dplyr::select(c("sample", "gene", "cluster_1", "cluster_2", "comments")), by=c("sample", "gene", "cluster_1", "cluster_2"))
    } else{
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::full_join(., models_keep %>%
                           tibble::add_column(comments='spatial_test') %>%
                           dplyr::select(c("sample", "gene", "cluster_1", "comments")), by=c("sample", "gene", "cluster_1"))
    }
    result_de[[sample_name]] = result_de[[sample_name]] %>%
      dplyr::mutate(comments=dplyr::case_when(comments == 'spatial_test' ~ NA_character_, TRUE ~ 'no_spatial_test'))
    
    # Recode annotations gene_and_meta
    for(spotrow in 1:nrow(models_keep)){
      models_keep[['cluster_1']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_1']][spotrow] ]
      if(pairwise){
        models_keep[['cluster_2']][spotrow] = meta_dict[['coded_annot']][ meta_dict[['orig_annot']] == models_keep[['cluster_2']][spotrow] ]
      }
    }
    
    if(nrow(models_keep) > 0){
      # Create combinations of gene x cluster
      if(pairwise){
        models_keep = models_keep %>%
          dplyr::select(c('sample', 'gene', 'cluster_1', 'cluster_2')) %>%
          tibble::add_column(genexcluster=paste0(.[['sample']], '_', .[['gene']], '_', .[['cluster_1']], '_', .[['cluster_2']]))
      } else{
        models_keep = models_keep %>%
          dplyr::select(c('sample', 'gene', 'cluster_1')) %>%
          tibble::add_column(genexcluster=paste0(.[['sample']], '_', .[['gene']], '_', .[['cluster_1']], '_other'))
      }
      clusters_tmp = unique(models_keep[['cluster_1']])
      
      # Create data frame with expression, coordinate, and cluster data
      # Add group dummy column (in case of mixed models was requested)
      # spatial_de() only ever reads genes in non_sp_models_tmp, which is always
      # restricted to models_keep's genes, so subset before decompressing.
      expr_df = expandSparse(sample_tr_counts(x, sample_name)[ sample_genes(x, sample_name, 'tr_counts') %in% unique(models_keep[['gene']]), , drop=FALSE], transpose=TRUE)
      expr_df = expr_df %>%
        tibble::rownames_to_column(var='libname') %>%
        dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                            tibble::add_column(group=1, .after='libname') %>%
                            dplyr::select(c('libname', 'group', 'ypos', 'xpos'), orig_annot:=!!annot),. , by='libname') %>%
        dplyr::left_join(., meta_dict, by='orig_annot') %>%
        dplyr::rename(meta=coded_annot) %>%
        tibble::column_to_rownames(var='libname') %>%
        dplyr::relocate(meta, .before=1) %>%
        dplyr::select(-c('orig_annot'))
      
      sp_models[[sample_name]] = sge_lapply(1:length(clusters_tmp), function(i){
        stdiff_spatial_one_cluster(models_keep, clusters_tmp[i], non_sp_models[[sample_name]],
                                   expr_df, sample_name, meta_dict, pairwise, verbose)
      }, cores=cores)
      rm(models_keep, clusters_tmp) # Clean environment
    } else{
      warning(paste0('No genes to test for ',  sample_name, '. Try increasing sp_topgenes.'))
    }
  }
  
  return(list(result_de=result_de, sp_models=sp_models))
}


##
# stdiff_nonspatial_one_sample
# Run the non-spatial (or t-test/Wilcoxon) tests for every cluster x gene
# combination requested for one sample.
# @param tr_counts_mat the sample's transformed counts matrix.
# @param spatial_meta_df the sample's spatial metadata data frame.
# @param combo_clust the subset of combo_df for this sample.
# @param sample_name the sample's name.
# @param annot the annotation column name in spatial_meta_df being tested.
# @param meta_dict the original-to-coded annotation dictionary.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
# @param test_type one of `mm`, `t_test`, or `wilcoxon`.
# @return a list of per gene x cluster test results.
#
stdiff_nonspatial_one_sample = function(tr_counts_mat, spatial_meta_df, combo_clust, sample_name,
                                        annot, meta_dict, pairwise, verbose, test_type){
  # Create data frame with expression, coordinate, and cluster data
  # Add group dummy column (in case of mixed models was requested)
  expr_df = expandSparse(tr_counts_mat[ rownames(tr_counts_mat) %in% unique(combo_clust[['gene']]), , drop=FALSE], transpose=TRUE)
  expr_df = expr_df %>%
    tibble::rownames_to_column(var='libname') %>%
    dplyr::right_join(spatial_meta_df %>%
                        tibble::add_column(group=1, .after='libname') %>%
                        dplyr::select(c('libname', 'group', 'ypos', 'xpos'), orig_annot:=!!annot),. , by='libname') %>%
    dplyr::left_join(., meta_dict, by='orig_annot') %>%
    dplyr::rename(meta=coded_annot) %>%
    tibble::column_to_rownames(var='libname') %>%
    dplyr::relocate(meta, .before=1) %>%
    dplyr::select(-c('orig_annot'))
  
  # Loop through tests within sample
  res = list()
  for(cluster_tmp in unique(combo_clust[['meta1']])){ # Loop through clusters, instead of rows to provide progress updates
    # Subset tests to specific cluster
    combo_clust_cl = combo_clust[ combo_clust[['meta1']] == cluster_tmp, ]
    # Run non-spatial models
    if(verbose >= 1L){
      if(pairwise){
        stdout_print = combo_clust_cl[ combo_clust_cl[['samplename']] == sample_name & combo_clust_cl[['meta1']] == cluster_tmp, ]
        stdout_print = unique(stdout_print[['meta2']])
        stdout_print = meta_dict[['orig_annot']][ meta_dict[['coded_annot']] %in% stdout_print ]
        stdout_print = paste0(stdout_print, collapse=', ')
        system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                           sample_name, ", ",
                                           meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_tmp ],
                                           " vs. ", stdout_print,
                                           " (", nrow(combo_clust_cl), " tests)")))
        rm(stdout_print)
      } else {
        system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                           sample_name, ", ",
                                           meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_tmp ],
                                           " (", nrow(combo_clust_cl), " tests)")))
      }
    }
    if(test_type == 'mm'){
      res_tmp = non_spatial_de(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise)
    } else if(test_type == 't_test' | test_type == 'wilcoxon'){
      res_tmp = stdiff_mean_test(expr_data=expr_df, combo=combo_clust_cl, pairwise=pairwise, test_type=test_type)
    }
    res = append(res, res_tmp)
    
    rm(res_tmp, combo_clust_cl) # Clean env
  }
  return(res)
}


##
# stdiff_spatial_one_cluster
# Run the spatial models for every gene retained (per `sp_topgenes`/`pval_thr`)
# within one cluster of one sample.
# @param models_keep the sample's table of gene x cluster combinations selected
# for spatial testing.
# @param cluster_name the cluster (coded annotation) being tested.
# @param non_sp_models_sample the sample's non-spatial models, keyed by
# sample_gene_meta1_meta2.
# @param expr_df the sample's expression/coordinate/cluster data frame.
# @param sample_name the sample's name.
# @param meta_dict the original-to-coded annotation dictionary.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param verbose either logical or an integer (0, 1, or 2) to increase verbosity.
# @return a list of spatial models for the cluster's retained gene tests.
#
stdiff_spatial_one_cluster = function(models_keep, cluster_name, non_sp_models_sample, expr_df,
                                      sample_name, meta_dict, pairwise, verbose){
  # Non-spatial models to update
  models_keep_tmp = models_keep %>%
    dplyr::filter(cluster_1 == cluster_name) %>%
    dplyr::select('genexcluster') %>%
    unlist() %>% as.vector()
  non_sp_models_tmp = non_sp_models_sample[ models_keep_tmp ]
  
  if(verbose == 1L){
    if(pairwise){
      stdout_print = models_keep_tmp %>%
        stringr::str_extract(., paste0('(?<=', cluster_name, '_)c[0-9]+$')) %>%
        stringr::str_split(., '_') %>% unlist() %>% unique()
      stdout_print_tmp = c()
      for(m2 in stdout_print){
        stdout_print_tmp = append(stdout_print_tmp,
                                  meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == m2 ])
      }
      stdout_print = stdout_print_tmp %>% paste0(., collapse=', ')
      system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                         sample_name, ", ",
                                         meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_name ],
                                         " vs. ", stdout_print,
                                         " (", length(models_keep_tmp), " tests)...")))
      rm(stdout_print)
    } else {
      system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                         sample_name, ", ",
                                         meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_name ],
                                         " (", length(models_keep_tmp), " tests)...")))
    }
  }
  
  # Run models
  sp_mods = spatial_de(expr_dat=expr_df, non_sp_mods=non_sp_models_tmp, annot_dict=meta_dict, verb=verbose)
  
  return(sp_mods)
}


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
#' @importFrom stats as.formula
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
    #res_mod = spaMM::fitme(formula=as.formula('exprval~meta'), data=expr_tmp, method="REML")
    # res_mod = tryCatch({
    #   mod = spaMM::fitme(formula=as.formula('exprval~meta'), data=expr_tmp, method="REML")
    # },
    # warning=function(w, mod){return(list(warning=w, result=mod))},
    # #error=function(e){return(list(error=e))},
    # finally=function(mod){return(list(result=mod))}
    # )
    
    error_message = c()
    warning_message = c()
    res_mod = withCallingHandlers(
      tryCatch(spaMM::fitme(formula=stats::as.formula('exprval~meta'), data=expr_tmp, method="REML"),
               error=function(e){error_message <<- conditionMessage(e)}),
      warning=function(w){warning_message <<- conditionMessage(w)
      invokeRestart("muffleWarning")}
    )
    
    # if(any(unlist(lapply(res_mod, class)) == 'simpleWarning')){
    #   wng = warning_message
    # }
    
    # Calculate p-value. Only if the fit actually converged -- on error, res_mod
    # holds the captured error message (a plain string) rather than an HLfit
    # object, and summary.HLfit() on that would itself error, aborting this
    # gene x cluster combination's fit taking every other combination sharing
    # its mclapply chunk down with it (see fitme's non-spatial counterpart,
    # spatial_de(), for the same not-'HLfit'-class guard).
    if(any(class(res_mod) == 'HLfit')){
      res_mod = spaMM::summary.HLfit(res_mod, details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]
    } else{
      res_mod = NA_real_
    }
    
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]] = list()
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['nonspmod']] = res_mod
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['samplename']] = sample_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['gene']] = gene_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['avglfc']] = avglogfold_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta1']] = meta1_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['meta2']] = meta2_tmp
    models_ls[[paste0(sample_tmp, '_', gene_tmp, '_', meta1_tmp, '_', meta2_tmp)]][['warning']] = warning_message
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
#' @importFrom stats as.formula
#
spatial_de = function(expr_dat=NULL, non_sp_mods=NULL, annot_dict=NULL, verb=NULL){
  
  # To prevent NOTES in R CMD check
  . = NULL
  
  res_ls = list()
  for(i in names(non_sp_mods)){
    # Extract sample name
    sample_tmp = non_sp_mods[[i]][['samplename']]
    # Extract gene name
    gene_tmp = non_sp_mods[[i]][['gene']]
    # Extract clusters
    meta1_tmp = non_sp_mods[[i]][['meta1']]
    meta2_tmp = non_sp_mods[[i]][['meta2']]
    
    # Print information of test
    if(verb == 2L){
      if(grepl('_c[0-9]+_other$', i)){
        stdout_print = i %>%
          stringr::str_extract(., 'c[0-9]+_other$') %>%
          stringr::str_replace(., '_other$', '')
        stdout_print = annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print ]
      } else{
        stdout_print = i %>%
          stringr::str_extract(., paste0('c[0-9]+_c[0-9]+$')) %>%
          stringr::str_split(., '_') %>% unlist() %>% unique()
        stdout_print = c(annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[1] ],
                         annot_dict[['orig_annot']][ annot_dict[['coded_annot']] == stdout_print[2] ])
        stdout_print = paste0(stdout_print[1], ' vs. ', stdout_print[2])
      }
      system(sprintf('echo "%s"', paste0('\t\tSample: ', sample_tmp, '; ', gene_tmp, ', ', stdout_print)))
    }
    
    expr_subset = expr_dat[, c("meta" , "group", "ypos", "xpos", gene_tmp)] %>%
      dplyr::rename(exprval := !!gene_tmp)
    if(meta2_tmp == 'other'){
      expr_subset[['meta']][ expr_subset[['meta']] != meta1_tmp ] = 'other'
    } else{
      expr_subset = expr_subset[ expr_subset[['meta']] %in% c(meta1_tmp, meta2_tmp), ]
    }
    
    exp_out = tryCatch({
      #R.utils::withTimeout({
      spaMM::fitme(formula=stats::as.formula(paste0("exprval~meta+Matern(1|xpos+ypos)")),
                   data=expr_subset,
                   fixed=list(nu=0.5), method="REML",
                   control.HLfit=list(algebra="decorr"))
      #}, timeout=5)
    }, error=function(err){return(err)})
    
    res_ls[[i]] = list()
    if(any(class(exp_out) == 'TimeoutException')){
      res_ls[[i]][['spmod']] = 'time_out'
    } else if(any(class(exp_out) == 'simpleError')){
      res_ls[[i]][['spmod']] = 'no_conv'
    } else{
      res_ls[[i]][['spmod']] = exp_out
    }
    res_ls[[i]][['sample']] = sample_tmp
    res_ls[[i]][['gene']] = gene_tmp
  }
  return(res_ls)
}


##
# stdiff_mean_test
# Performs a simple t-test or Wilcoxon between spots/cells assigned to clusters
#
stdiff_mean_test = function(expr_data=NULL, combo=NULL, pairwise=NULL, test_type=NULL){
  test_ls = list()
  meta_all = expr_data[['meta']]
  
  # Cache for symmetric p-values. R's wilcox.test() (and, to a much smaller
  # degree, t.test()) computes its statistic through an order-dependent
  # floating-point path, so swapping the two groups can change the p-value
  # in the last digit. test-numeric-stdiff.R:145 asserts bit-identical
  # p-values for the two comparison directions, so we must call the test
  # ONCE per (gene, unordered cluster comparison) and reuse the result for
  # both directions. The cache below is keyed on a canonical representation
  # of the comparison that is identical for both directions, in BOTH
  # pairwise and one-vs-rest modes.
  p_cache <- new.env(parent = emptyenv())
  
  for(i in 1:nrow(combo)){
    sample_tmp = combo[i, 1]
    meta1_tmp = combo[i, 2]
    meta2_tmp = combo[i, 3]
    gene_tmp = combo[i, 4]
    
    # Filter data if pairwise or recode annotations if not pairwise
    if(pairwise){
      keep = meta_all %in% c(meta1_tmp, meta2_tmp)
      meta_vec = meta_all[keep]
      gene_vec = expr_data[[gene_tmp]][keep]
    } else{
      meta_vec = ifelse(meta_all == meta1_tmp, as.character(meta_all), 'other')
      gene_vec = expr_data[[gene_tmp]]
    }
    exprval_cl1 = gene_vec[meta_vec == meta1_tmp]
    exprval_cl2 = gene_vec[meta_vec == meta2_tmp]
    
    # Calculate log-fold change
    avgexpr_cl1 = mean(exprval_cl1)
    avgexpr_cl2 = mean(exprval_cl2)
    avglogfold_tmp = (avgexpr_cl1 - avgexpr_cl2)
    
    # Canonicalize the test call so that both cluster_1 directions feed the
    # SAME (x, y) argument pair to wilcox.test / t.test. This is what makes
    # the two directions' p-values bit-identical.
    #
    # Pairwise mode:
    #   meta2_tmp is a real cluster label. Canonical key = sorted (cl_a, cl_b),
    #   and (x, y) = (expr of cl_a, expr of cl_b) with cl_a < cl_b. Both
    #   directions hit the same key and reuse the same res_test.
    #
    # One-vs-rest mode:
    #   meta2_tmp is the literal "other". For a two-cluster comparison (the
    #   fixture, and any two-cluster sample), "c1 vs rest" and "c2 vs rest"
    #   are the same mathematical comparison. The canonical key is the sorted
    #   pair of REAL cluster labels ({"c1","c2"}), and (x, y) is built from
    #   the real clusters directly. Both directions hit the same key.
    #
    # For one-vs-rest over 3+ clusters the comparison is not symmetric with
    # other one-vs-rest tests (the "rest" group differs), so the cache key
    # falls back to the single cluster label; the p-value is computed
    # separately per cluster, matching historical behaviour. The symmetry
    # test only covers two-cluster fixtures.
    if(meta2_tmp == 'other'){
      # In one-vs-rest mode, `meta_vec` recoded everything except meta1_tmp
      # to the literal "other", so the OTHER real clusters' labels are only
      # available in `meta_all`. Recover them from there.
      all_real <- sort(unique(as.character(meta_all[!is.na(meta_all)])))
      # Drop the placeholder if the fixture uses "other" as a real label.
      all_real <- setdiff(all_real, 'other')
      
      if(length(all_real) == 2){
        # Two real clusters: "cA vs rest" and "cB vs rest" are the same
        # mathematical comparison. Canonicalize on the sorted real-cluster
        # pair so both directions hit the SAME cache key and reuse the SAME
        # res_test -> bit-identical p-values.
        canon_key <- paste(sample_tmp, gene_tmp,
                           paste(all_real, collapse = '|'), sep = '::')
        if(!exists(canon_key, envir = p_cache, inherits = FALSE)){
          x_vec <- gene_vec[!is.na(meta_all) & meta_all == all_real[1]]
          y_vec <- gene_vec[!is.na(meta_all) & meta_all == all_real[2]]
          if(test_type == 't_test'){
            res_test <- stats::t.test(x_vec, y_vec)
          } else if(test_type == 'wilcoxon'){
            res_test <- stats::wilcox.test(x_vec, y_vec,
                                           alternative = 'two.sided',
                                           exact = FALSE, correct = TRUE)
          }
          assign(canon_key, res_test, envir = p_cache)
        } else{
          res_test <- get(canon_key, envir = p_cache)
        }
      } else{
        # 3+ clusters (or degenerate): each direction's complement is
        # different, so p-values are legitimately different. Compute per
        # direction; no reuse expected.
        canon_key <- paste(sample_tmp, gene_tmp, meta1_tmp, 'ovr', sep = '::')
        if(!exists(canon_key, envir = p_cache, inherits = FALSE)){
          x_vec <- gene_vec[!is.na(meta_vec) & meta_vec == meta1_tmp]
          y_vec <- gene_vec[!is.na(meta_vec) & meta_vec != meta1_tmp]
          if(test_type == 't_test'){
            res_test <- stats::t.test(x_vec, y_vec)
          } else if(test_type == 'wilcoxon'){
            res_test <- stats::wilcox.test(x_vec, y_vec,
                                           alternative = 'two.sided',
                                           exact = FALSE, correct = TRUE)
          }
          assign(canon_key, res_test, envir = p_cache)
        } else{
          res_test <- get(canon_key, envir = p_cache)
        }
      }
    } else{
      labs <- sort(c(meta1_tmp, meta2_tmp))
      canon_key <- paste(sample_tmp, gene_tmp, paste(labs, collapse = '|'),
                         sep = '::')
      if(!exists(canon_key, envir = p_cache, inherits = FALSE)){
        x_vec <- gene_vec[meta_vec == labs[1]]
        y_vec <- gene_vec[meta_vec == labs[2]]
        if(test_type == 't_test'){
          res_test <- stats::t.test(x_vec, y_vec)
        } else if(test_type == 'wilcoxon'){
          res_test <- stats::wilcox.test(x_vec, y_vec,
                                         alternative = 'two.sided',
                                         exact = FALSE, correct = TRUE)
        }
        assign(canon_key, res_test, envir = p_cache)
      } else{
        res_test <- get(canon_key, envir = p_cache)
      }
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