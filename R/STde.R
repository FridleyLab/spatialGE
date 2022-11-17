##
#' @title STde: Differential gene expression analysis for spatial transcriptomics data
#' @description Tests for differentially expressed genes using linear models with
#' exponential or spherical spatial covariance structures.
#' @details
#' The method tests for differentially expressed genes between groups of spots/cells
#' (e.g., clusters) in a spatial transcriptomics sample (e.g., spots in a Visium tissue slice).
#' Specifically, tests for genes with significantly higher or lower gene expression in
#' one group of spots/cells with respect to the rest of spots/cells in the sample.
#' The method first runs linear models on the genes to detect differentially expressed
#' genes. Then on a subset of differentially expressed genes, linear models with spherical
#' and exponential covariance structure are run using the x,y coordinates of each spot/sample.
#' The function outputs the Akaike (AIC) and Bayesian Information Criterion (BIC) for
#' the spherical and exponential models to aid in selection of differentially expressed
#' genes under the spatial model with the best fit.
#'
#' @param x an STlist
#' @param samples an integer indicating the spatial samples to be included in the DE tests.
#' Numbers follow the order in `names(x@counts)`. Sample names are also allowed.
#' If NULL, performs tests in all samples.
#' @param annot a column name in `x@spatial_meta` containing the groups/clusters to
#' be tested.
#' @param ws the spatial weight used in STclust
#' @param ks the k value used in STclust, or `dtc` (default) for dynamic tree cut.
#' @param deepSplit the deepSplit value if used in STclust.
#' @param topgenes an integer indicating the top variable genes to select from each sample
#' based on variance (default=2000).
#' @param pval_thr Cut off of adjusted p-values to define differentially expressed genes from
#' non-spatial linear models. Default=0.05.
#' @param pval_adj Method to adjust p-values from non-spatial models. Defaults to
#' `bonferroni`. Other options as available from `p.adjust`.
#' @param sp_topgenes Proportion of differentially expressed genes from non-spatial
#' linear models (and controlled by `pval_thr`) to use in differential gene expression
#' analysis with spatial linear models. Default=0.2.
#' @param cores Number of cores to use in parallelization. If `NULL`, the number of
#' cores is to use is detected automatically.
#' @return a data frame with results of differential gene expression analysis.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # cluster_p <- plot_STclusters(melanoma, samples=c(2,3), visium=F)
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom nlme lme.formula
#
#
STde = function(x=NULL, samples=NULL, annot=NULL, ws=NULL, ks='dtc', deepSplit=NULL,
                topgenes=2000, pval_thr=0.05, pval_adj=NULL, sp_topgenes=0.2, cores=NULL){
  # Record time
  zero_t = Sys.time()

  # Stop function if topgenes and sp_topgenes are not valid
  if((round(topgenes, 0) <= 0) | (sp_topgenes < 0 | sp_topgenes > 1)){
    stop('topgenes or sp_topgenes do not contain valid values.')
  }

  # Define which samples to plot
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

  # Define columns to plot
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
      stop('Specify one or several k values to test, or use ks=\'dtc\' (default).')
    }
  }

  # Check that the meta data column exists
  if(length(grep(paste0('^', annot, '$', collapse='|'),  colnames(x@spatial_meta[[1]]))) == 0){
    stop('No metadata column or clustering parameters were specified. Or specified parameters do not exist in metadata.')
  } else if(length(grep(paste0('^', annot, '$', collapse='|'),  colnames(x@spatial_meta[[1]]))) > 1){
    stop('Only one clustering solution or metadata column can be tested at a time.')
  }

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

  # Get top variable genes
  genes = list()
  for(sample_name in samples){
    genes[[sample_name]] = append(genes[[sample_name]],
                                  x@gene_meta[[sample_name]] %>%
                                    dplyr::arrange(desc(vst.variance.standardized)) %>%
                                    dplyr::slice_head(n=topgenes) %>%
                                    dplyr::select(gene) %>%
                                    unlist() %>%
                                    as.vector()
    )
  }
  rm(topgenes) # Clean environment

  # Print metadata assignment to be tested
  if(grepl('stclust_', annot[1])){
    spw_print = stringr::str_extract(annot[1], '(?<=spw)0\\.[0-9]+')
    if(grepl('_dspl', annot[1])){
      cut_print = stringr::str_extract(annot[1], '(?<=_dspl)[0-9]') %>% paste0('dtc deepSplit=', ., ')...')
    } else{
      cut_print = stringr::str_extract(annot[1], '(?<=_k)[0-9]') %>% paste0('k=', ., ')...')
    }
    test_print = paste0('Testing STclust assignment (w=', spw_print, ', ', cut_print)
    rm(spw_print, cut_print) # Clean environment
  } else{
    test_print = paste0('Testing metadata: ', annot[1], '...')
  }
  cat(crayon::yellow(paste0(test_print, '\n')))
  rm(test_print) # Clean environment

  # Create data frames for each sample with expression and coordinate data
  expr_dfs = list()
  for(sample_name in samples){
    # Add group dummy column and select relevant columns
    # Add character to annotation column that it's treated a factor in case it's an integer
    expr_dfs[[sample_name]] = x@tr_counts[[sample_name]] %>%
      expandSparse() %>%
      tibble::rownames_to_column(var='gene') %>%
      dplyr::filter(gene %in% genes[[sample_name]]) %>%
      tibble::column_to_rownames(var='gene') %>%
      t() %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='libname') %>%
      dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                          tibble::add_column(group=1, .after='libname') %>%
                          dplyr::select(libname, group, ypos, xpos, meta:=!!annot[1]),. , by='libname') %>%
      dplyr::mutate(meta=paste0('m', meta)) %>%
      tibble::column_to_rownames(var='libname')
  }

  # Run non-spatial models
  start_t = Sys.time()
  cat(crayon::yellow('Running non-spatial tests...\n'))
  non_sp_de = list()
  non_sp_models = list()
  for(sample_name in samples){
    # Create a table with combinations of samples, genes, and annotations to test using parallelization
    combo = expand.grid(meta=unique(expr_dfs[[sample_name]][['meta']]),
                        gene=colnames(expr_dfs[[sample_name]] %>%
                                        dplyr::select(-group, -ypos, -xpos, -meta)),
                        stringsAsFactors=F) %>%
      dplyr::arrange(meta)

    # Run models in parallel and get DE results
    non_sp_de[[sample_name]] = non_spatial_de(expr_df=expr_dfs[[sample_name]], combo=combo, pval_adj=pval_adj)

    # Save models
    models_keep = non_sp_de[[sample_name]][[1]] %>%
      dplyr::filter(adj_p_val < pval_thr) %>%
      dplyr::arrange(adj_p_val) %>%
      dplyr::slice_head(prop=sp_topgenes) %>%
      dplyr::select(gene, cluster) %>%
      tidyr::unite('x', gene, cluster) %>%
      unlist() %>%
      as.vector()
    non_sp_models[[sample_name]] = non_sp_de[[sample_name]][[2]][models_keep]
    rm(models_keep) # Clean environment

    # Extract DE results and add column with sample ID
    non_sp_de[[sample_name]] = non_sp_de[[sample_name]][[1]] %>%
      tibble::add_column(sample=sample_name, .before=1)
  }
  end_t = difftime(Sys.time(), start_t, units='min')
  cat(crayon::green(paste0('Non-spatial tests completed (', round(end_t, 2), ' min).\n')))

  # Run spatial tests on subset of significant genes
  start_t = Sys.time()
  if(sp_topgenes > 0){
    cat(crayon::yellow('Running spatial tests...\n'))
    sp_models = list()
    # Run models
    for(sample_name in samples){
      sp_models[[sample_name]] = spatial_de(expr_df=expr_dfs[[sample_name]], non_sp_mods=non_sp_models[[sample_name]], sp_topgenes=sp_topgenes)
    }
  }
  rm(expr_dfs) # Clean environment
  end_t = difftime(Sys.time(), start_t, units='min')
  cat(crayon::green(paste0('Spatial tests completed (', round(end_t, 2), ' min).\n')))

  sp_de = tibble::tibble()
  for(sample_name in samples){
    # Compile results of spatial models and merge with non-spatial results
    for(i in names(sp_models[[sample_name]])){
      # Check that model results exists and extract values, gene names, and clusters tested
      if(sp_models[[sample_name]][[i]][['sph']][[1]][[1]] != 'no_conv'){
        sph_summ = summary(sp_models[[sample_name]][[i]][['sph']])
        meta_sph_tmp = names(sp_models[[sample_name]][[i]][['sph']]$coefficients$fixed) %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_sph_tmp = gsub(paste0('_', meta_sph_tmp, '$'), '', i)
      } else{
        sph_summ = list(AIC=-9999, BIC=-9999, tTable=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
        meta_sph_tmp = sp_models[[sample_name]][[i]][['sph']][[2]] %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_sph_tmp = gsub(paste0('_', meta_sph_tmp, '$'), '', i)
      }
      if(sp_models[[sample_name]][[i]][['exp']][[1]][[1]] != 'no_conv'){
        exp_summ = summary(sp_models[[sample_name]][[i]][['exp']])
        meta_exp_tmp = names(sp_models[[sample_name]][[i]][['exp']]$coefficients$fixed) %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_exp_tmp = gsub(paste0('_', meta_exp_tmp, '$'), '', i)
      } else{
        exp_summ = list(AIC=-9999, BIC=-9999, tTable=data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999)))
        meta_exp_tmp = sp_models[[sample_name]][[i]][['exp']][[2]] %>% grep('metam', ., value=T) %>% gsub('metam', '', .)
        gene_exp_tmp = gsub(paste0('_', meta_exp_tmp, '$'), '', i)
      }

      # Make sure genes and clusters from both exponential and spherical models match
      if(meta_exp_tmp == meta_sph_tmp & gene_exp_tmp == gene_sph_tmp){
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
    # Merge spatial and non-spatial DE results
    non_sp_de[[sample_name]] = dplyr::full_join(non_sp_de[[sample_name]], sp_de, by=c('sample', 'gene', 'cluster'))
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  cat(crayon::green(paste0('STde completed in ', round(end_t, 2), ' min.\n')))

  return(non_sp_de)
}


# Helpers ----------------------------------------------------------------------

##
# non_spatial_de
# Run non-spatial linear models for a set of gene x cluster combinations
# @param expr_df a data frame (rows=spots/cells) with normalized expression data
# as well as coordinates (x,y), and a group variable for random effects.
# @param combo a data frame with two columns (meta andd gene) containing all
# combinations of gene x cluster to test.
# @param pval_adj a character string with the method for p-value adjustment. Passed
# to p.adjust.
# @return a list containing a dataframe wit DE results and a list of lme models
#
non_spatial_de = function(expr_df=NULL, combo=NULL, pval_adj=NULL){
  # Paralellize lme models
  if(is.null(cores)){
    cores = count_cores(length(non_sp_mods))
  }
  result_ls = parallel::mclapply(seq_along(1:nrow(combo)), function(i){
    res_df = list()
    meta_tmp = as.vector(unlist(combo[i, 1]))
    gene_tmp = as.vector(unlist(combo[i, 2]))
    expr_df = expr_df %>%
      dplyr::select(group, meta, !!gene_tmp) %>%
      dplyr::mutate(meta=factor(dplyr::case_when(meta != meta_tmp ~ 'other',
                                                 TRUE ~ meta), levels=c('other', meta_tmp)))

    # Workaround to avoid dashes in gene names as nlme does not handle well
    genenodash_tmp = gsub('\\-', '_', gene_tmp)
    colnames(expr_df) = gsub('\\-', '_', colnames(expr_df))

    # Create non-spatial model
    res_df[[gene_tmp]] = nlme::lme(fixed=as.formula(paste0("`", genenodash_tmp, "`", '~meta')), data=expr_df, random=~1|group, method='REML')
    avgexpr_ref = mean(expr_df %>% dplyr::filter(meta != 'other') %>% dplyr::select(!!genenodash_tmp) %>% unlist())
    avgexpr_other = mean(expr_df %>% dplyr::filter(meta == 'other') %>% dplyr::select(!!genenodash_tmp) %>% unlist())
    res_df[['test']] = tibble::tibble(gene=gene_tmp,
                                      avc_log2fc=(avgexpr_ref - avgexpr_other),
                                      cluster=gsub('metam', '', rownames(summary(res_df[[gene_tmp]])[['tTable']])[2]),
                                      p_val=summary(res_df[[gene_tmp]])[['tTable']][2, 5])
    return(res_df)
  }, mc.cores=cores)
  rm(expr_df) # Clean environment

  # Collapse in single data frame
  if(nrow(combo) == length(result_ls)){
    result_df = tibble::tibble()
    models_ls = list()
    for(i in 1:nrow(combo)){
      result_df = rbind(result_df, result_ls[[i]][['test']])
      # Keep model
      models_ls[[paste0(result_ls[[i]][['test']][1, c('gene', 'cluster')], collapse='_')]] = result_ls[[i]][[1]]
    }
    if(is.null(pval_adj)){
      pval_adj = 'bonferroni'
    }
    result_df[['adj_p_val']] = p.adjust(result_df[['p_val']], method=pval_adj)
    result_ls = result_ls[['models']]
  }
  return(list(result_df, models_ls))
}


##
# spatial_de
# Run spatial linear models with exponential and spherical covariance structures
# @param expr_df a data frame (rows=spots/cells) with normalized expression data
# as well as coordinates (x,y), and a group variable for random effects.
# @param non_sp_mods a list of non-spatial models to update with covariance structures
# @param sp_topgenes an integer indicating the number of top DE genes to apply spatial models
# @return a list of spatial models
#
spatial_de = function(expr_df=NULL, non_sp_mods=NULL, sp_topgenes=NULL){
  # Paralellize lme models
  if(is.null(cores)){
    cores = count_cores(length(non_sp_mods))
  }
  sp_mods = parallel::mclapply(names(non_sp_mods), function(i){
    res_ls = list() # Initiates list to storage results
    # Extract cluster and gene identifiers
    meta_tmp = rownames(summary(non_sp_mods[[i]])[['tTable']])[2] %>%
      gsub('metam', '', .)
    gene_tmp = i %>%
      gsub(paste0('_', meta_tmp), '', .)
    # Workaround to avoid dashes in gene names as nlme does not handle well
    genenodash_tmp = gsub('\\-', '_', gene_tmp)

    expr_tmp = expr_df %>%
      dplyr::select(group, ypos, xpos, meta, !!gene_tmp) %>%
      dplyr::mutate(meta=factor(dplyr::case_when(meta != paste0('m', meta_tmp) ~ 'other',
                                                 TRUE ~ meta), levels=c('other', paste0('m', meta_tmp))))

    colnames(expr_tmp) = gsub('\\-', '_', colnames(expr_tmp))

    # Rewrite the lme command call so that it works on the new environment
    non_sp_mods[[i]][['call']] = str2lang('nlme::lme.formula(fixed = as.formula(paste0("`", genenodash_tmp, "`", "~meta")), data = expr_tmp, random = ~1 | group, method = "REML")')

    # Run spherical model and catch error if no convergence
    sph_out = tryCatch({
      update(non_sp_mods[[i]], correlation=nlme::corSpher(form=~xpos+ypos, nugget=T, metric='euclidean'), method='REML', control=nlme::lmeControl(returnObject=F, msVerbose=F, gradHess=F))
    }, error=function(err){return(err)})

    if(class(sph_out) == 'simpleError'){
      res_ls[['sph']][[1]] = 'no_conv'
      res_ls[['sph']][[2]] = names(non_sp_mods[[i]]$coefficients$fixed)
    } else{
      res_ls[['sph']] = sph_out
    }

    # Run exponential model and catch error if no convergence
    exp_out = tryCatch({
      update(non_sp_mods[[i]], correlation=nlme::corExp(form=~xpos+ypos, nugget=T, metric='euclidean'), method='REML', control=nlme::lmeControl(returnObject=F, msVerbose=F, gradHess=F))
    }, error=function(err){return(err)})

    if(class(exp_out) == 'simpleError'){
      res_ls[['exp']][[1]] = 'no_conv'
      res_ls[['exp']][[2]] = names(non_sp_mods[[i]]$coefficients$fixed)
    } else{
      res_ls[['exp']] = exp_out
    }
    return(res_ls)
  }, mc.cores=cores)

  if(length(non_sp_mods) == length(sp_mods)){
    names(sp_mods) = names(non_sp_mods) # Transfer names of non-spatial list to  spatial list
  } else{
    stop('Looks like the number of non-spatial models do not match the number of spatial models')
  }
  return(sp_mods)
}

