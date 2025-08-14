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
#' using the same parameters (`w`, `k`, `deepSplit`). Otherwise, the assignments are
#' specified by indicating one of the column names in `x@spatial_meta`. The function
#' uses `spaMM::fitme` and is computationally expensive even on HPC environments.
#' To run the STdiff using the non-spatial approach (faster), set `sp_topgenes=0`.
#'
#' @param x an STlist
#' @param samples an integer indicating the spatial samples to be included in the DE tests.
#' Numbers follow the order in `names(x@counts)`. Sample names are also allowed.
#' If NULL, performs tests on all samples
#' @param annot a column name in `x@spatial_meta` containing the groups/clusters to
#' be tested. Required if `k` and `w` are empty.
#' @param w the spatial weight used in STclust. Required if `annot` is empty.
#' @param k the k value used in STclust, or `dtc` for dynamicTreeCut clusters. Required if `annot` is empty.
#' @param deepSplit the deepSplit value if used in STclust. Required if `k='dtc'`.
#' @param topgenes an integer indicating the top variable genes to select from each sample
#' based on variance (default=5000). If NULL, all genes are selected.
#' @param pval_thr cut-off of adjusted p-values to define differentially expressed genes from
#' non-spatial linear models. A proportion of genes (`sp_topgenes`) under this cut-off
#' will be applied the spatial models. Default=0.05
#' @param pval_adj Method to adjust p-values. Defaults to `FDR`. Other options as
#' available from `p.adjust`
#' @param test_type one of `mm`, `t_test`, or `wilcoxon`. Specifies the type of
#' test performed.
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
#' @param verbose either logical or an integer (0, 1, or 2) to increase verbosity
#' @return a list with one data frame per sample with results of differential gene
#' expression analysis
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom rlang :=
#' @importFrom stats p.adjust
#
#
STdiff = function(x=NULL, samples=NULL, annot=NULL, w=NULL, k=NULL, deepSplit=NULL,
                  topgenes=5000, pval_thr=0.05, pval_adj='fdr', test_type='mm', sp_topgenes=0.2,
                  clusters=NULL, pairwise=FALSE, verbose=1L, cores=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()

  # Convert user inputs to expected inputs
  topgenes = as.integer(ceiling(topgenes))
  pval_thr = as.double(pval_thr)
  sp_topgenes = as.double(sp_topgenes)

  # Set verbosity level
  verbose = as.integer(verbose)
  if(!is.integer(verbose) | !(verbose %in% c(0L, 1L, 2L))){
    verbose = 1L
  }

  # Check that at least two clusters have been specified if pairwise
  if(pairwise & !is.null(clusters)){
    if(length(clusters) < 2){
      stop('If pairwise tests requested, at least two clusters are required.')
    }
  }

  # Stop function if topgenes and sp_topgenes are not valid
  if((round(topgenes, 0) <= 0) | (sp_topgenes < 0 | sp_topgenes > 1)){
    stop('topgenes or sp_topgenes contain invalid values.')
  }

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

  # Define column to test
  if(is.null(annot)){
    # If annot=NULL, must specify w and k
    if(!is.null(w) & !is.null(k)){
      if(length(w) != 1 | length(k) != 1){
        stop('Please specify a single value for w and a single value for k.')
      }
      # Construct column name from w and k
      annot = paste0('stclust_spw', as.character(w))
      if(k == 'dtc'){
        if(is.null(deepSplit)){
          stop('If k=\"dtc\", then specify deepSplit.')
        } else if(is.logical(deepSplit)){
          if(deepSplit){
            annot = paste0(annot, '_dsplTrue')
          } else{
            annot = paste0(annot, '_dsplFalse')
          }
        } else{
          annot = paste0(annot, '_dspl', deepSplit)
        }
      } else{
        if(!is.numeric(k)){
          stop('Specify a valid k value.')
        }
        annot = paste0(annot, '_k', as.character(k))
      }
    } else{
      stop('If no specific annotation is specified, please specify both w and k (STclust parameters).')
    }
  } else{
    if(length(annot) == 0){
      stop('If w and k are not specified, one annotation column from @spatial_meta should be specified.')
    } else if(length(annot) > 1){
      stop('Only one annotation column from @spatial_meta can be tested at a time. ')
    }
  }

  # Check that the meta data column exists (if not in specific sample, remove sample)
  samples_tmp = samples
  for(i in samples){
    if(!(annot %in% colnames(x@spatial_meta[[i]]))){
      samples_tmp = grep(i, samples_tmp, value=T, invert=T)
      if(verbose){
        cat(paste0('Skipping ', i, '. Annotation not available for this sample.\n'))
      }
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
      #x@gene_meta[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
      x@gene_meta[[i]] = Seurat_FindVariableFeatures(x@counts[[i]]) %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::select('gene', 'vst.variance.standardized') %>%
        dplyr::left_join(x@gene_meta[[i]], ., by='gene')
    }
  }

  # Print metadata assignment to be tested
  if(grepl('stclust_', annot[1])){
    spw_print = stringr::str_extract(annot[1], '(?<=spw)0\\.?[0-9]*')
    if(grepl('_dspl', annot[1])){
      cut_print = stringr::str_extract(annot[1], '(?<=_dspl)[\\.0-9FalseTrue]*') %>% paste0('dtc deepSplit=', ., ')...\n')
    } else{
      cut_print = stringr::str_extract(annot[1], '(?<=_k)[0-9]') %>% paste0('k=', ., ')...\n')
    }
    test_print = paste0('Testing STclust assignment (w=', spw_print, ', ', cut_print)
    rm(spw_print, cut_print) # Clean environment
  } else{
    test_print = paste0('Testing metadata: ', annot[1], '...\n')
  }

  # Print metadata to be tested
  if(verbose){
    cat(test_print)
  }
  rm(test_print) # Clean environment

  # Extract spot/cell annotations for each sample
  metas = tibble::tibble()
  for(sample_name in samples){
    meta_tmp = as.character(unique(x@spatial_meta[[sample_name]][[annot]]))
    metas = dplyr::bind_rows(metas, tibble::tibble(samplename=sample_name, orig_annot=meta_tmp))
    rm(meta_tmp) # Clean env
  }

  # Find variable genes for each sample
  genes = tibble::tibble()
  for(sample_name in samples){
    # Get top variable genes
    genes_tmp = x@gene_meta[[sample_name]] %>% dplyr::arrange(dplyr::desc(vst.variance.standardized))
    if(!is.null(topgenes)){
      genes_tmp = genes_tmp %>% dplyr::slice_head(n=topgenes)
    }
    genes_tmp = as.vector(genes_tmp[['gene']])
    genes = dplyr::bind_rows(genes, tibble::tibble(samplename=sample_name, gene=genes_tmp))
    rm(genes_tmp) # Clean env
  }

  # Merge sample names, genes, and annotations in the same data frame
  gene_and_meta = dplyr::full_join(genes, metas, by='samplename', multiple='all')
  rm(genes, metas) # Clean env

  # Create "dictionary" with coded annotations (to avoid potentially problematic characters)
  meta_dict = tibble::tibble(orig_annot=unique(gene_and_meta[['orig_annot']]),
                             coded_annot=paste0('c', seq(1, length(unique(gene_and_meta[['orig_annot']])))))

  # Append recoded annotations in gene_and_meta
  gene_and_meta = gene_and_meta %>% dplyr::left_join(., meta_dict, by='orig_annot') %>% dplyr::rename(meta=coded_annot)

  # Create a table with unique combinations of genes and annotations to test using parallelization
  combo_df = prepare_stdiff_combo(to_expand=gene_and_meta, user_clusters=clusters, pairwise=pairwise, verbose=verbose)
  rm(gene_and_meta) # Clean env

  ######## ######## ######## BEGIN NON-SPATIAL TESTS ######## ######## ########
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

  # Get annotations to run in parallel
  # clusters_tmp = unique(combo_df[['meta1']])
  # if(pairwise){
  #   clusters_tmp = unique(append(clusters_tmp, unique(combo_df[['meta2']])))
  # }

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
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(unique(combo_df[['samplename']])))
  } else{
    cores = ceiling(cores)
  }

  # Paralellize non-spatial tests
  non_sp_models = parallel::mclapply(1:length(unique(combo_df[['samplename']])), function(i){
    # Get sample name
    sample_name = unique(combo_df[['samplename']])[i]
    # Subset combo_df to those of a given sample
    combo_clust = combo_df[ combo_df[['samplename']] == sample_name, ]

    # Create data frame with expression, coordinate, and cluster data
    # Add group dummy column (in case of mixed models was requested)
    expr_df = expandSparse(x@tr_counts[[sample_name]])
    expr_df = expr_df[ rownames(expr_df) %in% unique(combo_clust[['gene']]), ]
    expr_df = t(expr_df) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='libname') %>%
      dplyr::right_join(x@spatial_meta[[sample_name]] %>%
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
  }, mc.cores=cores)

  # Sometimes the mixed models option produces a warning
  # <simpleWarning in .check_y(family, y, BinomialDen): var(response) = 0, which may cause errors.>
  # if(any(names(non_sp_models) == 'warning')){
  #   non_sp_models = non_sp_models[['value']]
  # }

  # Put sample names to results of parallelezation
  names(non_sp_models) = unique(combo_df[['samplename']])

  # Summarize DE analyses
  result_de = tibble::tibble()
  for(i in names(non_sp_models)){
    for(mod in names(non_sp_models[[i]])){
      sample_tmp = non_sp_models[[i]][[mod]][['samplename']]
      gene_tmp = non_sp_models[[i]][[mod]][['gene']]
      avglfc_tmp = non_sp_models[[i]][[mod]][['avglfc']]
      cluster_1_tmp = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == non_sp_models[[i]][[mod]][['meta1']] ]))

      # Compile results in a row
      df_tmp = tibble::tibble(sample=sample_tmp,
                              gene=gene_tmp,
                              avg_log2fc=avglfc_tmp,
                              cluster_1=cluster_1_tmp)

      if(pairwise){
        #if(test_type == 'mm'){
        # cluster_2_tmp = unique(non_sp_models[[i]][[mod]][['nonspmod']][['data']][['meta']]) %>% grep(i, ., value=T, invert=T)
        #cluster_2_tmp = non_sp_models[[i]][[mod]][['meta2']]
        #} else if(test_type == 't_test' | test_type == 'wilcoxon'){
        cluster_2_tmp = non_sp_models[[i]][[mod]][['meta2']]
        #}
        cluster_2_tmp = as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == cluster_2_tmp ]))
        df_tmp[['cluster_2']] = cluster_2_tmp
        rm(cluster_2_tmp) # Clean env
      }

      if(test_type == 'mm'){
        # df_tmp[['mm_p_val']] = spaMM::summary.HLfit(non_sp_models[[i]][[mod]][['nonspmod']],
        #                                             details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]
        df_tmp[['mm_p_val']] = non_sp_models[[i]][[mod]][['nonspmod']]
      } else if(test_type == 't_test'){
        df_tmp[['ttest_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      } else if(test_type == 'wilcoxon'){
        df_tmp[['wilcox_p_val']] = non_sp_models[[i]][[mod]][['mean_test']][['p.value']]
      }

      # Add results to list of data frames
      result_de = dplyr::bind_rows(result_de, df_tmp)

      rm(sample_tmp, gene_tmp, avglfc_tmp, cluster_1_tmp, df_tmp) # Clean env
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
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, cluster_2, adj_p_val, dplyr::desc(avg_log2fc))
    } else{
      non_sp_de_tmp[[i]] = non_sp_de_tmp[[i]] %>% dplyr::arrange(cluster_1, adj_p_val, dplyr::desc(avg_log2fc))
    }
  }
  result_de = non_sp_de_tmp
  rm(non_sp_de_tmp) # Clean env

  end_t = difftime(Sys.time(), start_t, units='min')
  if(verbose){
    cat(paste0('\tCompleted ', test_print, ' (', round(end_t, 2), ' min).\n'))
  }

  ######## ######## ######## BEGIN SPATIAL TESTS ######## ######## ########     ##### MADE BIG CHANGES TO NON-SPATIAL... NEED TO CHECK SPATIAL TESTS
  if(test_type == 'mm' & sp_topgenes > 0){
    start_t = Sys.time()
    if(verbose){
      cat('\tRunning spatial tests...\n')
    }
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
        expr_df = expandSparse(x@tr_counts[[sample_name]])
        expr_df = t(expr_df) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var='libname') %>%
          dplyr::right_join(x@spatial_meta[[sample_name]] %>%
                              tibble::add_column(group=1, .after='libname') %>%
                              dplyr::select(c('libname', 'group', 'ypos', 'xpos'), orig_annot:=!!annot),. , by='libname') %>%
          dplyr::left_join(., meta_dict, by='orig_annot') %>%
          dplyr::rename(meta=coded_annot) %>%
          tibble::column_to_rownames(var='libname') %>%
          dplyr::relocate(meta, .before=1) %>%
          dplyr::select(-c('orig_annot'))

        sp_models[[sample_name]] = parallel::mclapply(1:length(clusters_tmp), function(i){ ####### MCLAPPLY INSTEAD OF PBMCLAPPLY (PBMCLAPPY MAY HAVE ISSUES WITH TRYCATCH)
          #sp_models[[sample_name]] = bettermc::mclapply(1:length(clusters_tmp), function(i){  ####### TRY BETTERMC INSTEAD OF PARALLEL
          #sp_models[[sample_name]] = pbmcapply::pbmclapply(1:length(clusters_tmp), function(i){ ####### TEST LOOP INSTEAD MCLAPPLY
          #for(i in 1:length(clusters_tmp)){  ####### TEST LOOP INSTEAD MCLAPPLY
          # Non-spatial models to update
          models_keep_tmp = models_keep %>%
            dplyr::filter(cluster_1 == clusters_tmp[i]) %>%
            dplyr::select('genexcluster') %>%
            unlist() %>% as.vector()
          non_sp_models_tmp = non_sp_models[[sample_name]][ models_keep_tmp ]

          if(verbose == 1L){
            if(pairwise){
              stdout_print = models_keep_tmp %>%
                stringr::str_extract(., paste0('(?<=', clusters_tmp[i], '_)c[0-9]+$')) %>%
                stringr::str_split(., '_') %>% unlist() %>% unique()
              stdout_print_tmp = c()
              for(m2 in stdout_print){
                stdout_print_tmp = append(stdout_print_tmp,
                                          meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == m2 ])
              }
              stdout_print = stdout_print_tmp %>% paste0(., collapse=', ')
              system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                                 sample_name, ", ",
                                                 meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == clusters_tmp[i] ],
                                                 " vs. ", stdout_print,
                                                 " (", length(models_keep_tmp), " tests)...")))
              rm(stdout_print)
            } else {
              system(sprintf('echo "%s"', paste0("\t\tSample: ",
                                                 sample_name, ", ",
                                                 meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == clusters_tmp[i] ],
                                                 " (", length(models_keep_tmp), " tests)...")))
            }
          }

          # Run models
          sp_mods = spatial_de(expr_dat=expr_df, non_sp_mods=non_sp_models_tmp, annot_dict=meta_dict, verb=verbose)

          return(sp_mods) ####### TEST LOOP INSTEAD MCLAPPLY
        }, mc.cores=cores) ####### TEST LOOP INSTEAD MCLAPPLY
        #sp_models[[sample_name]] = sp_mods  ####### TEST LOOP INSTEAD MCLAPPLY
        #} ####### TEST LOOP INSTEAD MCLAPPLY
        rm(models_keep, clusters_tmp) # Clean environment
      } else{
        warning(paste0('No genes to test for ',  sample_name, '. Try increasing sp_topgenes.'))
      }
    }
    end_t = difftime(Sys.time(), start_t, units='min')
    if(verbose){
      cat(paste0('\n\tCompleted spatial mixed models (', round(end_t, 2), ' min).\n'))
    }
    rm(non_sp_models) # Clean environment

    # Compile results of spatial models and merge with non-spatial results
    sp_de = tibble::tibble()
    for(i in names(sp_models)){
      for(cl in 1:length(sp_models[[i]])){
        for(test_name in names(sp_models[[i]][[cl]])){
          # Check that model results exists and extract values, gene names, and clusters tested
          if(any(class(sp_models[[i]][[cl]][[test_name]][['spmod']]) == 'HLfit')){
            exp_summ = spaMM::summary.HLfit(sp_models[[i]][[cl]][[test_name]][['spmod']], details=c(p_value='Wald'), verbose=F)
            exp_summ = exp_summ[['beta_table']]
            exp_warn = NA_character_
          } else if(sp_models[[i]][[cl]][[test_name]][['spmod']] == 'time_out'){
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'time_limit'
          } else if(sp_models[[i]][[cl]][[test_name]][['spmod']] == 'no_conv'){
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'no_convergence'
          } else{
            exp_summ = data.frame(c(NA,NA), c(NA,NA), c(NA,NA), c(NA, -9999))
            exp_warn = 'unknown_na'
          }
          meta1_exp_tmp = stringr::str_replace(test_name, paste0(i, '_', sp_models[[i]][[cl]][[test_name]][['gene']], '_'), '') %>%
            stringr::str_replace(., '_other$', '')
          if(pairwise){
            meta2_exp_tmp = stringr::str_extract(meta1_exp_tmp, 'c[0-9]+$')
            meta1_exp_tmp = stringr::str_extract(meta1_exp_tmp, '^c[0-9]+')
          }
          # Prepare table with results
          # Replace original annotations
          tibble_tmp = tibble::tibble(sample=i,
                                      gene=sp_models[[i]][[cl]][[test_name]][['gene']],
                                      cluster_1=as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta1_exp_tmp ])),
                                      exp_p_val=exp_summ[[2, 4]],
                                      comments_spatial=exp_warn)
          if(pairwise){
            tibble_tmp = tibble_tmp %>%
              tibble::add_column(cluster_2=as.vector(unlist(meta_dict[['orig_annot']][ meta_dict[['coded_annot']] == meta2_exp_tmp ])), .after='cluster_1')
          }
          sp_de = dplyr::bind_rows(sp_de, tibble_tmp)

          rm(list=grep('meta1_exp_tmp|meta2_exp_tmp|exp_summ|exp_warn', ls(), value=T))
        }
      }
    }

    # Merge spatial and non-spatial DE results
    for(sample_name in names(result_de)){
      if(pairwise){
        # Adjust spatial p-values
        tibble_tmp = sp_de %>%
          dplyr::select(-c('exp_p_val')) %>%
          dplyr::filter(sample == sample_name) %>%
          dplyr::left_join(., sp_de %>%
                             dplyr::select(-c('comments_spatial')) %>%
                             dplyr::filter(sample == sample_name) %>%
                             dplyr::filter(!is.na(exp_p_val) & exp_p_val != -9999) %>%
                             tibble::add_column(exp_adj_p_val=p.adjust(.[['exp_p_val']], method=pval_adj)), by=c('sample', 'gene', 'cluster_1', 'cluster_2')) %>%
          dplyr::relocate(comments_spatial, .after='exp_adj_p_val')
        result_de[[sample_name]] = dplyr::full_join(result_de[[sample_name]], tibble_tmp, by=c('sample', 'gene', 'cluster_1', 'cluster_2'))
        rm(tibble_tmp) # Clean env
      } else{
        tibble_tmp = sp_de %>%
          dplyr::select(-c('exp_p_val')) %>%
          dplyr::filter(sample == sample_name) %>%
          dplyr::left_join(., sp_de %>%
                             dplyr::select(-c('comments_spatial')) %>%
                             dplyr::filter(sample == sample_name) %>%
                             dplyr::filter(!is.na(exp_p_val) & exp_p_val != -9999) %>%
                             tibble::add_column(exp_adj_p_val=p.adjust(.[['exp_p_val']], method=pval_adj)), by=c('sample', 'gene', 'cluster_1')) %>%
          dplyr::relocate(comments_spatial, .after='exp_adj_p_val')
        result_de[[sample_name]] = dplyr::full_join(result_de[[sample_name]], tibble_tmp, by=c('sample', 'gene', 'cluster_1'))
        rm(tibble_tmp) # Clean env
      }
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::mutate(comments=dplyr::case_when(is.na(comments) ~ comments_spatial, TRUE ~ comments)) %>%
        dplyr::select(-c('comments_spatial')) %>%
        dplyr::relocate(comments, .after='exp_adj_p_val')

      result_de[[sample_name]] = result_de[[sample_name]] %>%
        #dplyr::mutate(exp_adj_p_val=dplyr::case_when(exp_p_val == -9999 ~ -9999, TRUE ~ exp_adj_p_val)) %>%
        dplyr::mutate(spatialmod=dplyr::case_when(is.na(exp_p_val) ~ '2', TRUE ~ '1')) # To put on top of table the genes with spatial models

      if(pairwise){
        result_de[[sample_name]] = result_de[[sample_name]] %>%
          dplyr::arrange(spatialmod, cluster_1, cluster_2, exp_adj_p_val)
      } else{
        result_de[[sample_name]] = result_de[[sample_name]] %>%
          dplyr::arrange(spatialmod, cluster_1, exp_adj_p_val)
      }
      result_de[[sample_name]] = result_de[[sample_name]] %>%
        dplyr::select(-c('spatialmod'))

      #rm(tibble_tmp) # Clean env
    }
    rm(sp_de) # Clean environment

  }   ##### CLOSES IF TO ENTER SPATIAL TESTS CODE

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STdiff completed in ', round(end_t, 2), ' min.\n'))
  }
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

    # Calculate p-value
    res_mod = spaMM::summary.HLfit(res_mod, details=c(p_value="Wald"), verbose=F)[['beta_table']][2, 4]

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
# prepare_stdiff_combo
# @return data frame with combinations
#
#prepare_stdiff_combo = function(to_expand=NULL, clusters=NULL, annot_dict=NULL, pairwise=NULL){
prepare_stdiff_combo = function(to_expand=NULL, user_clusters=NULL, pairwise=NULL, verbose=NULL){
  combo_df = tibble::tibble()
  for(sample_name in unique(to_expand[['samplename']])){
    # Extract annotations for a given sample
    to_expand_subset = to_expand %>% dplyr::filter(samplename == sample_name) #%>%
    #  dplyr::select('meta') %>% unlist() %>% unique()

    # Define clusters to test if NULL, or subset to those requested by user
    if(is.null(user_clusters)){
      annots_tmp = to_expand_subset %>% dplyr::select('meta') %>% unlist() %>% unique()
    } else{
      annots_tmp = to_expand_subset %>% dplyr::filter(orig_annot %in% user_clusters) %>% dplyr::select('meta') %>% unlist() %>% unique()
    }
    #annots_tmp = annots_tmp[annots_tmp %in% unlist(annot_dict[['coded_annot']][ annot_dict[['orig_annot']] %in% clusters_tmp ]) ]
    if(length(annots_tmp) >= 2){
      if(pairwise){
        # Get unique combinations of each two annotation, without repeating (e.g., c1 vs c2 and c2 vs c1)
        # Get combinations for each cluster
        combo_meta = tibble::tibble()
        annots_tmp2 = annots_tmp
        for(cl_tmp in annots_tmp){
          combo_meta = dplyr::bind_rows(combo_meta,
                                        expand.grid(cl_tmp,
                                                    grep(cl_tmp, annots_tmp2, value=T, invert=T),
                                                    stringsAsFactors=F))
          annots_tmp2 = grep(cl_tmp, annots_tmp2, value=T, invert=T)
        }
        names(combo_meta) = c('meta1', 'meta2')
        rm(cl_tmp, annots_tmp2) #Clean env

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
      rm(to_expand_subset, annots_tmp) # Clean env
      combo_df = dplyr::bind_rows(combo_df,
                                  combo_meta %>%
                                    tibble::add_column(samplename=sample_name, .before=1))
    } else{
      if(verbose >= 1L){
        cat(paste0('\t\tSkipping sample ', sample_name, '. Less than two clusters to compare.\n'))
      }
    }
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
      res_test = stats::t.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
                               expr_tmp %>% dplyr::filter(meta == meta2_tmp) %>% dplyr::select('exprval') %>% unlist())
    } else if(test_type == 'wilcoxon'){
      # Create non-spatial model
      res_test = stats::wilcox.test(expr_tmp %>% dplyr::filter(meta == meta1_tmp) %>% dplyr::select('exprval') %>% unlist(),
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

