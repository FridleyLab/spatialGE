# Assemble stage: compiles per-worker model results into the tables STdiff() returns.

##
# stdiff_assemble_nonspatial_results
# Compile the per-sample non-spatial models into one result table per sample,
# with adjusted p-values.
# @param non_sp_models the per-sample non-spatial models (named by sample), as
# returned by `stdiff_run_nonspatial_tests`.
# @param meta_dict the original-to-coded annotation dictionary.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param test_type one of `mm`, `t_test`, or `wilcoxon`.
# @param pval_adj the p-value adjustment method (passed to `p.adjust`).
# @return a named list (by sample) of DE result tibbles.
#
stdiff_assemble_nonspatial_results = function(non_sp_models, meta_dict, pairwise, test_type, pval_adj){
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

  return(result_de)
}


##
# stdiff_assemble_spatial_results
# Compile the per-sample spatial models into a results table, then merge it
# into (and adjust p-values against) the non-spatial result_de tables.
# @param result_de the per-sample non-spatial DE result tables (already
# labeled with which gene x cluster combinations were spatially tested, by
# `stdiff_run_spatial_tests`).
# @param sp_models the per-sample, per-cluster spatial models, as returned by
# `stdiff_run_spatial_tests`.
# @param meta_dict the original-to-coded annotation dictionary.
# @param pairwise whether tests are pairwise or one-vs-rest.
# @param pval_adj the p-value adjustment method (passed to `p.adjust`).
# @return a named list (by sample) of DE result tibbles with spatial columns merged in.
#
stdiff_assemble_spatial_results = function(result_de, sp_models, meta_dict, pairwise, pval_adj){
  # To prevent NOTES in R CMD check
  . = NULL

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

  return(result_de)
}
