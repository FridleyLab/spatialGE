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

  # Validate, and accept a Seurat object by coercing it (see assert_stlist).
  x = assert_stlist(x)

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
      raise_err(err_code='error0075')
    }
  }

  # Stop function if topgenes and sp_topgenes are not valid
  if((round(topgenes, 0) <= 0) | (sp_topgenes < 0 | sp_topgenes > 1)){
    raise_err(err_code='error0076')
  }

  # Define samples using names (convert indexes to names if necessary)
  samples = resolve_samples(x, samples, slot='spatial_meta', numeric_handling='na_omit',
                            char_transform='exact', post_check='any_in',
                            message='None of the requested samples are present in the STlist.')

  # Define column to test
  if(is.null(annot)){
    # If annot=NULL, must specify w and k
    if(!is.null(w) & !is.null(k)){
      if(length(w) != 1 | length(k) != 1){
        raise_err(err_code='error0077')
      }
      # Construct column name from w and k
      annot = paste0('stclust_spw', as.character(w))
      if(k == 'dtc'){
        if(is.null(deepSplit)){
          raise_err(err_code='error0078')
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
          raise_err(err_code='error0079')
        }
        annot = paste0(annot, '_k', as.character(k))
      }
    } else{
      raise_err(err_code='error0080')
    }
  } else{
    if(length(annot) == 0){
      raise_err(err_code='error0081')
    } else if(length(annot) > 1){
      raise_err(err_code='error0082')
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
      raise_err(err_code='error0083')
    }
  }
  samples = samples_tmp
  rm(samples_tmp) # Clean env

  # Validate the annotation's type and number of groups before anything is fit.
  # STdiff treats the column as group labels and joins on it internally, so a
  # numeric column (a perfectly ordinary metadata column such as total_counts)
  # previously surfaced only much later, from dplyr:
  #   Can't join `x$orig_annot` with `y$orig_annot` due to incompatible types.
  #   i `x$orig_annot` is a <double>
  # which names neither STdiff nor the argument to change.
  for(i in samples){
    annot_vals = x@spatial_meta[[i]][[annot]]
    if(is.numeric(annot_vals)){
      raise_err(err_code='error0070', samplename=i)
    }
    if(length(unique(stats::na.omit(as.character(annot_vals)))) < 2){
      raise_err(err_code='error0071', samplename=i)
    }
  }

  # Calculate standardized variance for all samples (regardless if sample is not in requested test)
  x = stdiff_ensure_vst(x)

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

  # Build table of gene x cluster combinations to test, and the annotation dictionary
  combo_prep = stdiff_build_combo_table(x, samples, annot, topgenes, clusters, pairwise, verbose)
  combo_df = combo_prep[['combo_df']]
  meta_dict = combo_prep[['meta_dict']]
  rm(combo_prep) # Clean env

  ######## ######## ######## BEGIN NON-SPATIAL TESTS ######## ######## ########
  nonspatial_run = stdiff_run_nonspatial_tests(x, combo_df, clusters, meta_dict, annot, pairwise, verbose, test_type, cores)
  non_sp_models = nonspatial_run[['non_sp_models']]
  cores = nonspatial_run[['cores']]
  rm(nonspatial_run) # Clean env

  result_de = stdiff_assemble_nonspatial_results(non_sp_models, meta_dict, pairwise, test_type, pval_adj)

  ######## ######## ######## BEGIN SPATIAL TESTS ######## ######## ########     ##### MADE BIG CHANGES TO NON-SPATIAL... NEED TO CHECK SPATIAL TESTS
  if(test_type == 'mm' & sp_topgenes > 0){
    start_t = Sys.time()
    if(verbose){
      cat('\tRunning spatial tests...\n')
    }
    spatial_run = stdiff_run_spatial_tests(x, result_de, non_sp_models, meta_dict, annot, pairwise, verbose, sp_topgenes, pval_thr, cores)
    result_de = spatial_run[['result_de']]
    sp_models = spatial_run[['sp_models']]
    rm(spatial_run) # Clean env
    end_t = difftime(Sys.time(), start_t, units='min')
    if(verbose){
      cat(paste0('\n\tCompleted spatial mixed models (', round(end_t, 2), ' min).\n'))
    }
    rm(non_sp_models) # Clean environment

    result_de = stdiff_assemble_spatial_results(result_de, sp_models, meta_dict, pairwise, pval_adj)
    rm(sp_models) # Clean environment
  }   ##### CLOSES IF TO ENTER SPATIAL TESTS CODE

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STdiff completed in ', round(end_t, 2), ' min.\n'))
  }
  return(result_de)
}
