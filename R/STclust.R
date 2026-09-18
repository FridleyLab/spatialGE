##
#' @title STclust: Detect clusters of spots/cells
#' @description Perform unsupervised spatially-informed clustering on the spots/cells of a
#' ST sample
#' @details
#' The function takes an STlist and calculates euclidean distances between cells or spots
#' based on the x,y spatial locations, and the expression of the top variable genes
#' (`Seurat::FindVariableFeatures`). The resulting distances are weighted by
#' applying 1-`ws` to the gene expression distances and `ws` to the spatial distances.
#' Hierarchical clustering is performed on the sum of the weighted distance matrices.
#' The `STclust` method allows for identification of tissue niches/domains that are
#' spatially cohesive.
#'
#' @param x an STlist with normalized expression data
#' @param samples a vector with strings or a vector with integers indicating the samples
#' to run STclust
#' @param ws a double (0-1) indicating the weight to be applied to spatial distances.
#' Defaults to 0.025
#' @param dist_metric the distance metric to be used. Defaults to 'euclidean'. Other
#' options are the same as in `wordspace::dist.matrix`
#' @param linkage the linkage method applied to hierarchical clustering. Passed to
#' `hclust` and defaults to 'ward.D'
#' @param ks the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied
#' @param topgenes the number of genes with highest spot-to-spot expression variation. The
#' variance is calculated via `Seurat::FindVariableFeatures`.
#' @param deepSplit a logical or integer (1-4), to be passed to `cutreeDynamic` and
#' control cluster resolution
#' @param cores an integer indicating the number of cores to use in parallelization (Unix only)
#' @param verbose either logical or an integer (0, 1, or 2) to increase verbosity
#' @return an STlist with cluster assignments
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- transform_data(melanoma)
#'   melanoma <- STclust(melanoma, ws=c(0, 0.025))
#'   STplot(melanoma, ws=0.025, samples='ST_mel1_rep2', ptsize=1)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#' @importFrom stats as.dist complete.cases cutree dist hclust prcomp sd na.omit
#
STclust = function(x=NULL, samples=NULL, ws=0.025, dist_metric='euclidean', linkage='ward.D2', ks='dtc', topgenes=2000, deepSplit=FALSE, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Record time
  zero_t = Sys.time()
  if(verbose){
    cat(paste0('STclust started...\n'))
  }

  # Clustering method to use. Set because other methods will be supported in future versions
  clmethod = 'hclust'

  # Force ws and topgenes as numeric in case entered as character
  ws = as.double(ws)
  topgenes = as.integer(topgenes)

  # Do not allow weights higher than 1
  if(any(ws < 0) | any(ws > 1)){
    stop('Please select a spatial weight between 0 and 1.')
  }

  # Check to ensure number of ks is acceptable
  if(is.numeric(ks)){
    ks = as.integer(ks)
    if(length(ks) == 1 & ks[1] < 2){
      raise_err(err_code='error0016')
    } else if(any(ks < 2)){
      warning('Refusing to generate < 2 clusters. Skipping any k < 2.')
      ks = ks[ks >= 2]
    }
  }

  # Test if an STList has been input.
  x = assert_stlist(x)

  # Check data has been normalized. Otherwise raise error
  assert_transformed(x, err_code='error0007')

  # Define samples using names (convert indexes to names if necessary)
  samples = resolve_samples(x, samples, slot='spatial_meta', numeric_handling='na_omit',
                            char_transform='exact', post_check='any_in',
                            message='None of the requested samples are present in the STlist.')

  # Define number of cores for parallelization of tests
  cores = resolve_cores(length(samples), cores, coerce='ceiling')

  # Identify variable genes (Seurat's VST)
  x = calculate_vst(x=x, samples=samples, cores=cores)

  # Subset variable genes
  trcounts_df = sge_lapply(samples, function(i){
    stclust_topgenes_one_sample(sample_gene_meta(x, i), sample_tr_counts(x, i), topgenes)
  }, cores=cores)
  names(trcounts_df) = samples

  # Parallelize clustering
  res_ls = sge_lapply(samples, function(i){
    stclust_one_sample(trcounts_df[[i]], sample_spatial_meta(x, i), dist_metric, ws, clmethod,
                       ks, deepSplit, linkage, verbose, i)
  }, cores=cores)
  names(res_ls) = samples

  rm(trcounts_df) # Clean env

  # Add results to STlist
  if(verbose){
    cat('Updating STlist with results...\n')
  }
  lapply(samples, function(i){
    for(w in 1:length(ws)){
      if(any(colnames(x@spatial_meta[[i]])[-c(1:5)] %in% colnames(res_ls[[i]][[w]]))){
        col_names = intersect(colnames(x@spatial_meta[[i]])[-1], colnames(res_ls[[i]][[w]]))
        x@spatial_meta[[i]] <<- x@spatial_meta[[i]] %>% dplyr::select(-!!col_names)
      }
      x@spatial_meta[[i]] <<- x@spatial_meta[[i]] %>% dplyr::full_join(., res_ls[[i]][[w]], by='libname')
    }
  })

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose){
    cat(paste0('STclust completed in ', round(end_t, 2), ' min.\n'))
  }

  return(x)
}


# Helpers ----------------------------------------------------------------------

##
# @title stclust_topgenes_one_sample: Subset a sample's counts to its top
# variable genes
# @description The per-sample unit of work dispatched by `STclust()`'s
# "Subset variable genes" `sge_lapply` call, extracted so it is directly
# callable/benchmarkable for one sample without going through the whole STlist.
#
# @param gene_meta_df the sample's gene metadata tibble (must contain
# `gene` and `vst.variance.standardized` columns).
# @param tr_counts_mat the sample's transformed count matrix.
# @param topgenes the number of top variable genes to keep.
# @return the transformed count matrix subset to the top variable genes.
#
stclust_topgenes_one_sample = function(gene_meta_df, tr_counts_mat, topgenes){
  # To prevent NOTES in R CMD check
  . = NULL

  topgenenames_tmp = gene_meta_df %>%
    dplyr::arrange(dplyr::desc(vst.variance.standardized)) %>%
    dplyr::slice_head(n=topgenes) %>%
    dplyr::select(gene) %>%
    unlist() %>%
    as.vector()

  # Get transformed counts
  trcounts_df_tmp = tr_counts_mat[rownames(tr_counts_mat) %in% topgenenames_tmp, ]

  return(trcounts_df_tmp)
}


##
# @title stclust_one_sample: Run spatially-weighted hierarchical clustering
# for one sample
# @description The per-sample unit of work dispatched by `STclust()`'s
# "Parallelize clustering" `sge_lapply` call, extracted so it is directly
# callable/benchmarkable for one sample without going through the whole STlist.
#
# @param trcounts_mat the sample's transformed counts, subset to top variable genes.
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param dist_metric the distance metric, passed to `calculate_dist_matrices`.
# @param ws vector of spatial weights.
# @param clmethod the clustering method (only 'hclust' is currently supported).
# @param ks the range of k values, or 'dtc' to apply `cutreeDynamic`.
# @param deepSplit passed to `get_hier_clusters_dtc`.
# @param linkage the linkage method, passed to `hclust`.
# @param verbose either logical or an integer (0, 1, or 2); values > 1 print
# a per-sample progress message.
# @param sample_name the sample's name, used only in the progress message.
# @return a list of data frames with cluster assignments, one per spatial weight.
#
stclust_one_sample = function(trcounts_mat, spatial_meta_df, dist_metric, ws, clmethod,
                              ks, deepSplit, linkage, verbose, sample_name){
  # Calculate scaled expression and spatial distance matrices
  scaled_dists = calculate_dist_matrices(expr_dat=trcounts_mat, coord_dat=spatial_meta_df, dist_metric=dist_metric)

  # Calculate weighted distance matrices
  weighted_dists = calculate_weighted_dist(scaled_dists=scaled_dists, ws=ws)

  rm(scaled_dists) # Clean env

  # Identify clustering method... ONLY HIERARCHICAL CLUSTERING IMPLEMENTED SO FAR
  if(clmethod == 'hclust'){

    # Apply dtc or split to k
    if(as.character(ks[1]) == 'dtc'){
      # Hierarchical clustering using DynamicTreeCLusters
      hierclusters_ls = get_hier_clusters_dtc(weighted_dists=weighted_dists, ws=ws, deepSplit=deepSplit, linkage=linkage)
    } else if(is.numeric(ks)){
      # Hierarchical clustering using range of Ks
      hierclusters_ls = get_hier_clusters_ks(weighted_dists=weighted_dists, ws=ws, ks=ks, linkage=linkage)
    } else{
      stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
    }

  } else{
    stop('Currently, only spatially-informed hierarchical clustering is supported.')
  }

  if(verbose > 1L){
    system(sprintf('echo "%s"', paste0("\tClustering completed for ", sample_name, "...")))
  }

  return(hierclusters_ls)
}


##
# calculate_dist_matrices: Calculate distance matrices and scale to 1
# @param expr_dat an sparse matrix with gene expression, genes in rows, spots/cells in columns
# @param coord_dat a data frame with three colums: 'libname', 'xpos', 'ypos'
# @param dist_metric a string indicating the type of distance to calculate
# @return a list with two matrices (scaled expression and scaled coordinates distance matrices)
#
calculate_dist_matrices = function(expr_dat=NULL, coord_dat=NULL, dist_metric=NULL){
  a = Matrix::t(expr_dat)
  b = coord_dat[, c('libname', 'xpos', 'ypos')] %>% tibble::column_to_rownames(var='libname')
  b = as.matrix(b[match(rownames(b), rownames(a)), ])

  # Get distance matrices
  da = wordspace::dist.matrix(a, method=dist_metric)
  db = dist(b, upper=T, diag=T, method=dist_metric)
  dam = as.matrix(da)
  dbm = as.matrix(db)

  rm(a, b, da, db) # Clean env

  # Scale matrices
  dam = dam/max(dam)
  dbm = dbm/max(dbm)

  return(list(scale_exp=dam, scale_coord=dbm))
}

##
# get_weighted_dist: Calculated weighted matrixes using spatial weights
# @param scaled_dists a list with two matrices, a scaled gene expression matrix and a scaled coordinate matrix
# @param ws vector with requested spatial weights
# @return a list of weighted matrices
#
calculate_weighted_dist = function(scaled_dists=NULL, ws=NULL){
  weight_mtx_ls = lapply(1:length(ws), function(w){
    weight_d = ws[w]
    weight_g = 1-weight_d

    # Create vector of weights for Reduce
    weight_ls = c(weight_g, weight_d)
    dmxs = list(scaled_dists[[1]], scaled_dists[[2]])

    # Apply weight element-wise
    weight_mtx = Reduce('+', Map('*', dmxs, weight_ls))

    return(weight_mtx)
  })

  return(weight_mtx_ls)
}

##
# get_clusters_dtc: Performs hierarchical clustering followed by DynamicTreeCuts
# @param weighted_dists a list of distance matrices (NOT dist objects) for each spatial weight
# @param ws a vector with spatial weights
# @param deepSplit a logical or number indicating whether to use deepSplit in DTC
# @param linkage a string with the linkage method to use
# @return a list of data frames with spot/cells cluster assignments for each weight
#
get_hier_clusters_dtc = function(weighted_dists=NULL, ws=NULL, deepSplit=NULL, linkage=NULL){
  grp_df_ls = lapply(1:length(ws), function(w){
    # Construct column name to be put in `spatial_meta` based on weight and deepSplit
    if(is.logical(deepSplit)){
      dspl = 'False'
      if(deepSplit){
        dspl = 'True'
      }
    } else{
      dspl = deepSplit
    }
    col_name = paste0('stclust_spw', ws[w], '_dspl', dspl)

    # Run hierarchical clustering
    hierclusters = hclust(as.dist(weighted_dists[[w]]), method=linkage)

    # Use DynamicTreeClusters
    grp_df = dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weighted_dists[[w]], deepSplit=deepSplit, verbose=F)

    # Create data frame with cluster assignments
    grp_df = tibble::tibble(libname=colnames(weighted_dists[[w]]), !!col_name:=as.factor(as.vector(grp_df)))
    # Convert zeroes to NAs
    grp_df[[col_name]][grp_df[[col_name]] == 0] = NA

    return(grp_df)
  })

  return(grp_df_ls)
}

##
# get_clusters_ks: Performs hierarchical clustering followed by cluster assignment (cuttree)
# @param weighted_dists a list of distance matrices (NOT dist objects) for each spatial weight
# @param ws a vector with spatial weights
# @param ks a vector with k values detect
# @param linkage a string with the linkage method to use
# @return a list of data frames with spot/cells cluster assignments for each weight
#
get_hier_clusters_ks = function(weighted_dists=NULL, ws=NULL, ks=NULL, linkage=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  grp_df_ls = lapply(1:length(ws), function(w){
    grp_df = tibble::tibble(libname=colnames(weighted_dists[[w]]))

    # Run hierarchical clustering once per weight: weighted_dists[[w]] and
    # linkage do not change across k, only the dendrogram cut does.
    hierclusters = hclust(as.dist(weighted_dists[[w]]), method=linkage)

    for(k in ks){
      # Construct column name to be put in `spatial_meta` based on weight and k
      col_name = paste0('stclust_spw', ws[w], '_k', k)

      # Cut the dendrogram
      grp_df_tmp = cutree(hierclusters, k=k)
      # Create data frame with cluster assignments
      grp_df_tmp = tibble::tibble(libname=colnames(weighted_dists[[w]]), !!col_name:=as.factor(as.vector(grp_df_tmp)))

      grp_df = grp_df %>% dplyr::left_join(., grp_df_tmp, by='libname')

      rm(grp_df_tmp, col_name) # Clean env
    }
    rm(hierclusters) # Clean env
    return(grp_df)
  })

  return(grp_df_ls)
}

