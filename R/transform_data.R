##
#' @title transform_data: Transformation of spatial transcriptomics data
#' @description Applies data transformation methods to spatial transcriptomics
#' samples within an STlist
#' @details
#' This function takes an STlist with raw counts and performs data transformation.
#' The user has the option to select between log transformation after library size
#' normalization (`method='log'`), or SCTransform (`method='sct'`). In the case of
#' logarithmic transformation, a scaling factor (10^4 by default) is applied. The
#' function uses parallelization using "forking" (not available in Windows OS).
#' Note that the method `sct` returns a matrix with less genes as filtering is
#' done for low expression genes.
#'
#' @param x an STlist with raw count matrices.
#' @param method one of `log` or `sct`. If `log`, log-normalization is performed.
#' If `sct`, then the SCTransform method is applied by calling `sctransform::vst`
#' @param scale_f the scale factor used in logarithmic transformation
#' @param sct_n_regr_genes the number of genes to be used in the regression model
#' during SCTransform. The function `sctransform::vst` makes a random gene selection
#' based on this number
#' @param sct_min_cells The minimum number of spots/cells to be used in the regression
#' model fit by `sctransform::vst`
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems.
#' @return x an updated STlist with transformed counts.
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
#
transform_data = function(x=NULL, method='log', scale_f=10000, sct_n_regr_genes=3000, sct_min_cells=5, cores=NULL){

  # Test if an STList has been input.
  assert_stlist(x)

  # Reject an unrecognized method up front. Without this the if/else below falls
  # through without assigning tr_results, and the failure surfaces 10 lines later
  # as 'object tr_results not found' -- an internal name the caller cannot act on.
  if(!is.character(method) || length(method) != 1 || !(method %in% c('log', 'sct'))){
    raise_err(err_code='error0043')
  }

  # Detect transformation method
  # if(method == 'voom'){
  #   stop('voom normalization has been deprecated. Please use "log" or "sct".')
  #   #x@tr_counts = voom_norm(x)
  #   #x@misc[['transform']] = 'voom'
  # } else
    if(method == 'log'){
    # log-normalize counts and obtain spots/cells with zero counts
    tr_results = log_transf(x, scale_f=scale_f, cores=cores)
    x@misc[['transform']] = 'log'
  } else if(method == 'sct'){
    # Apply SCTransform counts and obtain spots/cells with zero counts
    tr_results = sct_transf(x, sct_n_regr_genes=sct_n_regr_genes, sct_min_cells=sct_min_cells, cores=cores)
    x@misc[['transform']] = 'sct'
  }

  # Put transformed counts in STlist
  for(i in names(x@counts)){
    if(class(tr_results[[i]][['counts']])[1] != "dgCMatrix"){
      x@tr_counts[[i]] = makeSparse(tr_results[[i]][['counts']])
    } else{
      x@tr_counts[[i]] = tr_results[[i]][['counts']]
    }
    # Remove zero counts spots/cells if any from raw counts
    if(!is.null(tr_results[[i]][['zero_size']])){
      zero_size = tr_results[[i]][['zero_size']]
      x@counts[[i]] = x@counts[[i]][, -zero_size]
      libnames_nonzero = sample_cells(x, i, 'counts')
      x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
        dplyr::filter(libname %in% libnames_nonzero)
      rm(zero_size, libnames_nonzero) # Clean environment
    }
  }

  # Calculate gene-wise standard deviations in parallel if possible.
  cores = resolve_cores(length(x@tr_counts), cores, coerce='as_integer')
  # X is the list of transformed matrices itself, so the closure captures
  # nothing at all.
  x@gene_meta = sge_lapply(x@tr_counts, function(tr_mat){
    gene_meta_one_sample(tr_mat)
  }, cores=cores, preschedule=F)
  # Copy names of the list
  names(x@gene_meta) = names(x@counts)

  return(x)
}

# Helpers ----------------------------------------------------------------------

##
# @title log_transf: Log-transformation of spatial arrays
# @description Applies the natural logarithm log(x+1) to counts after library
# size normalization and multiplication by scaling factor (by default 10^4).
# The function produces similar results to Seurat's LogNormalize.
#
# @details
# The function works on parallel using "forking" (not in Windows OS).
#
# @param x an STList with raw count matrices.
# @param scale_f integer, the scaling factor to use.
# @return x an updated STList with transformed counts.
#
#' @importFrom methods as is new
#
#
log_transf = function(x=NULL, scale_f=NULL, cores=NULL){
  # Test if an STList has been input.
  assert_stlist(x)

  # Detect usable cores
  cores = resolve_cores(length(x@counts), cores, coerce='as_integer')
  if(is.na(cores)){
    stop('Could not recognize number of cores requested')
  }

  # Perform log-transformation on parallel if possible.
  #
  # The payload is built OUTSIDE the call and passed as X, so the closure
  # captures only scale_f rather than the whole STlist. See the note on
  # sge_lapply payloads in R/utils_parallel.R for what this does and does not
  # buy under forking.
  payloads = Map(function(counts_mat, sample_name){
    list(counts=counts_mat, name=sample_name)
  }, x@counts, names(x@counts))

  log_counts = sge_lapply(payloads, function(p){
    log_transf_one_sample(p$counts, p$name, scale_f)
  }, cores=cores, preschedule=F)

  # Copy list names from raw counts to transformed counts
  names(log_counts) = names(x@counts)
  names(log_counts) = names(x@counts)

  return(log_counts)
}


##
# @title sct_transf: Applies SCTransform on spatial samples
# @description Applies Seurat's SCTransform method.
#
# @details
# The function works on parallel using "forking" (not in Windows OS).
#
# @param x an STList with raw count matrices.
# @return x an updated STList with transformed counts.
#
#' @importFrom methods as is new
#
#
sct_transf = function(x=NULL, sct_n_regr_genes=3000, sct_min_cells=5, cores=NULL){
  # Test if an STlist has been input.
  assert_stlist(x)

  # Detect usable cores
  cores = resolve_cores(length(x@counts), cores, coerce='as_integer')
  if(is.na(cores)){
    stop('Could not recognize number of cores requested')
  }

  # Payload built outside the call so the closure captures only the two scalar
  # parameters. See R/utils_parallel.R on what this buys.
  payloads = Map(function(counts_mat, sample_name){
    list(counts=counts_mat, name=sample_name)
  }, x@counts, names(x@counts))

  sct_counts = sge_lapply(payloads, function(p){
    sct_transf_one_sample(p$counts, p$name, sct_n_regr_genes, sct_min_cells)
  }, cores=cores, preschedule=F)

  # Copy list names from raw counts to transformed counts
  names(sct_counts) = names(x@counts)
  names(sct_counts) = names(x@counts)

  return(sct_counts)
}


##
# @title log_transf_one_sample: Log-transform a single sample's raw counts
# @description The per-sample unit of work dispatched by `log_transf()`'s
# `sge_lapply` call, extracted so it is directly callable/benchmarkable for
# one sample without going through the whole STlist or the parallel dispatch.
#
# @param counts_mat a raw count matrix (spots/cells in columns).
# @param sample_name the sample's name, used only in progress messages.
# @param scale_f the scaling factor to use.
# @return a list with `counts` (the transformed matrix) and `zero_size`
# (indices of zero-count spots/cells removed, or NULL).
#
log_transf_one_sample = function(counts_mat, sample_name, scale_f){
  # Show progress
  system(sprintf('echo "%s"', paste0("Log-transforming sample ", sample_name, "....")))

  # Runs entirely on the sparse matrix. This used to open with
  # as.matrix(counts_mat) and stay dense through sweep(), the scale_f multiply
  # and log1p -- three full dense genes x cells allocations -- after which
  # transform_data() called makeSparse() to convert the result back. So the path
  # was sparse -> dense -> sparse, per sample and, under mclapply, per fork.
  #
  # It is safe because log1p(0) == 0 exactly, so the transform maps structural
  # zeros to structural zeros and sparsity is preserved rather than destroyed and
  # rebuilt.
  df_tmp = counts_mat
  if(!inherits(df_tmp, 'dgCMatrix')){
    df_tmp = methods::as(methods::as(methods::as(df_tmp, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
  }

  # Calculate (spot/cell) library sizes.
  libsizes = Matrix::colSums(df_tmp)
  # Check that there are not zero-count spot/cells
  zero_size = NULL
  if(any(libsizes == 0)){
    zero_size = as.vector(which(libsizes == 0))
    df_tmp = df_tmp[, -zero_size, drop=FALSE]
    system(sprintf('echo "%s"', paste0(length(zero_size), " spots/cells with zero counts will be removed from sample ", sample_name, " and from the entire STlist.")))
    libsizes = libsizes[libsizes != 0]
  }

  # Divide by the column's library size, scale, then log-transform -- in that
  # order, matching the dense sweep()/multiply/log1p sequence exactly, so the
  # stored values are bit-identical rather than merely close.
  #
  # @p holds the column boundaries of a CSC matrix, so diff(@p) is the number of
  # stored entries per column and rep.int() expands that into the column index of
  # each entry in @x. One pass over the non-zeros, no dense intermediate.
  # unname() is load-bearing: Matrix::colSums() returns a vector named by
  # colnames, so libsizes[col_of_entry] would carry those names into @x, and a
  # dgCMatrix's @x must be an unnamed numeric vector. The values are unaffected,
  # which is exactly why this needed the golden gate to catch -- comparing
  # as.matrix() of both sides reports them as identical.
  if(length(df_tmp@x) > 0){
    col_of_entry = rep.int(seq_len(ncol(df_tmp)), diff(df_tmp@p))
    df_tmp@x = log1p(df_tmp@x / unname(libsizes)[col_of_entry] * scale_f)
  }

  # drop0() so the result matches what the previous dense -> makeSparse() route
  # produced: that path dropped every zero on conversion, whereas operating in
  # place would carry an explicitly stored zero through as an explicit
  # log1p(0) == 0. A no-op unless @counts already held explicit zeros.
  df_tmp = Matrix::drop0(df_tmp)

  result_list = list(counts=df_tmp, zero_size=zero_size)

  return(result_list)
}


##
# @title sct_transf_one_sample: Apply SCTransform to a single sample's raw counts
# @description The per-sample unit of work dispatched by `sct_transf()`'s
# `sge_lapply` call, extracted so it is directly callable/benchmarkable for
# one sample without going through the whole STlist or the parallel dispatch.
#
# @param counts_mat a raw count matrix (spots/cells in columns).
# @param sample_name the sample's name, used only in progress messages.
# @param sct_n_regr_genes the number of genes used in `sctransform::vst`'s
# regression model.
# @param sct_min_cells the minimum number of cells used in `sctransform::vst`'s
# regression model.
# @return a list with `counts` (the transformed matrix) and `zero_size`
# (indices of zero-count spots/cells removed, or NULL).
#
sct_transf_one_sample = function(counts_mat, sample_name, sct_n_regr_genes, sct_min_cells){
  # Show progress
  system(sprintf('echo "%s"', paste0("Applying SCTransform to sample ", sample_name, "....")))

  df_tmp = as.matrix(counts_mat)
  # Calculate (spot/cell) library sizes.
  libsizes = colSums(df_tmp)
  # Check that there are not zero-count spot/cells
  zero_size = NULL
  if(any(libsizes == 0)){
    zero_size = as.vector(which(libsizes == 0))
    df_tmp = df_tmp[, -zero_size]
    system(sprintf('echo "%s"', paste0(length(zero_size), " spots/cells with zero counts will be removed from sample ", sample_name, " and from the entire STlist.")))
    libsizes = libsizes[libsizes != 0]
  }

  # Apply SCTransform
  df_tmp = sctransform::vst(df_tmp, n_genes=sct_n_regr_genes, min_cells=sct_min_cells,
                            return_corrected_umi=T, verbosity=0)
  # log transform the corrected UMI counts from SCTransform
  df_tmp = log1p(df_tmp[['umi_corrected']])

  result_list = list(counts=df_tmp, zero_size=zero_size)

  return(result_list)
}


##
# @title gene_meta_one_sample: Compute per-gene mean/stdev for a single sample
# @description The per-sample unit of work dispatched by `transform_data()`'s
# gene-wise standard deviation `sge_lapply` call, extracted so it is directly
# callable/benchmarkable for one sample without going through the whole STlist.
#
# @param tr_counts_mat a transformed count matrix for one sample (genes in rows).
# @return a tibble with columns `gene`, `gene_mean`, `gene_stdevs`.
#
gene_meta_one_sample = function(tr_counts_mat){
  gene_stdevs_df = tibble::tibble(gene_mean=as.vector(Matrix::rowMeans(tr_counts_mat))) %>%
    tibble::add_column(gene_stdevs=as.vector(sparse_row_sds(tr_counts_mat))) %>%
    tibble::add_column(gene=rownames(tr_counts_mat), .before=1)
  return(gene_stdevs_df)
}

