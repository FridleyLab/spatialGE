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
##'
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#' unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#' # Generate the file paths to be passed to the STlist function
#' count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='counts')
#' coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='mapping')
#' clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                         full.names=TRUE, pattern='clinical')
#' # Create STlist
#' library('spatialGE')
#' melanoma <- STlist(rnacounts=count_files[c(1,2)],
#'                    spotcoords=coord_files[c(1,2)],
#'                    samples=clin_file) # Only first two samples
#' melanoma <- transform_data(melanoma, method='log')
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
transform_data = function(x=NULL, method='log', scale_f=10000, sct_n_regr_genes=3000, sct_min_cells=5, cores=NULL){

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
      libnames_nonzero = colnames(x@counts[[i]])
      x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
        dplyr::filter(libname %in% libnames_nonzero)
      rm(zero_size, libnames_nonzero) # Clean environment
    }
  }

  # Calculate gene-wise standard deviations in parallel if possible.
  cores = 1
  x@gene_meta = parallel::mclapply(seq_along(x@tr_counts), function(i){
    #gene_stdevs_df = apply(x@tr_counts[[i]], 1, sd, na.rm=T) # There shouldn't be any NAs in data transformed counts
    gene_stdevs_df = tibble::tibble(gene_mean=as.vector(Matrix::rowMeans(x@tr_counts[[i]]))) %>%
      tibble::add_column(gene_stdevs=as.vector(apply(x@tr_counts[[i]], 1, sd))) %>%
      tibble::add_column(gene=rownames(x@tr_counts[[i]]), .before=1)
    return(gene_stdevs_df)
  }, mc.cores=cores, mc.preschedule=F)
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
  if(is.null(x) | !is(x, "STlist")){
    stop("The input must be a STlist.")
  }

  # Detect usable cores
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(x@counts))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Perform log-transformation on parallel if possible.
  log_counts = parallel::mclapply(seq_along(x@counts), function(i){
    # Show progress
    system(sprintf('echo "%s"', paste0("Log-transforming sample ", names(x@counts)[i], "....")))

    df_tmp = as.matrix(x@counts[[i]])
    # Calculate (spot/cell) library sizes.
    libsizes = colSums(df_tmp)
    # Check that there are not zero-count spot/cells
    zero_size = NULL
    if(any(libsizes == 0)){
      zero_size = as.vector(which(libsizes == 0))
      df_tmp = df_tmp[, -zero_size]
      system(sprintf('echo "%s"', paste0(length(zero_size), " spots/cells with zero counts will be removed from sample ", names(x@counts)[i], " and from the entire STlist.")))
      libsizes = libsizes[libsizes != 0]
    }

    # Divide each count value by their respective column (spot) normalization factor.
    df_tmp = sweep(df_tmp, 2, libsizes, '/')
    # Then multiply by scaling factor
    df_tmp = df_tmp * scale_f
    # Apply log transformation to count data.
    df_tmp = log1p(df_tmp)

    # # Put back gene names to matrix and store in object.
    # df = tibble::as_tibble(df) %>%
    #   tibble::add_column(gene=rownames(x@counts[[i]]), .before=1)

    result_list = list(counts=df_tmp, zero_size=zero_size)

    return(result_list)
  }, mc.cores=cores, mc.preschedule=F)

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
  if(is.null(x) | !is(x, "STlist")){
    stop("The input must be a STlist.")
  }

  # Detect usable cores
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(x@counts))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Perform log-transformation on parallel if possible.
  sct_counts = parallel::mclapply(seq_along(x@counts), function(i){
    # Show progress
    system(sprintf('echo "%s"', paste0("Applying SCTransform to sample ", names(x@counts)[i], "....")))

    df_tmp = as.matrix(x@counts[[i]])
    # Calculate (spot/cell) library sizes.
    libsizes = colSums(df_tmp)
    # Check that there are not zero-count spot/cells
    zero_size = NULL
    if(any(libsizes == 0)){
      zero_size = as.vector(which(libsizes == 0))
      df_tmp = df_tmp[, -zero_size]
      system(sprintf('echo "%s"', paste0(length(zero_size), " spots/cells with zero counts will be removed from sample ", names(x@counts)[i], " and from the entire STlist.")))
      libsizes = libsizes[libsizes != 0]
    }

    # Apply SCTransform
    df_tmp = sctransform::vst(df_tmp, n_genes=sct_n_regr_genes, min_cells=sct_min_cells,
                              return_corrected_umi=T, verbosity=0)
    # log transform the corrected UMI counts from SCTransform
    df_tmp = log1p(df_tmp[['umi_corrected']])

    result_list = list(counts=df_tmp, zero_size=zero_size)

    return(result_list)
  }, mc.cores=cores, mc.preschedule=F)

  # Copy list names from raw counts to transformed counts
  names(sct_counts) = names(x@counts)
  names(sct_counts) = names(x@counts)

  return(sct_counts)
}

