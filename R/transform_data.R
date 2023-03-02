##
#' @title transform_data: Data transformation of spatial transcriptomics data
#' @description Applies data transformation methods to spatial transcriptomics
#' samples within an STlist
#' @details
#' This function takes an STlist with raw counts and applies logarithmic transformation
#' after library size normalization. This logarithmic transformation uses
#' a scaling factor (10^4) by default.  The function works on parallel using
#' "forking" (not available in Windows OS).
#'
#' @param x an STList with raw count matrices.
#' @param scale_f, the scale factor used in logarithmic transformation.
#' @return x an updated STList with transformed counts.
#'
#' @examples
#' # In this example, melanoma is an STlist.
#' # melanoma <- transform_data(melanoma)
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
#
transform_data = function(x=NULL, scale_f=10000){

  method = 'log'

  # Detect transformation method
  if(method == 'voom'){
    x@tr_counts = voom_norm(x)
    x@misc[['transform']] = 'voom'
  } else if(method == 'log'){
    # log-normalize counts and obtain spots/cells with zero counts
    log_results = log_transf(x, scale_f=scale_f)
    for(i in names(x@counts)){
      # Remove zero counts spots/cells if any
      if(!is.null(log_results[[i]][['zero_size']])){
        zero_size = log_results[[i]][['zero_size']]
        x@counts[[i]] = x@counts[[i]][, -zero_size]
        libnames = colnames(x@counts[[i]])
        x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
          dplyr::filter(libname %in% libnames)
        rm(zero_size, libnames) # Clean environment
      }
      x@tr_counts[[i]] = makeSparse(log_results[[i]][['counts']])
    }
    x@misc[['transform']] = 'log'
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
log_transf = function(x=NULL, scale_f=NULL){
  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STlist")){
    stop("The input must be a STlist.")
  }

  # Detect usable cores
  cores = count_cores(length(x@counts))

  # Perform log-transformation on parallel if possible.
  log_counts = parallel::mclapply(seq_along(x@counts), function(i){
    # Show progress
    system(sprintf('echo "%s"', crayon::yellow(paste0("Log-transforming sample ", names(x@counts)[i], "...."))))

    df_tmp = as.matrix(x@counts[[i]])
    # Calculate (spot/cell) library sizes.
    libsizes = colSums(df_tmp)
    # Check that there are not zero-count spot/cells
    zero_size = NULL
    if(any(libsizes == 0)){
      zero_size = as.vector(which(libsizes == 0))
      df_tmp = df_tmp[, -zero_size]
      system(sprintf('echo "%s"', crayon::red(paste0(length(zero_size), " spots/cells with zero counts will be removed from sample ", i, " and from the entire STlist."))))
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
# @title voom_norm: voom transformation of ST arrays
# @description Applies limma-voom transformation to spatial transcriptomics arrays data.
# @details
# This function takes an STList and returns roughly Gaussian count data matrices
# in two steps. In the first step, (edgeR) normalization factors are used to scale
# the counts from each library/spot. In the second step, limma-voom transformation
# is applied. The resulting normalized count matrix is stored in the voom_counts
# slot of the STList. The function also calculates gene-wise standard deviation from
# the transformed counts and stores them in the gene_stdev slot.
#
# The function works on parallel using "forking" (not in Windows OS).
#
# @param x an STList with raw count matrices.
# @return x an updated STList with transformed counts.
#
# @examples
# # In this example, melanoma is an STList.
# # melanoma <- voom_norm(melanoma)
#
#' @importFrom methods as is new
#
#
voom_norm = function(x=NULL){
  # Define number of available cores to use.
  cores = count_cores(length(x@counts))

  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STlist")){
    stop("The input must be a STlist.")
  }

  # Perform voom transformation on parallel if possible.
  voom_counts = parallel::mclapply(seq_along(x@counts), function(i){

    # Progress will probably show if cores=1
    cat(paste0("Normalizing spatial sample: ", names(x@counts[i]), "...\n"))

    # Expand sparse matrix
    df = expandSparse(x@counts[[i]])

    # Calculate edgeR normalization factors.
    norm_factors = edgeR::calcNormFactors(df, method="TMM", lib.size=NULL)

    # Divide each count value by their respective column (spot) normalization factor.
    df = sweep(df, 2, norm_factors, '/')

    # Apply voom transformation to count data.
    df_voom = limma::voom(df, design=NULL, lib.size=colSums(df), normalize.method="none", plot=F)

    # Put back gene names to matrix and store in object.
    df_voom = tibble::as_tibble(df_voom$E) %>%
      tibble::add_column(gene=rownames(x@counts[[i]]), .before=1)

    return(df_voom)

  }, mc.cores=cores, mc.preschedule=F)

  # Copy list names from raw counts to transformed counts
  names(voom_counts) = names(x@counts)

  return(voom_counts)
}
