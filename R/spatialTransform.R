##
#' @title spatialTransform: Data transformation of spatial transcriptomics data
#' @description Applies data transformation methods to spatia transcriptomics experiments of an STList
#' @details
#' This function takes an STList with raw counts and applies either logarithmic or
#' voom transformation after library size normalization. In the case of voom, library
#' size normalization applied is the trimmed mean of median values (TMM). If log,
#' counts are normalized using counts per spot. Logarithmic transformation also uses
#' a scaling factor (10^4) by default. The log, log(x+1) produces similar results to
#' Seurat's LogNormalize. The function works on parallel using "forking" (not in Windows OS).
#'
#' @param x an STList with raw count matrices.
#' @param method one of 'log' (default) or 'voom'.
#' @param scale_f, the scale factor used in logarithmic transformation.
#' @return x an updated STList with transformed counts.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # melanoma <- spatialTransform(melanoma)
#' # melanoma <- spatialTransform(melanoma, method='voom')
#'
#' @export
#
#
spatialTransform = function(x=NULL, method='log', scale_f=10000){

  require("magrittr")

  # Detect transformation method
  if(method == 'voom'){
    x@tr_counts = voom_norm(x)
    x@misc[['transform']] = 'voom'
  } else if(method == 'log'){
    x@tr_counts = log_transf(x, scale_f=scale_f)
    x@misc[['transform']] = 'log'
  }

  # Calculate gene-wise standard deviations in parallel if possible.
  cores = 1
  x@gene_var = parallel::mclapply(seq_along(x@tr_counts), function(i){
    gene_stdevs_df = apply(x@tr_counts[[i]][, -1], 1, sd, na.rm = T)
    #gene_stdevs_df = tibble::as_tibble_col(gene_stdevs_df, column_name='gene_stdevs') %>%
    gene_stdevs_df = tibble::tibble(gene_stdevs=gene_stdevs_df) %>%
      tibble::add_column(gene=x@tr_counts[[i]][['gene']], .before = 1)
    return(gene_stdevs_df)
  }, mc.cores=cores, mc.preschedule=F)
  # Copy names of the list
  names(x@gene_var) = names(x@counts)

  return(x)
}

# Helpers ----------------------------------------------------------------------

##
# @title log_transf: Log-transformation of spatial arrays
# @description Applies the natural logarithm [log(x+1)] to counts after library
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
#
log_transf = function(x=NULL, scale_f=scale_f){
  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STList")){
    stop("The input must be a STList.")
  }

  # Detect usable cores
  cores = count_cores(length(x@counts))

  # Perform voom transformation on parallel if possible.
  log_counts = parallel::mclapply(seq_along(x@counts), function(i){

    # Progress will probably show if cores=1
    cat(paste0("Log-transforming spatial sample: ", names(x@counts[i]), "...\n"))

    # Expand sparse matrix
    df = expandSparse(x@counts[[i]])

    # Calculate (spot) library sizes. Then, add 1 to each library size.
    libsizes = colSums(df)
    # Check that there are not zero-count spots
    if(any(libsizes == 0)){
      stop('Please, remove spots containing zero reads...')
    }

    # Remove gene names from data
    #df = x@counts[[i]][, -1]
    # Divide each count value by their respective column (spot) normalization factor.
    df = sweep(df, 2, libsizes, '/')
    # Then multiply by scaling factor
    df = df * scale_f
    # Apply log transformation to count data.
    df = log1p(df)

    # Put back gene names to matrix and store in object.
    df = tibble::as_tibble(df) %>%
      tibble::add_column(gene=rownames(x@counts[[i]]), .before=1)

    return(df)
  }, mc.cores=cores, mc.preschedule=F)

  # Copy list names from raw counts to transformed counts
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
#
voom_norm = function(x=NULL){
  # Define number of available cores to use.
  cores = count_cores(length(x@counts))

  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STList")){
    stop("The input must be a STList.")
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
