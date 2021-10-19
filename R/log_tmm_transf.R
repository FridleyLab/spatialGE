##
#' @title log_transf:
#' @description Applies
#' @details
#' This function takes
#'
#'
#'
#' The function works on parallel using "forking" (not in Windows OS).
#'
#' @param x an STList with raw count matrices.
#' @param cores integer, the number of cores to use in calculations
#' @return x an updated STList with transformed counts.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' #
#'
#' @export
#
#
log_tmm_transf = function(x=NULL, cores=1){

  require("magrittr")

  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STList")){
    stop("The input must be a STList.")
  }

  # Perform voom transformation on parallel if possible.
  x@log_counts = parallel::mclapply(seq_along(x@counts), function(i){

    # Progress will probably show if cores=1
    if(cores == 1){
      cat(paste0("Log-transforming spatial sample: ", names(x@counts[i]), "...\n"))
    }

    # Calculate edgeR normalization factors.
    norm_factors = edgeR::calcNormFactors(x@counts[[i]][, -1], method="TMM", lib.size=NULL)

    # Divide each count value by their respective column (spot) normalization factor.
    df = sweep(x@counts[[i]][, -1], 2, norm_factors, '/')

    # Apply log (1+x) transformation to count data.
    df = log1p(df)

    # Put back gene names to matrix and store in object.
    df = tibble::as_tibble(df) %>%
      tibble::add_column(x@counts[[i]][, 1], .before=1)

    return(df)

  }, mc.cores=cores, mc.preschedule=T)

  # Copy list names from raw counts to transformed counts
  names(x@log_counts) = names(x@counts)

  # Calculate gene-wise standard deviations in parallel if possible.
  x@log_stdev = parallel::mclapply(seq_along(x@counts), function(i){

    gene_stdevs_df = apply(x@log_counts[[i]][, -1], 1, sd)
    gene_stdevs_df = tibble::as_tibble_col(gene_stdevs_df, column_name='log_stdevs') %>%
      tibble::add_column(x@counts[[i]][, 1], .before=1)

    return(gene_stdevs_df)

  }, mc.cores=cores)

  # Copy list names from raw counts to log-transformed standard deviation list
  names(x@log_stdev) = names(x@counts)

  return(x)

}

