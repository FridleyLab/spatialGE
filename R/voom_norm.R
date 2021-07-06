##
#' @title voom_norm: voom-transformation of ST arrays
#' @description Applies limma-voom transformation to spatial transcriptomics arrays data.
#' @details
#' This function takes an STList and returns roughly Gaussian count data matrices
#' in two steps. In the first step, (edgeR) normalization factors are used to scale
#' the counts from each library/spot. In the second step, limma-voom transformation
#' is applied. The resulting normalized count matrix is stored in the voom_counts
#' slot of the STList. The function also calculates gene-wise standard deviation from
#' the transformed counts and stores them in the gene_stdev slot.
#'
#' The function works on parallel using "forking" (not in Windows OS).
#'
#' @param x, an STList with raw count matrices.
#' @return x, an STList with normalized counts.
#' @export
#
#
voom_norm = function(x=NULL){

  require("magrittr")
  require("parallel")

  # Define number of available cores to use.
  cores = 1
  if(.Platform$OS.type == 'unix'){
    avail_cores = parallel::detectCores()
    if(avail_cores > (length(x@counts) + 1)){
      cores = (length(x@counts) + 1)
    } else if( (avail_cores <= (length(x@counts) + 1)) && avail_cores > 1){
      cores = avail_cores - 1
    }
  }

  # Test if an STList has been input.
  if(is.null(x) | !is(x, "STList")){
    stop("The input must be a STList.")
  }

  # Perform voom transformation on parallel if possible.
  x@voom_counts = mclapply(seq_along(x@counts), function(i){

    # Progress will probably show if cores=1
    cat(paste0("Normalizing spatial array #", i, "...\n"))

    # Calculate edgeR normalization factors.
    norm_factors = edgeR::calcNormFactors(x@counts[[i]][-1], method="TMM", lib.size=NULL)

    # Create new object to store the size-corrected counts.
    df = c()

    # Divide each count value by their respective column (sample) normalization factor.
    for (raw_col in names(x@counts[[i]][-1])) {
      df = dplyr::bind_cols(df, tibble::as_tibble(x@counts[[i]][raw_col]/norm_factors[raw_col]))
    }

    # Apply voom transformation to count data.
    df_voom = limma::voom(df, design = NULL, lib.size = colSums(df), normalize.method = "none", plot = F)

    # Put back gene names to matrix and store in object.
    df_voom = tibble::as_tibble(df_voom$E) %>%
      tibble::add_column(x@counts[[i]][, 1], .before = 1)

    return(df_voom)

  }, mc.cores=cores)

  # Calculate gene-wise standard deviations in parallel if possible.
  x@gene_stdev = mclapply(seq_along(x@counts), function(i){

    gene_stdevs_df = apply(x@voom_counts[[i]][, -1], 1, sd, na.rm = T)
    gene_stdevs_df = tibble::as_tibble_col(gene_stdevs_df, column_name='gene_stdevs') %>%
      tibble::add_column(x@counts[[i]][, 1], .before = 1)

    return(gene_stdevs_df)

  }, mc.cores=cores)

  return(x)

}

