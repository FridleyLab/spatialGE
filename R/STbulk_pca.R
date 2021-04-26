##
#' @title STBulk_pca
#' @description Perform and plot a PCA of simulated bulk RNA-Seq ST data.
#' @details
#' This function takes an STList and optionally the name of a clinical variable
#' and performs PCA on simulated bulk RNA-Seq from ST data. The counts are simulated
#' by summing all counts from an array for a given gene.
#'
#' @param x, an STList.
#' @param clinvar, a string indicating the name of the variable in the clinical
#' data frame.
#' @export
#
#
STbulk_pca <- function(x=NULL, clinvar=NULL) {

  require('magrittr')
  require('ggplot2')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  if(!is.null(clinvar)){
    clinvar_vals <- as.character(x@clinical[[clinvar]])
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='var_vals')
  }

  # Create data frame to store "bulk counts".
  bulkcounts_df <- x@counts[[1]][, 1]

  # Loop through count matrices in STList
  for(i in 1:length(x@counts)){

    bulk_expr <- rowSums(x@counts[[i]][, -1])
    bulk_expr <- tibble::as_tibble_col(bulk_expr,
                                       column_name=paste0("st_", i))
    bulk_expr <- bulk_expr %>%
      tibble::add_column(x@counts[[i]][, 1], .before=1)

    bulkcounts_df <- dplyr::inner_join(bulkcounts_df,
                                       bulk_expr, by=c('gene'))
    }

  norm_factors <- edgeR::calcNormFactors(bulkcounts_df[, -1], method='TMM',
                                         lib.size=NULL)

  # Create new object to store the size-corrected counts.
  bulkvoom_df <- c()

  # Divide each count value by their respective column (sample) normalization
  # factor.
  for(raw_col in names(bulkcounts_df[, -1])){
    bulkvoom_df <- dplyr::bind_cols(bulkvoom_df, tibble::as_tibble(
      bulkcounts_df[raw_col] / norm_factors[raw_col]))
  }

  bulkvoom <- limma::voom(bulkcounts_df[, -1], design=NULL,
                          lib.size=colSums(bulkcounts_df[, -1]),
                          normalize.method='none', plot=F)

  bulkvoom_df <- tibble::as_tibble(bulkvoom$E, .name_repair='minimal')
  bulkvoom_df <- bulkvoom_df %>%
    tibble::add_column(bulkcounts_df[, 1], .before=1)


  bulkvoom_mx <- t(as.matrix(bulkvoom_df[, -1]))

  pca_expr <- prcomp(bulkvoom_mx, scale=TRUE)

  pca_tbl <- tibble::as.tibble(pca_expr$x)
  pca_tbl <- pca_tbl %>%
    tibble::add_column(clinvar_vals, .before=1)

  ggplot(pca_tbl) +
    geom_point(aes(x=PC1, y=PC2, col=var_vals)) +
    theme_classic()

}
