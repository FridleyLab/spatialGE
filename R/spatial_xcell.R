##
# This function applies xCell to the stored normalized matrices in order to obtain
# cell scores for each of the library/spots. The results are stored as untransformed,
# and square-root transformed scores. Stroma/Tumor scores are stored in a different
# slot.
#
# @param x, a STList with normalized count matrices.
# @return x, an updated STList with xCell scores.
#
#
# require('tidyverse')
# require('xCell')
# require('janitor')
require('magrittr')
require('xCell')

spatial_xcell <- function(x=NULL){

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Test that normalized count matrices are available.
  if(length(x@voom_counts) == 0){
    stop("The STList does not have normalized counts.")
  }

  for(i in 1:length(x@voom_counts)){

    # Creates list for a given spatial array.
    x@cell_deconv[[i]] <- list()

    # Store name of deconvolution method.
    x@cell_deconv[[i]]$deconv_method <- 'xCell'

    # Extract count data.
    gene_names <- x@voom_counts[[i]]$gene
    df <- as.matrix(x@voom_counts[[i]][, -1])
    rownames(df) <- gene_names

    # # Filter out genes with low variance.
    # gene_var_mask <- y@gene_stdev$gene_stdevs <= 1.9
    # discarded_genes <- sum(gene_var_mask)
    # df <- df[!gene_var_mask, ]

    # Show progress.
    cat(paste0("Applying xCell to spatial array #", i, "...\n"))

    # Perform xCell analysis and make cell names as column (no rownames).
    # NOTE: Applies `janitor` to clean cell names.
    df_xcell <- xCellAnalysis(df, rnaseq=T)
    cell_names <- rownames(df_xcell) %>%
      janitor::make_clean_names(.)
    df_xcell <- tibble::as_tibble(df_xcell) %>%
      tibble::add_column(cell_names, .before=1)

    # Remove immune, stroma, and microenvironment xCell scores.
    scores_xCell <- grep("immune_score|stroma_score|microenvironment_score",
                          unlist(df_xcell[, 1]))
    df_xcell_NoPurityScores <- df_xcell[-scores_xCell, ]
    df_xcell_PurityScores <- df_xcell[scores_xCell, ]

    # Store untransformed and transformed scores.
    x@cell_deconv[[i]]$deconv_matrix <- df_xcell_NoPurityScores
    df_xcell_sqrt <- sqrt(df_xcell_NoPurityScores[, -1])
    df_xcell_sqrt <- df_xcell_sqrt %>%
      tibble::add_column(df_xcell_NoPurityScores[, 1], .before=1)
    x@cell_deconv[[i]]$transf_deconv_matrix <- df_xcell_sqrt

    # Store transformed tumor/stroma scores.
    df_xcell_purity_sqrt <- sqrt(df_xcell_PurityScores[, -1])
    df_xcell_purity_sqrt <- df_xcell_purity_sqrt %>%
      tibble::add_column(df_xcell_PurityScores[, 1], .before=1)
    x@cell_deconv[[i]]$transf_tumorstroma <- df_xcell_purity_sqrt

    # Calculate cell means and standard deviations, and store in object.
    cell_stdevs <- apply(as.data.frame(
      x@cell_deconv[[i]]$transf_deconv_matrix[, -1]), 1, sd, na.rm=T)
    cell_means <- rowMeans(x@cell_deconv[[i]]$transf_deconv_matrix[, -1], na.rm=T)
    cell_stdevs_df <- tibble::tibble(cell_means, cell_stdevs) %>%
      tibble::add_column(df_xcell_NoPurityScores[, 1], .before=1)
    names(cell_stdevs_df)[1] <- 'cell'
    x@cell_stdev[[i]] <- cell_stdevs_df
  }

  return(x)

}
