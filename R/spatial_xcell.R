
require('tidyverse')
require('xCell')
require('janitor')

spatial_xcell <- function(x=NULL){
  # Store name of deconvolution method.
  x@cell_deconv$deconv_method <- 'xCell'

  # Extract count data.
  gene_names <- x@voom_counts$gene
  df <- as.matrix(x@voom_counts[, -1])
  rownames(df) <- gene_names

  # # Filter out genes with low variance.
  # gene_var_mask <- y@gene_stdev$gene_stdevs <= 1.9
  # discarded_genes <- sum(gene_var_mask)
  # df <- df[!gene_var_mask, ]

  # Perform xCell analysis and make cell names as column (no rownames).
  # NOTE: Applies `janitor` to clean cell names.
  df_xcell <- xCellAnalysis(df, rnaseq=T)
  cell_names <- rownames(df_xcell) %>%
    make_clean_names(.)
  df_xcell <- as_tibble(df_xcell) %>%
    add_column(cell_names, .before=1)

  # Remove immune, stroma, and microenvironment xCell scores.
  scores_remove <- grep("immune_score|stroma_score|microenvironment_score",
                         unlist(df_xcell[, 1]))
  df_xcell_NoPurityScores <- df_xcell[-scores_remove, ]

  # Store untransformed and transformed scores.
  x@cell_deconv$deconv_matrix <- df_xcell_NoPurityScores
  df_xcell_sqrt <- sqrt(df_xcell_NoPurityScores[, -1])
  df_xcell_sqrt <- df_xcell_sqrt %>%
    add_column(df_xcell_NoPurityScores[, 1], .before=1)
  x@cell_deconv$transf_deconv_matrix <- df_xcell_sqrt

  # Calculate cell means and standard deviations, and store in object.
  cell_stdevs <- apply(as.data.frame(x@cell_deconv$transf_deconv_matrix[, -1]), 1, sd, na.rm=T)
  cell_means <- rowMeans(x@cell_deconv$transf_deconv_matrix[, -1], na.rm=T)
  cell_stdevs_df <- tibble(cell_means, cell_stdevs) %>%
    add_column(df_xcell_NoPurityScores[, 1], .before=1)
  names(cell_stdevs_df)[1] <- 'cell'
  x@cell_stdev <- cell_stdevs_df

  return(x)
}

