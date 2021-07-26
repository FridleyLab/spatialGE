##
#' @title bulk_pca: Performs PCA on "bulk" spatial array
#' @description Perform and plot a PCA of simulated bulk RNA-Seq from spatial
#' transcriptomics data.
#' @details
#' This function takes an STList, and optionally the name of a clinical or sample-associated
#' variable, and performs PCA on simulated bulk RNA-Seq from ST data. The counts are simulated
#' by summing all counts from an array for each gene.
#'
#' @param x an STList.
#' @param clinvar a string indicating the name of the variable in the clinical data.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # bulk_pca(melanoma, clinvar='gender')
#'
#' @export
#
#
bulk_pca <- function(x=NULL, clinvar=NULL, color_pal="muted") {

  if(!is.null(clinvar) && !(clinvar %in% colnames(x@clinical))){
    stop('Variable not in sample data. Verify that input matches variable name in sample data')
  }

  require('magrittr')
  require('ggplot2')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Stop function if only one sample provided.
  if(length(x@counts) < 2){
    stop('Refusing to plot a single observation PCA!')
  }

  # Extract clinical data from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  # Also get labels for PCA points.
  if(!is.null(clinvar)){
    clinvar_vals <- as.character(x@clinical[[clinvar]])
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name=as.character(clinvar))
    pca_labs <- as.character(x@clinical[[1]])
  }else{
    clinvar <- 'ST_sample'
    #clinvar_vals <- as.character(1:length(x@counts))
    clinvar_vals = names(x@counts)
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='ST_sample')
    #pca_labs <- as.character(1:length(x@counts))
    pca_labs = as.character(clinvar_vals[[1]])
  }

  # Create data frame to store "bulk counts".
  bulkcounts_df <- x@counts[[1]][, 1]

  # Loop through count matrices in STList
  for(i in 1:length(x@counts)){

    bulk_expr <- rowSums(x@counts[[i]][, -1])
    bulk_expr <- tibble::as_tibble_col(bulk_expr, column_name=paste0("st_", i))
    bulk_expr <- bulk_expr %>%
      tibble::add_column(x@counts[[i]][, 1], .before=1)

    bulkcounts_df <- dplyr::inner_join(bulkcounts_df, bulk_expr, by=c('gene'))
    }

  # Get normalization factors from the "bulk" libraries.
  norm_factors <- edgeR::calcNormFactors(bulkcounts_df[, -1], method='TMM', lib.size=NULL)

  # Create new object to store the size-corrected counts.
  normcounts_df <- c()

  # Divide each count value by their respective column (sample) normalization factor.
  for(raw_col in names(bulkcounts_df[, -1])){
    normcounts_df <- dplyr::bind_cols(normcounts_df,
                                      tibble::as_tibble(bulkcounts_df[raw_col] / norm_factors[raw_col]))
  }

  # Apply limma-voom to lib-size normalized "bulk" libraries.
  bulkvoom <- limma::voom(normcounts_df, design=NULL,
                          lib.size=colSums(bulkcounts_df[, -1]),
                          normalize.method='none', plot=F)

  # Extract voom-normlized matrix and put back gene names.
  bulkvoom_df <- tibble::as_tibble(bulkvoom$E, .name_repair='minimal')
  bulkvoom_df <- bulkvoom_df %>%
    tibble::add_column(bulkcounts_df[, 1], .before=1)

  # Turn voom-counts to a transposed matrix.
  bulkvoom_mx <- t(as.matrix(bulkvoom_df[, -1]))

  # Perform PCA on transoposed matrix.
  pca_expr <- prcomp(bulkvoom_mx, scale=TRUE)

  # Get PCA coordinates and add clinical variable data.
  pca_tbl <- tibble::as_tibble(pca_expr$x)
  pca_tbl <- pca_tbl %>%
    tibble::add_column(clinvar_vals, .before=1)


  # Get number of categories from selected variable.
  n_cats <- nlevels(as.factor(pca_tbl[[clinvar]]))

  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)

  # Assocuate colors with categories.
  names(cat_cols) <- levels(as.factor(pca_tbl[[clinvar]]))

  # Define shapes of points according to clinical variable.
  cat_shapes <- (16:25)[1:n_cats]
  names(cat_shapes) <- levels(as.factor(pca_tbl[[clinvar]]))

# NOTE: INSTEAD OF NUMBERS, WOULD BE GREAT TO HAVE SAMPLE ID PLOTTED
  ggplot(pca_tbl) +
    geom_point(aes(x=PC1, y=PC2,
                   shape=get(clinvar),
                   color=get(clinvar)),
               size=5) +
    ggrepel::geom_text_repel(aes(x=PC1, y=PC2, label=pca_labs)) +
    scale_color_manual(clinvar, values=cat_cols) +
    scale_shape_manual(clinvar, values=cat_shapes) +
    ggtitle('PCA of "bulk" ST RNA-Seq libraries') +
    coord_fixed() +
    theme_bw()

}
