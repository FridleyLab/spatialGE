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
#' @param tr_method one of 'log' or 'voom'. The data transformation to use.
#' @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # bulk_pca(melanoma, clinvar='gender')
#'
#' @export
#'
#' @import Matrix
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
#
bulk_pca <- function(x=NULL, clinvar=NULL, color_pal="muted", tr_method='log', ptsize=5) {

  #require('magrittr')
  #require('Matrix')
  #require('ggplot2')

  if(!is.null(clinvar) && !(clinvar %in% colnames(x@clinical))){
    stop('Variable not in sample data. Verify that input matches variable name in sample data')
  }

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
    clinvar_vals = names(x@counts)
    clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='ST_sample')
    pca_labs = as.character(clinvar_vals[[1]])
  }

  # Create data frame to store "bulk counts".
  bulkcounts_df <- tibble::tibble(gene=rownames(x@counts[[1]]))

  # Loop through count matrices, get gene-wise sums, and collate samples
  for(i in 1:length(x@counts)){
    expanded_mtx = expandSparse(x@counts[[i]])

    bulk_expr <- rowSums(expanded_mtx) %>%
      tibble::as_tibble_col(., column_name=paste0("st_", i)) %>%
      tibble::add_column(gene=rownames(x@counts[[i]]), .before=1)

    bulkcounts_df <- dplyr::inner_join(bulkcounts_df, bulk_expr, by='gene')
  }

  # Apply voom or log transform
  if(tr_method == 'voom'){
    # Get normalization factors from the "bulk" libraries.
    norm_factors <- edgeR::calcNormFactors(bulkcounts_df[, -1], method='TMM', lib.size=NULL)
    # Create new object to store the size-corrected counts.
    normcounts_df <- c()
    # Divide each count value by their respective column (sample) normalization factor.
    for(raw_col in names(bulkcounts_df[, -1])){
      normcounts_df <- dplyr::bind_cols(normcounts_df, tibble::as_tibble(bulkcounts_df[raw_col] / norm_factors[raw_col]))
    }
    # Apply limma-voom to lib-size normalized "bulk" libraries.
    tr_df <- limma::voom(normcounts_df, design=NULL, lib.size=colSums(bulkcounts_df[, -1]), normalize.method='none', plot=F)
    # Put back gene names to matrix and store in object.
    tr_df = tibble::as_tibble(tr_df$E) %>%
      tibble::add_column(gene=bulkcounts_df[[1]], .before=1)
  } else if(tr_method == 'log'){
    # Calculate (spot) library sizes. Then, add 1 to each library size.
    libsizes = colSums(bulkcounts_df[, -1])
    # Check that there are not zero-count spots
    if(any(libsizes == 0)){
      stop('Please, remove spots containing zero reads...')
    }
    # Divide each count value by their respective column (spot) normalization factor.
    normcounts_df = sweep(bulkcounts_df[, -1], 2, libsizes, '/')
    # Then multiply by scaling factor
    # df = df * scale_f
    # Apply log transformation to count data.
    tr_df = log1p(normcounts_df)
    # Put back gene names to matrix and store in object.
    tr_df = tibble::as_tibble(tr_df) %>%
      tibble::add_column(gene=bulkcounts_df[[1]], .before=1)
  }

  # Turn transformed counts to a transposed matrix.
  tr_df <- t(as.matrix(tr_df[, -1]))

  # Perform PCA on transoposed matrix.
  pca_expr <- prcomp(tr_df, scale=TRUE)

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
    geom_point(aes(x=PC1, y=PC2, shape=get(clinvar), color=get(clinvar)), size=ptsize) +
    ggrepel::geom_text_repel(aes(x=PC1, y=PC2, label=pca_labs)) +
    scale_color_manual(clinvar, values=cat_cols) +
    scale_shape_manual(clinvar, values=cat_shapes) +
    ggtitle('PCA of "bulk" ST samples') +
    coord_fixed() +
    theme_bw()

}
