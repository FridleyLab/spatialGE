##
#' @title pseudobulk_samples: Aggregates spot/cell counts into "pseudo bulk" samples
#' @description Aggregates counts into "pseudo bulk" samples
#' @details
#' This function takes an STlist and aggregates the spot/cell counts into "pseudo bulk"
#' counts by summing all counts from all cell/spots for each gene.
#'
#' @param x an STList.
#' @param max_var_genes number of most variable genes (standard deviation) to use in pseudobulk analysis
#'
#' @examples
#' # In this example, melanoma is an STlist.
#' # pseudobulk_samples(melanoma, max_var_genes=5000)
#'
#' @export
#'
#' @import Matrix
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#
#
pseudobulk_samples = function(x=NULL, max_var_genes=5000) {

  #require('magrittr')
  #require('Matrix')
  #require('ggplot2')

  tr_method='log'

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be an STList.")
  }

  # Stop function if only one sample provided.
  if(length(x@counts) < 3){
    stop('Refusing to make PCA and heatmap containing less than two samples!')
  }

  # Create data frame to store "bulk counts".
  bulkcounts_df <- tibble::tibble(gene=rownames(x@counts[[1]]))

  # Loop through count matrices, get gene-wise sums, and collate samples
  for(i in names(x@counts)){
    expanded_mtx = expandSparse(x@counts[[i]])

    bulk_expr <- rowSums(expanded_mtx) %>%
      tibble::as_tibble_col(., column_name=i) %>%
      tibble::add_column(gene=rownames(x@counts[[i]]), .before=1)

    bulkcounts_df <- dplyr::inner_join(bulkcounts_df, bulk_expr, by='gene')
  }

  # Apply voom or log transform
  # if(tr_method == 'voom'){
  #   # Get normalization factors from the "bulk" libraries.
  #   norm_factors <- edgeR::calcNormFactors(bulkcounts_df[, -1], method='TMM', lib.size=NULL)
  #   # Create new object to store the size-corrected counts.
  #   normcounts_df <- c()
  #   # Divide each count value by their respective column (sample) normalization factor.
  #   for(raw_col in names(bulkcounts_df[, -1])){
  #     normcounts_df <- dplyr::bind_cols(normcounts_df, tibble::as_tibble(bulkcounts_df[raw_col] / norm_factors[raw_col]))
  #   }
  #   # Apply limma-voom to lib-size normalized "bulk" libraries.
  #   tr_df <- limma::voom(normcounts_df, design=NULL, lib.size=colSums(bulkcounts_df[, -1]), normalize.method='none', plot=F)
  #   # Put back gene names to matrix and store in object.
  #   tr_df = tibble::as_tibble(tr_df$E) %>%
  #     tibble::add_column(gene=bulkcounts_df[[1]], .before=1)
  # } else
  if(tr_method == 'log'){
    # Calculate (spot) library sizes. Then, add 1 to each library size.
    libsizes = colSums(bulkcounts_df[, -1])
    # Check that there are not zero-count spots
    if(any(libsizes == 0)){
      stop('Please, remove samples containing zero reads...')
    }
    # Divide each count value by their respective column (spot) normalization factor.
    normcounts_df = sweep(bulkcounts_df[, -1], 2, libsizes, '/')
    # Then multiply by scaling factor
    # df = df * scale_f
    # Apply log transformation to count data.
    tr_df = log1p(normcounts_df * 100000)
    # Put back gene names to matrix and store in object.
    tr_df = tibble::as_tibble(tr_df) %>%
      tibble::add_column(gene=bulkcounts_df[[1]], .before=1) %>%
      tibble::column_to_rownames('gene')
  }

  # Turn transformed counts to a transposed matrix.
  tr_df <- t(as.matrix(tr_df))

  # Check that number of genes specified by user is lower than genes available across all samples
  min_genes_pseudobulk = ncol(tr_df)
  if(max_var_genes > min_genes_pseudobulk){
    warning(paste0('Not enough genes in the pseudobulk matrix. Setting max_var_genes to ', min_genes_pseudobulk, '.'))
    max_var_genes = min_genes_pseudobulk
  }
  rm(min_genes_pseudobulk) # Clean env

  # Get variable genes and subset
  vargenes = sort(apply(tr_df, 2, sd), decreasing=T)
  tr_df = tr_df[, match(names(vargenes[1:max_var_genes]), colnames(tr_df))]
  tr_df = scale(tr_df)

  # Save scaled pseudobulk matrix to STlist for plotting
  x@misc[['scaled_pbulk_mtx']] = tr_df

  # Perform PCA on transposed matrix and save PCA coordinates to STlist
  # Also save percentage of explained variance
  pca_expr = prcomp(tr_df, scale=F, center=F)
  pca_expl_var = pca_expr[['sdev']]^2 / sum(pca_expr[['sdev']]^2)
  names(pca_expl_var) = colnames(pca_expr[['x']])
  x@misc[['pbulk_pca']] = as.data.frame(pca_expr[['x']])
  x@misc[['pbulk_pca_var']] = pca_expl_var

  return(x)
}


##
#' @title pseudobulk_pca_plot
#'
#' @param x an STList with pseudobulk counts in the `@misc` slot
#' @param plot_meta a string indicating the name of the variable in the sample metadata.
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @param ptsize the size of the points in the PCA plot. Passed to `size` aesthetic from `ggplot2`.
#' @return a ggplot object
#'
#' @export
#
#
pseudobulk_pca_plot = function(x=NULL, color_pal='muted', plot_meta=NULL, pcx=1, pcy=2, ptsize=5){
  # Prepare meta data
  meta_df = prepare_pseudobulk_meta(x, plot_meta) %>%
    tibble::rownames_to_column('pca_labs')
  # Get PCA coordinates and add clinical variable data.
  pca_tbl = x@misc[['pbulk_pca']] %>%
    tibble::rownames_to_column(var='pca_labs') %>%
    dplyr::full_join(meta_df, by='pca_labs')
  # Get number of categories from selected variable.
  n_cats <- nlevels(as.factor(pca_tbl[[plot_meta]]))
  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)
  # Associate colors with categories.
  names(cat_cols) <- levels(as.factor(pca_tbl[[plot_meta]]))
  # Define shapes of points according to clinical variable.
  cat_shapes <- (16:25)[1:n_cats]
  names(cat_shapes) <- levels(as.factor(pca_tbl[[plot_meta]]))

  pcx = grep(paste0('^PC', pcx, '$'), colnames(pca_tbl), value=T)
  pcy = grep(paste0('^PC', pcy, '$'), colnames(pca_tbl), value=T)

  # Make plot
  pca_p = ggplot(pca_tbl) +
    geom_point(aes(x=.data[[pcx]], y=.data[[pcy]], shape=.data[[plot_meta]], color=.data[[plot_meta]]), size=ptsize) +
    ggrepel::geom_text_repel(aes(x=.data[[pcx]], y=.data[[pcy]], label=pca_labs)) +
    xlab(paste0(pcx, ' (', round(x@misc[['pbulk_pca_var']][pcx], 3) * 100, '%)')) +
    ylab(paste0(pcy, ' (', round(x@misc[['pbulk_pca_var']][pcy], 3) * 100, '%)')) +
    scale_color_manual(values=cat_cols) +
    scale_shape_manual(values=cat_shapes) +
    ggtitle('PCA of "pseudobulk" samples') +
    coord_fixed() +
    theme_bw()

  return(pca_p)
}


##
#' @title pseudobulk_heatmap
#'
#' @param x an STList with pseudobulk counts in the `@misc` slot
#' @param plot_meta a string indicating the name of the variable in the sample metadata.
#' @param hm_display_genes number of genes to display in heatmap. Must be equal or lower than `max_var_genes`
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector with colors with enough elements to plot categories.
#' @return a ComplexHeatmap object
#'
#' @export
#
#
pseudobulk_heatmap = function(x=NULL, color_pal='muted', plot_meta=NULL, hm_display_genes=30){
  # Prepare meta data
  meta_df = prepare_pseudobulk_meta(x, plot_meta)
  # Get number of categories from selected variable.
  n_cats = nlevels(as.factor(meta_df[[plot_meta]]))
  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)
  # Associate colors with categories.
  names(cat_cols) = levels(as.factor(meta_df[[plot_meta]]))
  # Make list to use in heatmap annotation
  cat_cols = list(cat_cols)
  names(cat_cols) = plot_meta

  # Create annotation object for ComplexHeatmap
  hm_ann = ComplexHeatmap::HeatmapAnnotation(df=meta_df, col=cat_cols)

  # Order samples according to annotation and subset genes
  hm_mtx = t(x@misc[['scaled_pbulk_mtx']])
  hm_mtx = hm_mtx[, match(rownames(meta_df), colnames(hm_mtx))]
  hm_mtx = hm_mtx[1:hm_display_genes, ]
  # Make heatmap
  hm_p = ComplexHeatmap::Heatmap(hm_mtx, show_row_dend=F,
                                 top_annotation=hm_ann,
                                 column_names_rot=60,
                                 heatmap_legend_param=list(title="Scaled mean\nexpression"),
                                 row_names_gp=grid::gpar(fontsize=8),
                                 show_column_names=T, column_title='Aggregated gene expression\n("pseudobulk")')

  hm_p = ComplexHeatmap::draw(hm_p, merge_legend=T, padding=unit(c(2, 10, 2, 2), "mm"))

  return(hm_p)
}


# Helpers ----------------------------------------------------------------------

##
# @title prepare_pseudobulk_meta
#
# @param x a STlist
# @param plot_meta a string indicating the name of the variable in the sample metadata.
# @return a data frame with sample-level metadata
#
#
prepare_pseudobulk_meta = function(x=NULL, plot_meta=NULL){
  if(!is.null(plot_meta) && !(plot_meta %in% colnames(x@sample_meta))){
    stop('Variable not in sample metadata. Verify that input matches variable name in sample data')
  }

  # Extract clinical/metadata from specified variable. If none specified, use the
  # array IDs from the first column of clinical data.
  # Also get labels for PCA points.
  if(!is.null(plot_meta)){
    clinvar_vals <- as.character(x@sample_meta[[plot_meta]])
    clinvar_vals <- tibble::tibble(!!plot_meta := clinvar_vals,
                                   sample_name=x@sample_meta[['samplename']]) %>%
      tibble::column_to_rownames('sample_name')
  } else{
    plot_meta <- 'ST_sample'
    clinvar_vals = names(x@counts)
    clinvar_vals <- tibble::tibble(ST_sample=clinvar_vals,
                                   sample_name=clinvar_vals) %>%
      tibble::column_to_rownames('sample_name')
  }
  return(clinvar_vals)
}

