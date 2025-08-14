##
#' @title pseudobulk_samples: Aggregates counts into "pseudo bulk" samples
#' @description Aggregates spot/cell counts into "pseudo bulk" samples for data exploration
#' @details
#' This function takes an STlist and aggregates the spot/cell counts into "pseudo bulk"
#' counts by summing all counts from all cell/spots for each gene. Then performs
#' Principal Component Analysis (PCA) to explore non-spatial sample-to-sample
#' variation
#'
#' @param x an STlist.
#' @param max_var_genes number of most variable genes (standard deviation) to use in
#' pseudobulk analysis
#' @param calc_umap logical, whether to calculate UMAP embeddings in addition to PCs
#' @return an STlist with appended pseudobulk counts and PCA coordinates
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- pseudobulk_samples(melanoma)
#'   pseudobulk_dim_plot(melanoma)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export pseudobulk_samples
#'
#' @import Matrix
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#'
pseudobulk_samples = function(x=NULL, max_var_genes=5000, calc_umap=FALSE){

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STlist')){
    raise_err(err_code='error0018')
  }

  # Stop function if only one sample provided.
  if(length(x@counts) < 4){
    raise_err(err_code='error0017')
  }

  # Need at least two variable genes to create 2 PCs
  if(max_var_genes < 2){
    raise_err(err_code='error0005')
  }

  #Loop through count matrices, get gene-wise sums, and collate samples
  bulkcounts_tmp = lapply(1:length(x@counts), function(i){
    bulk_expr = data.frame(V1=Matrix::rowSums(x@counts[[i]]), gene=rownames(x@counts[[i]]))
    rownames(bulk_expr) = NULL
    colnames(bulk_expr)[1] = names(x@counts)[i]
    return(bulk_expr)
  })
  bulkcounts_df = data.frame(gene=bulkcounts_tmp[[1]][['gene']])
  lapply(1:length(bulkcounts_tmp), function(i){
    bulkcounts_df <<- dplyr::full_join(bulkcounts_df, bulkcounts_tmp[[i]], by='gene')
    return(NULL)
  })
  rm(bulkcounts_tmp) # Clean env

  # Apply log transform
  # Calculate (spot) library sizes. Then, add 1 to each library size.
  libsizes = colSums(bulkcounts_df[, -1], na.rm=T)
  # Check that there are not zero-count spots
  if(any(libsizes == 0)){
    raise_err(err_code='error0019')
  }
  # Divide each count value by their respective column (spot) normalization factor.
  normcounts_df = sweep(bulkcounts_df[, -1], 2, libsizes, '/')
  # Apply log transformation to count data.
  tr_df = log1p(normcounts_df * 100000)
  # Put back gene names to matrix and store in object.
  tr_df = tibble::as_tibble(tr_df) %>%
    tibble::add_column(gene=bulkcounts_df[[1]], .before=1) %>%
    tibble::column_to_rownames('gene')

  rm(bulkcounts_df, normcounts_df, libsizes) # Clean env

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
  vargenes = as.vector(na.omit(names(vargenes[1:max_var_genes])))
  tr_df = tr_df[, match(vargenes, colnames(tr_df))]
  tr_df = scale(tr_df)

  rm(vargenes) # Clean env

  # Save scaled pseudobulk matrix to STlist for plotting
  x@misc[['scaled_pbulk_mtx']] = tr_df

  # Perform PCA on transposed matrix and save PCA coordinates to STlist
  # Also save percentage of explained variance
  pca_expr = prcomp(tr_df, scale=F, center=F)
  pca_expl_var = pca_expr[['sdev']]^2 / sum(pca_expr[['sdev']]^2)
  names(pca_expl_var) = colnames(pca_expr[['x']])

  x@misc[['pbulk_pca']] = as.data.frame(pca_expr[['x']])
  x@misc[['pbulk_pca_var']] = pca_expl_var

  if(calc_umap){
    pcs = 30
    if(ncol(pca_expr[['x']]) < 30){
      pcs = ncol(pca_expr[['x']])
    }
    if(nrow(pca_expr[['x']]) <= 4){
      n_neighbors = 2
    } else{
      n_neighbors = round(pcs/2, 0)
    }
    umap_obj = uwot::umap(pca_expr[['x']][, 1:pcs], pca=NULL, n_neighbors=n_neighbors)
    #rownames(umap_obj) = rownames(merged_cts_scl)
    colnames(umap_obj) = c('UMAP1', 'UMAP2')
    x@misc[['pbulk_umap']] = as.data.frame(umap_obj)
  }

  return(x)
}


##
#' @title pseudobulk_dim_plot: Plot PCA of pseudobulk samples
#' @description Generates a PCA plot after computation of "pseudobulk" counts
#' @details
#' Generates a Principal Components Analysis plot to help in initial data exploration of
#' differences among samples. The points in the plot represent "pseudobulk" samples.
#' This function follows after usage of `pseudobulk_samples`.
#'
#' @param x an STlist with pseudobulk PCA results in the `@misc` slot (generated by
#' `pseudobulk_samples`)
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector of color names or HEX values. Each color represents a category in the
#' variable specified in `plot_meta`
#' @param plot_meta a string indicating the name of the variable in the sample
#' metadata to color points in the PCA plot
#' @param dim one of `umap` or `pca`. The dimension reduction to plot
#' @param pcx integer indicating the principal component to plot in the x axis
#' @param pcy integer indicating the principal component to plot in the y axis
#' @param ptsize the size of the points in the PCA plot. Passed to the `size`
#' aesthetic from `ggplot2`
#' @return a ggplot object
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- pseudobulk_samples(melanoma)
#'   pseudobulk_dim_plot(melanoma, plot_meta='patient')
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export pseudobulk_dim_plot
#'
pseudobulk_dim_plot = function(x=NULL, color_pal='muted', plot_meta=NULL, dim='pca', pcx=1, pcy=2, ptsize=5){

  # Prepare meta data
  meta_df = prepare_pseudobulk_meta(x, plot_meta) %>%
    tibble::rownames_to_column('pca_labs')

  # If no sample metadata column was specified, use 'ST_sample' as generated by `prepare_pseudobulk_meta`
  if(is.null(plot_meta)){
    plot_meta = 'ST_sample'
  }

  # Select PCA or UMAP
  if(tolower(dim) == 'umap'){
    pca_tbl = x@misc[['pbulk_umap']]
  } else{
    pca_tbl = x@misc[['pbulk_pca']]
  }

  # Get PCA coordinates and add clinical variable data.
  pca_tbl = pca_tbl %>%
    tibble::rownames_to_column(var='pca_labs') %>%
    dplyr::full_join(meta_df, by='pca_labs')
  # Get number of categories from selected variable.
  n_cats <- nlevels(as.factor(pca_tbl[[plot_meta]]))
  # Create color palette.
  cat_cols = color_parse(color_pal, n_cats=n_cats)
  # Associate colors with categories.
  names(cat_cols) <- levels(as.factor(pca_tbl[[plot_meta]]))

  # Select columns (PCs/UMAPs) to plot
  if(tolower(dim) == 'umap'){
    pcx = 'UMAP1'
    pcy = 'UMAP2'
    pcx_lab = pcx
    pcy_lab = pcy
    title = 'UMAP'
  } else{
    pcx = grep(paste0('^PC', pcx, '$'), colnames(pca_tbl), value=T)
    pcy = grep(paste0('^PC', pcy, '$'), colnames(pca_tbl), value=T)
    pcx_lab = paste0(pcx, ' (', round(x@misc[['pbulk_pca_var']][pcx], 3) * 100, '%)')
    pcy_lab = paste0(pcy, ' (', round(x@misc[['pbulk_pca_var']][pcy], 3) * 100, '%)')
    title = 'PCA'
  }

  # Make plot
  pca_p = ggplot2::ggplot(pca_tbl) +
    ggplot2::geom_point(ggplot2::aes(x=.data[[pcx]], y=.data[[pcy]], color=.data[[plot_meta]]), size=ptsize) +
    ggrepel::geom_text_repel(ggplot2::aes(x=.data[[pcx]], y=.data[[pcy]], label=pca_labs)) +
    ggplot2::xlab(pcx_lab) +
    ggplot2::ylab(pcy_lab)

  if(plot_meta == 'pca_labs'){
    pca_p = pca_p + ggplot2::labs(color='Sample')
  }

  pca_p = pca_p +
    ggplot2::scale_color_manual(values=cat_cols) +
    #scale_shape_manual(values=cat_shapes) +
    ggplot2::ggtitle(paste0(title, ' of "pseudobulk" samples')) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw()

  return(pca_p)
}


##
#' @title pseudobulk_heatmap: Heatmap of pseudobulk samples
#' @description Generates a heatmap plot after computation of "pseudobulk" counts
#' @details
#' Generates a heatmap of transformed "pseudobulk" counts to help in initial data
#' exploration of differences among samples. Each column in the heatmap represents a
#' "pseudobulk" sample. Rows are genes, with the number of genes displayed controlled by
#' the `hm_display_genes` argument. This function follows after usage of `pseudobulk_samples`.
#'
#' @param x an STlist with pseudobulk counts in the `@misc` slot (generated by
#' `pseudobulk_samples`)
#' @param color_pal a string of a color palette from khroma or RColorBrewer, or a
#' vector of color names or HEX values. Each color represents a category in the
#' variable specified in `plot_meta`
#' @param plot_meta a string indicating the name of the variable in the sample
#' metadata to annotate heatmap columns
#' @param hm_display_genes number of genes to display in heatmap, selected based on
#' decreasing order of standard deviation across samples
#' @return a ggplot object
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- pseudobulk_samples(melanoma)
#'   hm <- pseudobulk_heatmap(melanoma, plot_meta='BRAF_status', hm_display_genes=30)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export pseudobulk_heatmap
#'
pseudobulk_heatmap = function(x=NULL, color_pal='muted', plot_meta=NULL, hm_display_genes=30){
  # Prepare meta data
  meta_df = prepare_pseudobulk_meta(x, plot_meta)

  # If no sample metadata column was specified, use 'ST_sample' as generated by `prepare_pseudobulk_meta`
  if(is.null(plot_meta)){
    plot_meta = 'ST_sample'
  }

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

  # Order genes based on coefficient of variation
  order_genes = sort(apply(hm_mtx, 1, function(i){
    mean_tmp = mean(unlist(i)) + 1e-20 # Add shift to avoid zero division
    sd_tmp = sd(unlist(i))
    cv_tmp = mean_tmp/sd_tmp
  }), decreasing=T)

  hm_mtx = hm_mtx[ names(order_genes)[1:hm_display_genes], ]

  # Make heatmap
  hm_p = ComplexHeatmap::Heatmap(hm_mtx, show_row_dend=F,
                                 top_annotation=hm_ann,
                                 column_names_rot=60,
                                 heatmap_legend_param=list(title="Scaled mean\nexpression"),
                                 row_names_gp=grid::gpar(fontsize=8),
                                 show_column_names=T, column_title='Aggregated gene expression\n("pseudobulk")')

  #hm_p = capture.output(ComplexHeatmap::draw(hm_p, merge_legend=T, padding=unit(c(2, 10, 2, 2), "mm")))
  hm_p = ComplexHeatmap::draw(hm_p, merge_legend=T, padding=unit(c(2, 10, 2, 2), "mm"))

  return(hm_p)
}


# Helpers ----------------------------------------------------------------------

##
# @title prepare_pseudobulk_meta
#
# @param x a STlist
# @param plot_meta a string indicating the name of the variable in the sample metadata
# @return a data frame with sample-level metadata
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
                                   sample_name=x@sample_meta[[1]]) %>%
      tibble::column_to_rownames('sample_name')
  } else{
    #plot_meta <- 'ST_sample'
    clinvar_vals = names(x@counts)
    clinvar_vals <- tibble::tibble(ST_sample=clinvar_vals,
                                   sample_name=clinvar_vals) %>%
      tibble::column_to_rownames('sample_name')
  }
  return(clinvar_vals)
}

