##
#' @title filter_data: Filters cells/spots, genes, or samples
#' @description Filtering of spots/cells, genes or samples, as well as count-based
#' filtering
#' @details
#' This function provides options to filter elements in an STlist. It can remove
#' cells/spots or genes based on raw counts (`x@counts`). Users can input an
#' regular expression to query gene names and calculate percentages (for example %
#' mtDNA genes). The function also can filter entire samples. Note that the function
#' removes cells/spots, genes, and/or samples in the raw counts, transformed counts,
#' spatial variables, gene variables, and sample metadata. Also note that the function
#' filters in the following order:
#' 1. Samples (`rm_tissue`)
#' 2. Spots (`rm_spots`)
#' 3. Genes (`rm_genes`)
#' 4. Genes matching `rm_genes_expr`
#' 5. Min and max counts
#'
#' @param x an STlist
#' @param spot_minreads the minimum number of total reads for a spot to be retained
#' @param spot_maxreads the maximum number of total reads for a spot to be retained
#' @param spot_mingenes the minimum number of non-zero counts for a spot to be retained
#' @param spot_maxgenes the maximum number of non-zero counts for a spot to be retained
#' @param spot_minpct the minimum percentage of counts for features defined by `spot_pct_expr` for a spot to be retained.
#' @param spot_maxpct the maximum percentage of counts for features defined by `spot_pct_expr` for a spot to be retained.
#' @param gene_minreads the minimum number of total reads for a gene to be retained
#' @param gene_maxreads the maximum number of total reads for a gene to be retained
#' @param gene_minspots he minimum number of spots with non-zero counts for a gene to be retained
#' @param gene_maxspots the maximum number of spots with non-zero counts for a gene to be retained
#' @param gene_minpct the minimum percentage of spots with non-zero counts for a gene to be retained
#' @param gene_maxpct the maximum percentage of spots with non-zero counts for a gene to be retained
#' @param samples samples (as in `names(x@counts)`) to perform filtering.
#' @param rm_tissue sample (as in `names(x@counts)`) to remove from STlist. Removes samples in `x@counts`, `x@tr_counts`, `x@spatial_meta`, `x@gene_meta`, and `x@sample_meta`
#' @param rm_spots vector of spot/cell IDs to remove. Removes spots/cells in `x@counts`, `x@tr_counts`, and `x@spatial_meta`
#' @param rm_genes vector of gene names to remove from STlist. Removes genes in `x@counts`, `x@tr_counts`, and `x@gene_meta`
#' @param rm_genes_expr a regular expression that matches genes to remove. Removes genes in `x@counts`, `x@tr_counts`, and `x@gene_meta`
#' @param spot_pct_expr a expression to use with `spot_minpct` and `spot_maxpct`. By default '^MT-'.
#' @return an STlist containing the filtered data
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
#'   zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
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
#'   melanoma <- filter_data(melanoma, spot_minreads=2000)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#'
#' @importFrom magrittr %>%
#
# Dispatches through a fixed sequence of named filter stages (see
# `filter_stage_*` in filter_data_stages.R), matching the order documented
# above: samples, spots, genes (by name), genes (by regex), min/max counts.
#
filter_data = function(x=NULL,
                       spot_minreads=0, spot_maxreads=NULL,
                       spot_mingenes=0, spot_maxgenes=NULL,
                       spot_minpct=0, spot_maxpct=NULL,
                       gene_minreads=0, gene_maxreads=NULL,
                       gene_minspots=0, gene_maxspots=NULL,
                       gene_minpct=0, gene_maxpct=NULL,
                       samples=NULL,
                       rm_tissue=NULL,
                       rm_spots=NULL,
                       rm_genes=NULL,
                       rm_genes_expr=NULL,
                       spot_pct_expr="^MT-"){

  # Check than an STList was provided
  assert_stlist(x, message='Please, specify an STList to filter.', check_class=FALSE)

  # 1. Remove entire samples/tissues
  x = filter_stage_samples(x, rm_tissue)

  # Define set of samples to work on
  samples = resolve_samples(x, samples, slot='counts', numeric_handling='plain',
                            char_transform='validate_all', post_check='none',
                            message='Samples in rm_tissue are not present in this STlist.')

  # 2. Remove spots/cells by cell IDs
  x = filter_stage_spots(x, samples, rm_spots)

  # 3. Remove genes by name
  x = filter_stage_genes(x, samples, rm_genes)

  # 4. Remove genes by regular expression
  x = filter_stage_genes_expr(x, samples, rm_genes_expr)

  # 5. Filter data sets using other parameters (min/max counts)
  x = filter_stage_counts(x, samples, spot_minreads, spot_maxreads, spot_mingenes, spot_maxgenes,
                          spot_minpct, spot_maxpct, gene_minreads, gene_maxreads, gene_minspots,
                          gene_maxspots, gene_minpct, gene_maxpct, spot_pct_expr)

  return(x)
}

