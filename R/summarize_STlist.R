##
#' @title summarize_STlist: Generates a data frame with summary statistics
#' @description Produces a data frame with counts per gene and counts per ROI/spot/cell
#' @details
#' The function creates a table with counts per gene and counts per region of interest (ROI),
#' spot, or cell in the samples stored in the STlist
#'
#' @param x an STlist
#' @return a data frame
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
#'                      samples=clin_file) # Only first two samples
#'   summarize_STlist(melanoma)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?.")
#'   return(NULL)
#' })
#' }
#'
#' @export summarize_STlist
#
summarize_STlist = function(x=NULL){
  assert_stlist(x, message='The input is not an STlist', check_null=FALSE)

  # Get statistics per sample
  df_tmp = lapply(names(x@counts), function(i){
    # Read the matrix once, via the accessor, and compute each column sum once.
    # This previously touched x@counts[[i]] five times and recomputed
    # Matrix::colSums() three times over for mean/min/max -- wasted work in
    # memory, and three separate file reads once the slot holds a manifest.
    mtx_tmp = sample_counts(x, i)
    col_counts = Matrix::colSums(mtx_tmp)
    col_genes = Matrix::colSums(mtx_tmp > 0)

    total_counts = sum(mtx_tmp)
    # Count stats
    mean_counts_spot = mean(col_counts, na.rm=T)
    min_counts_spot = min(col_counts)
    max_counts_spot = max(col_counts)
    # Gene stats (non-zero counts)
    mean_genes_spot = mean(col_genes)
    min_genes_spot = min(col_genes)
    max_genes_spot = max(col_genes)

    # Make a row with stats
    df_row_tmp = tibble::tibble(sample_name=i,
                                spotscells=sample_n_cells(x, i, 'counts'),
                                genes=sample_n_genes(x, i, 'counts'),
                                min_counts_per_spotcell=min_counts_spot,
                                mean_counts_per_spotcell=mean_counts_spot,
                                max_counts_per_spotcell=max_counts_spot,
                                min_genes_per_spotcell=min_genes_spot,
                                mean_genes_per_spotcell=mean_genes_spot,
                                max_genes_per_spotcell=max_genes_spot)
    return(df_row_tmp)
  })

  # Compile all rows into table
  df_tmp = dplyr::bind_rows(df_tmp)

  return(df_tmp)
}

