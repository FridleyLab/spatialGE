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
##' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- system.file("extdata", 'melanoma_thrane', package="spatialGE")
#' count_files <- list.files(data_files, full.names=TRUE, pattern='counts')
#' coord_files <- list.files(data_files, full.names=TRUE, pattern='mapping')
#' clin_file <- list.files(data_files, full.names=TRUE, pattern='clinical')
#' melanoma <- STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples
#' summarize_STlist(melanoma)
#'
#' @export summarize_STlist
#'
summarize_STlist = function(x=NULL){
  if(class(x) != 'STlist'){
    stop('The input is not an STlist')
  }

  # Get statistics per sample
  df_tmp = lapply(names(x@counts), function(i){
    total_counts = sum(x@counts[[i]])
    # Count stats
    mean_counts_spot = mean(Matrix::colSums(x@counts[[i]]), na.rm=T)
    min_counts_spot = min(Matrix::colSums(x@counts[[i]]))
    max_counts_spot = max(Matrix::colSums(x@counts[[i]]))
    # Gene stats (non-zero counts)
    mtx_mask = x@counts[[i]] > 0
    mean_genes_spot = mean(Matrix::colSums(mtx_mask))
    min_genes_spot = min(Matrix::colSums(mtx_mask))
    max_genes_spot = max(Matrix::colSums(mtx_mask))

    # Make a row with stats
    df_row_tmp = tibble::tibble(sample_name=i,
                                spotscells=ncol(x@counts[[i]]),
                                genes=nrow(x@counts[[i]]),
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

