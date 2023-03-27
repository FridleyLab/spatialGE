##
#' @title summarize_STlist: Generates a data frame with summary statistics
#' @description Produces a data frame with statistics of the samples in the STlist
#' @details The function creates a table with count and gene statistics for the samples
#' stored in the STlist
#'
#' @param x an STlist
#'
#' @export
#'
#
#
summarize_STlist = function(x=NULL){
  if(class(x) != 'STlist'){
    stop('The input is not an STlist')
  }

  df_tmp = tibble::tibble()
  for(i in names(x@counts)){
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

    df_tmp = dplyr::bind_rows(df_tmp,
                              tibble::tibble(sample_name=i,
                                             mean_counts_per_spotcell=mean_counts_spot,
                                             min_counts_per_spotcell=min_counts_spot,
                                             max_counts_per_spotcell=max_counts_spot,
                                             mean_genes_per_spotcell=mean_genes_spot,
                                             min_genes_per_spotcell=min_genes_spot,
                                             max_genes_per_spotcell=max_genes_spot,))
  }
  return(df_tmp)
}
