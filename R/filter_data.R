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

  # To prevent NOTES in R CMD check
  . = NULL

  # Check than an STList was provided
  if(is.null(x)){
    stop('Please, specify an STList to filter.')
  }

  # Remove entire samples/tissues
  if(!is.null(rm_tissue)){
    # Detect if sample entered as index number or sample name
    if(is.numeric(rm_tissue)){
      sample_id = names(x@counts)[rm_tissue]
    }  else if(is.character(rm_tissue)){
      sample_id = rm_tissue
    }

    # Remove samples
    if(all(sample_id %in% names(x@counts))){
      x@counts = x@counts[names(x@counts)[!(names(x@counts) %in% sample_id)]]
      x@spatial_meta = x@spatial_meta[names(x@spatial_meta)[!(names(x@spatial_meta) %in% sample_id)]]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts = x@tr_counts[names(x@tr_counts)[!(names(x@tr_counts) %in% sample_id)]]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta = x@gene_meta[names(x@gene_meta)[!(names(x@gene_meta) %in% sample_id)]]
      }
      if(nrow(x@sample_meta) != 0){
        x@sample_meta = x@sample_meta[!(x@sample_meta[[1]] %in% sample_id), ]
      }
    } else{
      stop('Samples in rm_tissue are not present in this STlist.')
    }

    rm(sample_id) # Clean env
  }

  # Define set of samples to work on
  if(is.null(samples)){
    samples = names(x@counts)
  }
  if(is.character(samples)){
    if(!all(samples %in% names(x@counts))){
      stop('Samples in rm_tissue are not present in this STlist.')
    }
  } else if(is.numeric(samples)){
    samples = names(x@counts)[samples]
  }

  # Remove spots/cells by cell IDs
  if(!is.null(rm_spots)){
    for(i in samples){
      x@counts[[i]] = x@counts[[i]][, !(colnames(x@counts[[i]]) %in% rm_spots)]
      x@spatial_meta[[i]] = x@spatial_meta[[i]][!(x@spatial_meta[[i]][['libname']] %in% rm_spots), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = x@tr_counts[[i]][, !(colnames(x@tr_counts[[i]]) %in% rm_spots)]
      }
    }
  }

  # Remove genes by name
  if(!is.null(rm_genes)){
    for(i in samples){
      x@counts[[i]] = x@counts[[i]][!(rownames(x@counts[[i]]) %in% rm_genes), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = x@tr_counts[[i]][!(rownames(x@tr_counts[[i]]) %in% rm_genes), ]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta[[i]] = x@gene_meta[[i]][!(x@gene_meta[[i]][['gene']] %in% rm_genes), ]
      }
    }
  }

  # Remove genes by regular expression
  if(!is.null(rm_genes_expr)){
    for(i in samples){
      x@counts[[i]] = x@counts[[i]][!grepl(rm_genes_expr, rownames(x@counts[[i]])), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = x@tr_counts[[i]][!grepl(rm_genes_expr, rownames(x@tr_counts[[i]])), ]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta[[i]] = x@gene_meta[[i]][!grepl(rm_genes_expr, rownames(x@gene_meta[[i]])),, ]
      }
    }
  }

  # Filter data sets using other parameters
  for(i in samples){
    # Decompress counts
    df_tmp = expandSparse(x@counts[[i]])

    # Get mask of non-zero counts
    nonzero_df = (df_tmp != 0)

    # Calculate total spot/cell counts and zero-count mask for number of genes per cell/spot
    col_total_reads = colSums(df_tmp)
    col_total_genes = colSums(nonzero_df)

    # Calculate percentage of counts from genes matching `spot_pct_expr`
    col_expr_reads = colSums(df_tmp[grepl(spot_pct_expr, rownames(df_tmp)), ])
    col_expr_percent = col_expr_reads/col_total_reads
    rm(col_expr_reads) # Clean env
    # Force NaNs tp zero, which result from zero as denominator (i.e., zero counts in the cell)
    col_expr_percent[is.na(col_expr_percent)] = 0

    # If no maximum counts set by user, then make it the max counts of each spot/cell
    if(is.null(spot_maxreads)){
      spot_maxreads_tmp = max(col_total_reads)
    } else{
      spot_maxreads_tmp = spot_maxreads
    }
    if(is.null(spot_maxgenes)){
      spot_maxgenes_tmp = max(col_total_genes)
    } else{
      spot_maxgenes_tmp = spot_maxgenes
    }
    if(is.null(spot_maxpct)){
      spot_maxpct_tmp = 1
    } else{
      spot_maxpct_tmp = spot_maxpct
    }

    # Get logical mask for spot filtering
    spot_mask = (col_total_reads >= spot_minreads & col_total_reads <= spot_maxreads_tmp) &
      (col_total_genes >= spot_mingenes & col_total_genes <= spot_maxgenes_tmp) &
      (col_expr_percent >= spot_minpct & col_expr_percent <= spot_maxpct_tmp)

    # Calculate total counts per gene
    row_total_reads = rowSums(df_tmp)

    # Calculate percentage of non-zero spots/cells per gene
    row_total_spots = rowSums(nonzero_df)
    row_expr_percent = row_total_spots/ncol(df_tmp)
    rm(nonzero_df) # Clean env

    # If no maximum counts set by user, then make it the max counts of each gene
    if(is.null(gene_maxreads)){
      gene_maxreads_tmp = max(row_total_reads)
    } else{
      gene_maxreads_tmp = gene_maxreads
    }
    if(is.null(gene_maxspots)){
      gene_maxspots_tmp = max(row_total_spots)
    } else{
      gene_maxspots_tmp = gene_maxspots
    }
    if(is.null(gene_maxpct)){
      gene_maxpct_tmp = 1
    } else{
      gene_maxpct_tmp = gene_maxpct
    }

    # Get logical mask for gene filtering
    gene_mask = (row_total_reads >= gene_minreads & row_total_reads <= gene_maxreads_tmp) &
      (row_total_spots >= gene_minspots & row_total_spots <= gene_maxspots_tmp) &
      (row_expr_percent >= gene_minpct & row_expr_percent <= gene_maxpct_tmp)

    # Perform filter
    df_tmp = df_tmp[gene_mask, spot_mask]

    # Subset data in STlist
    df_tmp = tibble::rownames_to_column(df_tmp, var='gene')
    if(ncol(df_tmp) > 1){ # At least one spot
      x@counts[[i]] = makeSparse(df_tmp)
      x@spatial_meta[[i]] = x@spatial_meta[[i]][x@spatial_meta[[i]]$libname %in% colnames(x@counts[[i]]), ] %>%
        dplyr::select(-total_counts, -total_genes) %>%
        dplyr::left_join(., tibble::tibble(libname=colnames(x@counts[[i]]),
                                           total_counts=Matrix::colSums(x@counts[[i]]),
                                           total_genes=Matrix::colSums(x@counts[[i]] != 0)), by='libname')
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = x@tr_counts[[i]][rownames(x@tr_counts[[i]]) %in% rownames(x@counts[[i]]),
                                            colnames(x@tr_counts[[i]]) %in% colnames(x@counts[[i]])]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta[[i]] = x@gene_meta[[i]] %>%
          dplyr::select(-gene_mean, -gene_stdevs) %>%
          dplyr::left_join(tibble::tibble(gene_mean=as.vector(Matrix::rowMeans(x@tr_counts[[i]])),
                                          gene_stdevs=as.vector(apply(x@tr_counts[[i]], 1, sd)),
                                          gene=rownames(x@tr_counts[[i]])),. , by='gene') %>%
          dplyr::relocate(gene_mean, gene_stdevs, .after='gene')
      }
    } else{
      warning(paste0('All spots from ', i, ' were removed. Removing sample from STlist...'))
      x@counts = x@counts[!grepl(i, names(x@counts))]
      x@spatial_meta = x@spatial_meta[!grepl(i, names(x@spatial_meta))]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts = x@tr_counts[!grepl(i, names(x@tr_counts))]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta = x@gene_meta[!grepl(i, names(x@gene_meta))]
      }
    }
  }
  return(x)
}

