# Filter stages for filter_data(), one function per step in the fixed
# filtering order documented in filter_data()'s @details: samples, spots,
# genes (by name), genes (by regular expression), then min/max counts. Each
# function runs one of filter_data()'s original blocks verbatim.

##
# filter_stage_samples
# STEP 1: remove entire samples/tissues (`rm_tissue`).
# @param x an STlist.
# @param rm_tissue sample (as in `names(x@counts)`) to remove from STlist.
# @return the (possibly updated) STlist.
#
filter_stage_samples = function(x, rm_tissue){
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
  return(x)
} # CLOSE filter_stage_samples


##
# filter_stage_spots
# STEP 2: remove spots/cells by cell IDs (`rm_spots`).
# @param x an STlist.
# @param samples the samples (as in `names(x@counts)`) to perform filtering on.
# @param rm_spots vector of spot/cell IDs to remove.
# @return the (possibly updated) STlist.
#
filter_stage_spots = function(x, samples, rm_spots){
  if(!is.null(rm_spots)){
    for(i in samples){
      x@counts[[i]] = x@counts[[i]][, !(colnames(x@counts[[i]]) %in% rm_spots)]
      x@spatial_meta[[i]] = x@spatial_meta[[i]][!(x@spatial_meta[[i]][['libname']] %in% rm_spots), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = x@tr_counts[[i]][, !(colnames(x@tr_counts[[i]]) %in% rm_spots)]
      }
    }
  }
  return(x)
} # CLOSE filter_stage_spots


##
# filter_stage_genes
# STEP 3: remove genes by name (`rm_genes`).
# @param x an STlist.
# @param samples the samples (as in `names(x@counts)`) to perform filtering on.
# @param rm_genes vector of gene names to remove from STlist.
# @return the (possibly updated) STlist.
#
filter_stage_genes = function(x, samples, rm_genes){
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
  return(x)
} # CLOSE filter_stage_genes


##
# filter_stage_genes_expr
# STEP 4: remove genes matching a regular expression (`rm_genes_expr`).
# @param x an STlist.
# @param samples the samples (as in `names(x@counts)`) to perform filtering on.
# @param rm_genes_expr a regular expression that matches genes to remove.
# @return the (possibly updated) STlist.
#
filter_stage_genes_expr = function(x, samples, rm_genes_expr){
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
  return(x)
} # CLOSE filter_stage_genes_expr


##
# filter_stage_counts
# STEP 5: filter spots/cells and genes using min/max total counts, min/max
# number of non-zero genes/spots, and min/max percentage of counts/spots
# matching `spot_pct_expr`. Also drops a sample entirely if all its spots
# were removed.
# @param x an STlist.
# @param samples the samples (as in `names(x@counts)`) to perform filtering on.
# @param spot_minreads,spot_maxreads,spot_mingenes,spot_maxgenes,spot_minpct,spot_maxpct,gene_minreads,gene_maxreads,gene_minspots,gene_maxspots,gene_minpct,gene_maxpct,spot_pct_expr
# see `filter_data()`.
# @return the (possibly updated) STlist.
#
filter_stage_counts = function(x, samples, spot_minreads, spot_maxreads, spot_mingenes, spot_maxgenes,
                               spot_minpct, spot_maxpct, gene_minreads, gene_maxreads, gene_minspots,
                               gene_maxspots, gene_minpct, gene_maxpct, spot_pct_expr){
  # To prevent NOTES in R CMD check
  . = NULL

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
} # CLOSE filter_stage_counts
