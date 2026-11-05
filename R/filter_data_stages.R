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
      x@counts[[i]] = sample_counts(x, i)[, !(sample_cells(x, i, 'counts') %in% rm_spots)]
      x@spatial_meta[[i]] = x@spatial_meta[[i]][!(x@spatial_meta[[i]][['libname']] %in% rm_spots), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = sample_tr_counts(x, i)[, !(sample_cells(x, i, 'tr_counts') %in% rm_spots)]
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
      x@counts[[i]] = sample_counts(x, i)[!(sample_genes(x, i, 'counts') %in% rm_genes), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = sample_tr_counts(x, i)[!(sample_genes(x, i, 'tr_counts') %in% rm_genes), ]
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
      x@counts[[i]] = sample_counts(x, i)[!grepl(rm_genes_expr, sample_genes(x, i, 'counts')), ]
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = sample_tr_counts(x, i)[!grepl(rm_genes_expr, sample_genes(x, i, 'tr_counts')), ]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta[[i]] = x@gene_meta[[i]][!grepl(rm_genes_expr, x@gene_meta[[i]][['gene']]), ]
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
    # Every aggregate below runs on the sparse matrix directly.
    #
    # This used to open with expandSparse(x@counts[[i]]) -- a full dense
    # data.frame -- and then build `nonzero_df`, a SECOND full dense logical
    # matrix, so peak footprint was genes x cells x (8 + 4) bytes per sample and,
    # under mclapply, per fork. That is ~200 MB for a 200 x 50000 Xenium sample.
    # Every consumer before the filter is a row or column aggregate, so none of
    # the densification was needed.
    #
    # Exact rather than approximate: @counts holds integer-valued counts, every
    # aggregate here is a sum of them, values sit far below 2^53 so double
    # arithmetic is exact, and skipping stored zeros cannot change a sum
    # (x + 0 == x). Confirmed by the golden gate passing with no regeneration.
    mtx_tmp = sample_counts(x, i)

    # Non-zero mask stays sparse: dgCMatrix != 0 yields an lgCMatrix, not a
    # dense logical matrix. Note Matrix::colSums on an lgCMatrix returns integer
    # where base colSums on a dense logical matrix returns double -- harmless
    # here, since col_total_genes/row_total_spots only feed comparisons and
    # max(), and neither is stored in the output (contrast
    # stlist_process.R's deliberate as.double(), where the type IS contract).
    nonzero_tmp = (mtx_tmp != 0)

    # Calculate total spot/cell counts and zero-count mask for number of genes per cell/spot
    col_total_reads = Matrix::colSums(mtx_tmp)
    col_total_genes = Matrix::colSums(nonzero_tmp)

    # Calculate percentage of counts from genes matching `spot_pct_expr`.
    # Logical mask, never character indexing, per fidelity rule 9.
    col_expr_reads = Matrix::colSums(mtx_tmp[grepl(spot_pct_expr, rownames(mtx_tmp)), , drop=FALSE])
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
    row_total_reads = Matrix::rowSums(mtx_tmp)

    # Calculate percentage of non-zero spots/cells per gene
    row_total_spots = Matrix::rowSums(nonzero_tmp)
    row_expr_percent = row_total_spots/ncol(mtx_tmp)
    rm(nonzero_tmp) # Clean env

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

    # Perform filter. drop=FALSE matters: the previous data.frame version
    # collapsed to a vector when exactly one spot survived, which then broke
    # rownames_to_column() -- so filter_data() could not produce a single-spot
    # sample at all. Staying sparse removes that edge case rather than
    # preserving it.
    mtx_tmp = mtx_tmp[gene_mask, spot_mask, drop=FALSE]

    # Subset data in STlist. Already a dgCMatrix, so no makeSparse() round trip
    # through a dense matrix. The previous path went sparse -> dense -> data.frame
    # -> as.matrix -> sparse, which also meant a duplicate gene name raised an
    # incidental error from data.frame(); that error is gone now. It was never a
    # designed check (fidelity rule 9 records it as incidental), and @counts
    # tolerates duplicate rownames, but it does mean duplicate genes now pass
    # through here silently.
    if(ncol(mtx_tmp) > 0){ # At least one spot
      x@counts[[i]] = mtx_tmp
      x@spatial_meta[[i]] = x@spatial_meta[[i]][x@spatial_meta[[i]]$libname %in% sample_cells(x, i, 'counts'), ] %>%
        dplyr::select(-total_counts, -total_genes) %>%
        dplyr::left_join(., tibble::tibble(libname=sample_cells(x, i, 'counts'),
                                           total_counts=Matrix::colSums(sample_counts(x, i)),
                                           total_genes=Matrix::colSums(sample_counts(x, i) != 0)), by='libname')
      if(!rlang::is_empty(x@tr_counts)){
        x@tr_counts[[i]] = sample_tr_counts(x, i)[sample_genes(x, i, 'tr_counts') %in% sample_genes(x, i, 'counts'),
                                                sample_cells(x, i, 'tr_counts') %in% sample_cells(x, i, 'counts')]
      }
      if(!rlang::is_empty(x@gene_meta)){
        x@gene_meta[[i]] = x@gene_meta[[i]] %>%
          dplyr::select(-gene_mean, -gene_stdevs) %>%
          dplyr::left_join(tibble::tibble(gene_mean=as.vector(Matrix::rowMeans(sample_tr_counts(x, i))),
                                          gene_stdevs=as.vector(sparse_row_sds(sample_tr_counts(x, i))),
                                          gene=sample_genes(x, i, 'tr_counts')),. , by='gene') %>%
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
