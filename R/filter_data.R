##
#' @title filter_data: Removes spots, genes, or samples
#' @description Provides methods to remove spots and genes or entire samples.
#' @details
#' This function provides methods to filter an STList. It can remove spots or genes
#' based on absolute counts or percentages. Users can input a expression to be matched
#' with gene names to calculate percentages (for example %mtDNA genes).
#' Note that the function removes counts raw counts and spatial information. Transformed
#' counts or other results within the STList remain unaffected.
#'
#' @param x an STList
#' @param spot_minreads the minimum number of total reads for a spot to be retained
#' @param spot_maxreads the maximum number of total reads for a spot to be retained
#' @param spot_mingenes the minimum number of non-zero counts for a spot to be retained
#' @param spot_maxgenes the maximum number of non-zero counts for a spot to be retained
#' @param spot_minpct the minimum percentage of counts for features defined by `spot_pctexpr`
#' for a spot to be retained.
#' @param spot_maxpct the maximum percentage of counts for features defined by `spot_pctexpr`
#' for a spot to be retained.
#' @param gene_minreads the minimum number of total reads for a gene to be retained
#' @param gene_maxreads the maximum number of total reads for a gene to be retained
#' @param gene_minspots he minimum number of spots with non-zero counts for a gene to be retained
#' @param gene_maxspots the maximum number of spots with non-zero counts for a gene to be retained
#' @param gene_minpct the minimum percentage of spots with non-zero counts for a gene to be retained
#' @param gene_maxpct the maximum percentage of spots with non-zero counts for a gene to be retained
#' @param who vector of integers (as in `x@counts`) or sample names to perform filtering
#' @param rm_tissue vector of integers (as in `x@counts`) or sample names to remove from STList
#' @param rm_genes vector of gene names to remove from STList
#' @param rm_partialgenes a expression to match gene names. Useful to remove entire gene classes (for example mtDNA '^MT-')
#' @param spot_pctExpr a expression to use with `spot_minpct`. By default '^MT-'.
#'
#' @export
#'
#' @importFrom magrittr %>%
#
#
filter_data = function(x=NULL,
                       spot_minreads=NULL, spot_maxreads=NULL,
                       spot_mingenes=NULL, spot_maxgenes=NULL,
                       spot_minpct=NULL, spot_maxpct=NULL,
                       gene_minreads=NULL, gene_maxreads=NULL,
                       gene_minspots=NULL, gene_maxspots=NULL,
                       gene_minpct=NULL, gene_maxpct=NULL,
                       who=NULL,
                       rm_tissue=NULL,
                       rm_genes=NULL,
                       rm_partialgenes=NULL,
                       spot_pctExpr="^MT-"){

  #require('magrittr')

  # Check than an STList was provided
  if(is.null(x)){
    stop('Please, specify an STList to filter.')
  }

  # Set minimums for arguments if NULL
  if(is.null(spot_minreads)){
    spot_minreads = 0
  }
  if(is.null(spot_mingenes)){
    spot_mingenes = 0
  }
  if(is.null(spot_minpct)){
    spot_minpct = 0
  }
  if(is.null(gene_minreads)){
    gene_minreads = 0
  }
  if(is.null(gene_minspots)){
    gene_minspots = 0
  }
  if(is.null(gene_minpct)){
    gene_minpct = 0
  }

  # Remove entire samples/tissues
  if(!is.null(rm_tissue)){
    if(is.numeric(rm_tissue)){
      samplename = names(x@counts)[rm_tissue]
      x@counts = x@counts[-rm_tissue]
      x@coords = x@coords[-rm_tissue]
      if(nrow(x@clinical) != 0){
        x@clinical = x@clinical[!(x@clinical[[1]] %in% samplename), ]
      }
    }  else if(is.character(rm_tissue)){
      x@counts = x@counts[!grepl(paste0(rm_tissue, collapse="|"), names(x@counts))]
      x@coords = x@coords[!grepl(paste0(rm_tissue, collapse="|"), names(x@counts))]
      if(nrow(x@clinical) != 0){
        x@clinical = x@clinical[!(x@clinical[[1]] %in% rm_tissue), ]
      }
    } else{
      stop("Could not find the specified samples in the STList.")
    }
  }

  # Define set of samples to work on
  if(is.null(who)){
    who = names(x@counts)
  }else if(is.character(who)){
    who = grep(paste0('^', who, '$', collapse="|"), names(x@counts))
  }else if(is.numeric(who)){
    who = names(x@counts[who])
  }

  # Remove genes by name
  if(!is.null(rm_genes)){
    for(i in who){
      x@counts[[i]] = x@counts[[i]][-grep(paste0('^', rm_genes, '$', collapse='|'), rownames(x@counts[[i]])), ]
    }
  }

  # Remove genes by partial name match
  if(!is.null(rm_partialgenes)){
    for(i in who){
      x@counts[[i]] = x@counts[[i]][-grep(paste0(rm_genes, collapse='|'), rownames(x@counts[[i]])), ]
    }
  }

  # Filter data sets using other parameters
  for(i in who){
    # Decompress counts and create mask of zero-counts
    df = expandSparse(x@counts[[i]])

    # Get maximums for spot filter arguments if NULL
    col_total_reads = colSums(df)
    nonzero_df = (df != 0)
    col_expr_reads = colSums(df[grep(spot_pctExpr, rownames(df)), ])
    col_expr_percent = col_expr_reads/col_total_reads
    if(is.null(spot_maxreads)){
      spot_maxreads_tmp = max(col_total_reads)
    } else{
      spot_maxreads_tmp = spot_maxreads
    }
    col_total_genes = colSums(nonzero_df)
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

    # Get maximums for gene filter arguments if NULL
    row_total_reads = rowSums(df)
    row_total_spots = rowSums(nonzero_df)
    row_expr_percent = row_total_spots/ncol(df)
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
    df = df[gene_mask, spot_mask]

    # Deposit data in STList
    df = tibble::rownames_to_column(df, var='gene')
    x@counts[[i]] = makeSparse(df)
    x@coords[[i]] = x@coords[[i]][x@coords[[i]]$libname %in% colnames(df), ]
  }
  return(x)
}

