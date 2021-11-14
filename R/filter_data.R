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
#' @param spot_minreads @param spot_maxreads the minimum (or maximum) number of total reads for a spot to be retained
#' @param spot_mingenes @param spot_maxgenes the minimum (or maximum) number of non-zero counts for a spot to be retained
#' @param spot_minpct @param spot_maxpctExpr the minimum (or maximum) percentage of counts for features defined by `spot_pctexpr`
#' for a spot to be retained. By default, mtDNA genes.
#' @param gene_minreads @param gene_maxreads the minimum (or maximum) number of total reads for a gene to be retained
#' @param gene_minspots @param gene_maxspots the minimum (or maximum) number of spots with non-zero counts for a gene to be retained
#' @param gene_minpct @param gene_maxpct the minimum (or maximum) percentage of spots with non-zero counts for a gene to be retained
#' @param who vector of integers (as in `x@counts`) or sample names to perform filtering
#' @param rm_tissue vector of integers (as in `x@counts`) or sample names to remove from STList
#' @param rm_genes vector of gene names to remove from STList
#' @param rm_partialgenes a expression to match gene names. Useful to remove entire gene classes (for example mtDNA '^MT-')
#' @param spot_pctExpr a expression to use with `spot_minpct`. By default '^MT-'.
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

  require('magrittr')

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

  # Define set of samples to work on
  if(is.null(who)){
    who = 1:length(x@counts)
  }
  if(is.character(who)){
    who = grep(paste0('^', who, '$', collapse="|"), names(x@counts))
  }

  # Remove entire samples/tissues
  if(!is.null(rm_tissue)){
    if(is.numeric(rm_tissue)){
      x@counts = x@counts[-rm_tissue]
      x@coords = x@coords[-rm_tissue]
    }  else if(is.character(rm_tissue)){
      x@counts = x@counts[!grepl(paste0(rm_tissue, collapse="|"), names(x@counts))]
      x@coords = x@coords[!grepl(paste0(rm_tissue, collapse="|"), names(x@counts))]
    } else{
      stop("Could not find the specified samples in the STList.")
    }
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
    nonzero_df = df != 0

    # SPOT-WISE FILTER
    # Get total read and gene counts per spot
    col_total_reads = colSums(df)
    col_total_genes = colSums(nonzero_df)
    # Get percentage of genes matching expression by spots
    col_expr_reads = colSums(df[grep(spot_pctExpr, rownames(df)), ])
    col_expr_percent = col_expr_reads/col_total_reads

    # GENE-WISE FILTER
    # Get total read and spot counts per gene
    row_total_reads = rowSums(df)
    row_total_spots = rowSums(nonzero_df)
    # Get percentage of spots for each gene
    row_expr_percent = row_total_spots/ncol(df)

    # Set maximums for SPOT-WISE filter
    if(is.null(spot_maxreads)){
      spot_maxreads = max(col_total_reads)
    }
    if(is.null(spot_maxgenes)){
      spot_maxgenes = max(col_total_genes)
    }
    if(is.null(spot_maxpct)){
      spot_maxpct = 1
    }
    # Perform SPOT-WISE filter
    spot_mask = (col_total_reads >= spot_minreads & col_total_reads <= spot_maxreads &
                   col_total_genes >= spot_mingenes & col_total_genes <= spot_maxgenes &
                   col_expr_percent >= spot_minpct & col_expr_percent <= spot_maxpct)

    df = df[, spot_mask]

    # Set maximums for GENE-WISE filter
    if(is.null(gene_maxreads)){
      gene_maxreads = max(row_total_reads)
    }
    if(is.null(gene_maxspots)){
      gene_maxspots = max(row_total_spots)
    }
    if(is.null(gene_maxpct)){
      gene_maxpct = 1
    }
    # Perform GENE-WISE filter
    spot_mask = (row_total_reads >= gene_minreads & row_total_reads <= gene_maxreads &
                   row_total_spots >= gene_minspots & row_total_spots <= gene_maxspots &
                   row_expr_percent >= gene_minpct & row_expr_percent <= gene_maxpct)

    df = df[spot_mask, ]

    x@counts[[i]] = tibble::rownames_to_column(df, var='gene')
    x@counts[[i]] = makeSparse(x@counts[[i]])
    x@coords[[i]] = x@coords[[i]][x@coords[[i]]$libname %in% colnames(df), ]
  }

  return(x)
}

