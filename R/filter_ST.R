##
# This function provides a method to remove spots based on the number of reads,
# number of genes, and percentage of reads coming from mtDNA.
#
#
filter_ST = function(x=NULL,
                     minreads=NULL, maxreads=NULL,
                     mingenes=NULL, maxgenes=NULL,
                     minMT=NULL, maxMT=NULL,
                     who=NULL){

  if(is.null(x)){
    stop('Please, specify an STList to filter.')
  }
  if(is.null(minreads)){
    minreads = 0
  }
  if(is.null(mingenes)){
    mingenes = 0
  }
  if(is.null(minMT)){
    minMT = 0
  }
  if(is.null(who)){
    who = 1:length(x@counts)
  }

  for(i in who){

    df = x@counts[[i]]
    total_reads = colSums(x@counts[[i]][, -1])
    nonzero_df = x@counts[[i]][, -1] != 0
    total_genes = colSums(nonzero_df)
    total_mtct = colSums(x@counts[[i]][grep("^MT-", x@counts[[i]][[1]]), -1])
    percent_mt = total_mtct/total_reads

    if(is.null(maxreads)){
      maxreads = max(total_reads)
    }
    if(is.null(maxgenes)){
      maxgenes = max(total_genes)
    }
    if(is.null(maxMT)){
      maxMT = max(percent_mt)
    }

    spot_mask = (total_reads >= minreads & total_reads <= maxreads &
                   total_genes >= mingenes & total_genes <= maxgenes &
                   percent_mt >= minMT & percent_mt <= maxMT)

    df = df[, c(TRUE, spot_mask)]

    x@counts[[i]] = df
    x@coords[[i]] = x@coords[[i]][(x@coords[[i]]$libname %in% colnames(x@counts[[i]])[-1]), ]
  }

  return(x)

}

