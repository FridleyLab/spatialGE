##
# @title spatial_xcell
# @description Applies xCell to ST data.
# @details
# This function applies xCell to the stored transformed matrices in order to obtain
# cell scores for each of the library/spots. The results are stored square-root
# transformed scores, and p-values. Standard deviations are calculated for each cell type.
# The function runs in parallel if unix system available.
#
# @param x, a STList with normalized count matrices.
# @return x, an updated STList with xCell scores.
#
#
spatial_xcell = function(x=NULL){

  require('magrittr')
  require("parallel")
  require('xCell') # Needs to be 'required' because of databases loaded by packages.

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Test that normalized count matrices are available.
  if(length(x@voom_counts) == 0){
    stop("The STList does not have transformed counts.")
  }

  # Define number of available cores to use.
  cores = 1
  if(.Platform$OS.type == 'unix'){
    avail_cores = parallel::detectCores()
    if(avail_cores > (length(x@voom_counts) + 1)){
      cores = (length(x@voom_counts) + 1)
    } else if( (avail_cores <= (length(x@voom_counts) + 1)) && avail_cores > 1){
      cores = avail_cores - 1
    }
  }

  # Estimate xCell cell type scores in parallel if possible.
  xcell_scores = mclapply(seq_along(x@voom_counts), function(i){

    # Extract count data.
    gene_names = x@voom_counts[[i]]$gene
    df = as.matrix(x@voom_counts[[i]][, -1])
    rownames(df) = gene_names

    # Show progress (probably will show when single core).
    cat(paste0("\nApplying xCell to spatial array #", i, "...\n"))

    # Perform xCell analysis.
    #    invisible(capture.output(
    df_xcell = xCellAnalysis(df, rnaseq=T, parallel.sz=1)
    #    ))

    return(df_xcell)

  }, mc.cores=cores)

  # Estimate p-values of xCell scores in parallel if possible.
  xcell_pvals = mclapply(seq_along(xcell_scores), function(i){

    # Show progress (probably will show when single core).
    cat(paste0("\nCalculating p-values for spatial array #", i, "...\n"))

    # Estimate p-values for cell type scores.
    df_pv = xCellSignifcanceBetaDist(xcell_scores[[i]], rnaseq=T)
    colnames(df_pv) = colnames(x@voom_counts[[i]][-1])

    return(df_pv)

  }, mc.cores=cores)

  # Create list to be stored in STList.
  x@cell_deconv[['xCell']] = list()

  # Applies `janitor` to clean cell type names of matrix of first ST array (they will
  # be applied to the others, assuming xCell produces the same cell types.)
  cell_names = rownames(xcell_scores[[1]]) %>%
    janitor::make_clean_names(.)

  # Store results in list.
  for(i in 1:length(xcell_scores)){

    # Remove immune, stroma, and microenvironment xCell scores.
    purity_names = grep("immune_score|stroma_score|microenvironment_score", cell_names)
    xcell_scores[[i]] = xcell_scores[[i]][-purity_names, ]

    # Replace cell type names.
    xcell_scores[[i]] = sqrt(xcell_scores[[i]]) %>%
      tibble::as_tibble(xcell_scores[[i]]) %>%
      tibble::add_column(cell_names=cell_names[-purity_names], .before=1)

    # NOTE: Make sure .name_repair is not doing something weird here...
    xcell_pvals[[i]] = tibble::as_tibble(xcell_pvals[[i]]) %>%
      tibble::add_column(cell_names=cell_names[-purity_names], .before=1)

    x@cell_deconv[['xCell']][[i]] = list()
    x@cell_deconv[['xCell']][[i]][['sqrt_scores']] = xcell_scores[[i]]
    x@cell_deconv[['xCell']][[i]][['pvals']] = xcell_pvals[[i]]

  }

  # Calculate cell means and standard deviations, and store in object in parallel (if possible).
  cell_stdevs = mclapply(seq_along(xcell_scores), function(i){

    df_stdev = apply(as.data.frame(xcell_scores[[i]][, -1]), 1, sd, na.rm=T)

    df_stdev = tibble::as_tibble_col(df_stdev, column_name='cell_stdevs') %>%
      tibble::add_column(cell=cell_names[-purity_names], .before=1)

  }, mc.cores=cores)

  # Put standard deviations in corresponding slot.
  for(i in 1:length(cell_stdevs)){
    x@cell_deconv[['xCell']][[i]][['cell_stdev']] = cell_stdevs[[i]]
  }

  return(x)

}
