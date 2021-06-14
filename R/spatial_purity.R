##
# @title spatial_purity
# @description Applies ESTIMATE to ST data.
# @details
# This function applies ESTIMATED to the stored normalized matrices in order to
# obtain purity scores for each of the library/spots.
#
# @param x, a STList with normalized count matrices.
# @return x, an updated STList with ESTIMATE purity scores.
#
#
spatial_purity <- function(x=NULL){

  require('magrittr')
  require('estimate') # Needs to be 'required' because of databases loaded by package.

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Test that normalized count matrices are available.
  if(length(x@voom_counts) == 0){
    stop("The STList does not have normalized counts.")
  }

  # Create slot to store ESTIMATE results.
  x@cell_deconv[['ESTIMATE']] <- list()

  for(i in 1:length(x@voom_counts)){

    # Create slot for spatial array within ESTIMATE slot.
    x@cell_deconv[['ESTIMATE']][[i]] <- list()

    # Extract count data.
    gene_names <- x@voom_counts[[i]]$gene
    df <- as.data.frame(x@voom_counts[[i]][, -1])
    rownames(df) <- gene_names

    # Write temporary file to store expression data. Required by ESTIMATE.
    tmp_expr <- tempfile(fileext = ".txt", pattern='voomExpr')
    write.table(df, tmp_expr, quote=F, sep="\t")

    # Show progress.
    cat(paste0("Estimating 'purity' scores for spatial array #", i, "...\n"))

    # Write temporary file to store filtered data set. Required by ESTIMATE.
    tmp_filterexpr <- tempfile(fileext = ".gct", pattern='estimateCommonGenes_')
    estimate::filterCommonGenes(input.f=tmp_expr, output.f=tmp_filterexpr, id="GeneSymbol")

    # Write temporary file to store filtered data set. Required by ESTIMATE.
    tmp_purout <- tempfile(fileext = ".gct", pattern='purityOut_')
    estimateScore(tmp_filterexpr, tmp_purout, platform="illumina")

    # Read results from temporary file.
    df_purity <- readr::read_delim(tmp_purout, delim="\t", skip=2,
                                   col_types=readr::cols())

    # Store scores.
    # NOTE: Scores are not transformed given that ESTIMATE scores seem to
    # follow a normal distribution.
    x@cell_deconv[['ESTIMATE']][[i]]$estimate_purity <- df_purity[, -2]

  }

  return(x)

}
