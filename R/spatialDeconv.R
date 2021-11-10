##
#' @title spatialDeconv: Deconvolution of ST data
#' @description Applies gene expression deconvolution methods to data from spatial
#' arrays. Produces tumor/stroma classes via model-based clustering of ESTIMATE purity scores.
#' @details
#' This function applies deconvolution methods (xCell and ESTIMATE, for now) to the stored
#' transformed count matrices in order to obtain cell scores for each of the library/spots.
#' The results are stored within the STList. The method ESTIMATE is automatically
#' applied regardless of the selected method.
#'
#' @param x an STList with normalized count matrices.
#' @param method the name of the method to deconvolute data. As of now, only xCell is supported.
#' If 'none', only ESTIMATE is applied.
#' @return x an updated STList with deconvolution scores.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # melanoma <- spatialDeconv(melanoma, method='xcell')
#'
#' @export
#
#
spatialDeconv = function(x=NULL, method='xcell'){

  method=tolower(method)

  # If no ESTIMATE scores available, then generate them
  if(is.null(x@cell_deconv[['ESTIMATE']])){
    cat('This will take some time...\n')
    # Get ESTIMATE scores
    estimate_scores <- spatial_purity(x)
    # Store ESTIMATE scores in STList
    x@cell_deconv[['ESTIMATE']] = list()
    for(i in 1:length(estimate_scores)){
      x@cell_deconv[['ESTIMATE']][[i]] = list()
      x@cell_deconv[['ESTIMATE']][[i]][['estimate_purity']] = estimate_scores[[i]]
    }

    #Generate clusters of purity scores
    purity_clusters = cluster_purity(x)
    # Store tumor/stroma classes in STList
    for(i in 1:length(purity_clusters)){
      x@cell_deconv[['ESTIMATE']][[i]][['purity_clusters']] = purity_clusters[[i]]
    }
  }

  if(method == 'none' | method == 'estimate'){
    cat('Only ESTIMATE was applied to the spatial arrays.')
  }

  if(method == 'xcell'){
    require('xCell') # Needs to be 'required' because of databases loaded by packages.
    x <- spatial_xcell(x)
  }

  return(x)

}


# Helpers ----------------------------------------------------------------------

##
# @title spatial_purity
# @description Applies ESTIMATE to ST data.
# @details
# This function applies ESTIMATED to the stored transformed matrices in order to
# obtain purity scores for each of the library/spots. The function runs on parallel
# if unix system available.
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

  # Test that transformed count matrices are available.
  if(length(x@tr_counts) == 0){
    stop("The STList does not have transformed counts.")
  }

  # Define number of available cores to use.
  cores = count_cores(length(x@tr_counts))

  # Get ESTIMATE scores in parallel if possible.
  estimate_scores = parallel::mclapply(seq_along(x@tr_counts), function(i){
    # Extract count data.
    df = as.matrix(x@tr_counts[[i]][, -1])
    rownames(df) = x@tr_counts[[i]]$gene

    # Write temporary file to store expression data. Required by ESTIMATE.
    tmp_expr <- tempfile(fileext = ".txt", pattern='trExpr_')
    write.table(df, tmp_expr, quote=F, sep="\t")

    # If single core, will probably show progress.
    system(sprintf('echo "%s"',
                   crayon::yellow(paste0("Estimating 'purity' scores for sample ", names(x@tr_counts)[i], "...\n"))))
    #cat(paste0("Estimating 'purity' scores for sample ", names(x@tr_counts)[i], "...\n"))

    # Write temporary file to store filtered data set. Required by ESTIMATE.
    tmp_filterexpr <- tempfile(fileext = ".gct", pattern='estimateCommonGenes_')
    estimate::filterCommonGenes(input.f=tmp_expr, output.f=tmp_filterexpr, id="GeneSymbol")

    # Write temporary file to store filtered data set. Required by ESTIMATE.
    tmp_purout <- tempfile(fileext = ".gct", pattern='purityOut_')
    estimate::estimateScore(tmp_filterexpr, tmp_purout, platform="illumina")

    # Read results from temporary file.
    df_purity <- readr::read_delim(tmp_purout, delim="\t", skip=2, col_types=readr::cols())

    # NOTE: Scores are not transformed given that ESTIMATE scores seem to follow a normal distribution.
    df_purity <- df_purity[3, -2]

    return(df_purity)

  }, mc.cores=cores, mc.preschedule=F)

  return(estimate_scores)
}

##
# @title cluster_purity
# @description Perform model-based clustering to define stromal/tumor clusters from
# ESTIMATE purity scores.
# @details
# Takes ESTIMATE tumor purity scores and performs model-based clustering. The function
# only evaluates k=2, to separate between likely tumor and stromal spots.
#
# @param x, an STList with ESTIMATE scores.
# @return x, an STList with tumor/stroma classifications.
#
#
cluster_purity <- function(x=NULL) {

  suppressMessages(
    require('mclust') #Mclust can't call mclust if the library is not loaded.
  )

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Detect usable cores
  cores = count_cores(length(x@cell_deconv[['ESTIMATE']]))

  # Iterate over spatial arrays.
  purity_clusters =  parallel::mclapply(seq_along(x@cell_deconv[['ESTIMATE']]), function(i){

    # Extract purity scores.
    estimate_df <- x@cell_deconv$ESTIMATE[[i]]$estimate_purity

    # Transpose values to cluster.
    estimate_mx <- t(as.matrix(estimate_df[1, -1]))

    # Perform clustering.
    BIC <- mclust::mclustBIC(estimate_mx, G=c(2), verbose=F)
    clust_mod <- mclust::Mclust(estimate_mx, x=BIC, G=c(2), verbose=F)

    # Store tumor/stroma classifications.
    clusters <- tibble::as_tibble_col(clust_mod$classification, column_name='cluster')
    clusters$cluster <- gsub('1', 'tumor', clusters$cluster)
    clusters$cluster <- gsub('2', 'stroma', clusters$cluster)
    clusters <- tibble::add_column(clusters, libname=names(clust_mod$classification), .before=1)

    return(clusters)
  }, mc.cores=cores, mc.preschedule=F)

  return(purity_clusters)
}


