##
# @title spatial_xcell
# @description Applies xCell to ST data.
# @details
# This function applies xCell to the stored normalized matrices in order to obtain
# cell scores for each of the library/spots. The results are stored as untransformed,
# square-root transformed scores, and p-values. Standard deviations are calculated
# for each cell type.
#
# @param x, a STList with normalized count matrices.
# @return x, an updated STList with xCell scores.
#
#
spatial_xcell <- function(x=NULL){

  require('magrittr')
  require('xCell') # Needs to be 'required' because of databases loaded by packages.

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Test that normalized count matrices are available.
  if(length(x@voom_counts) == 0){
    stop("The STList does not have normalized counts.")
  }

  # Get number of cores available for deconvolution.
  cores <- parallel::detectCores()

  # Create slot to store xCell results.
  x@cell_deconv[['xCell']] <- list()

  for(i in 1:length(x@voom_counts)){

    # Create slot for spatial array within xCell slot.
    x@cell_deconv[['xCell']][[i]] <- list()

    # Extract count data.
    gene_names <- x@voom_counts[[i]]$gene
    df <- as.matrix(x@voom_counts[[i]][, -1])
    rownames(df) <- gene_names

    # Show progress.
    cat(paste0("\nApplying xCell to spatial array #", i, "...\n"))

    # Perform xCell analysis and make cell names as column (no rownames).
    # NOTE1: Applies `janitor` to clean cell names.
    # NOTE 2: Before modifying names, it runs p-value estimation form xCell.

    invisible(capture.output(
      df_xcell <- xCellAnalysis(df, rnaseq=T, parallel.sz=cores)
    ))

    cat(paste0("\nCalculating p-values for spatial array #", i, "...\n"))

    # Extract values for a given cell type.
    # NOTE: This part was planned to get stromal spots...
    xcell_pval <- xCellSignifcanceBetaDist(df_xcell, rnaseq=T)
    colnames(xcell_pval) <- names(x@voom_counts[[i]][-1])

    cell_names <- rownames(df_xcell) %>%
      janitor::make_clean_names(.)
    df_xcell <- tibble::as_tibble(df_xcell) %>%
      tibble::add_column(cell_names, .before=1)

	#NOTE: Make sure .name_repair is not doing something weird here...
    xcell_pval <- tibble::as_tibble(xcell_pval)
    xcell_pval <- xcell_pval %>%
      tibble::add_column(cell_names[1:64], .before=1)
    colnames(xcell_pval) <- colnames(df_xcell)

    # Remove immune, stroma, and microenvironment xCell scores.
    scores_xCell <- grep("immune_score|stroma_score|microenvironment_score",
                          unlist(df_xcell[, 1]))
    df_xcell_NoPurityScores <- df_xcell[-scores_xCell, ]
#    df_xcell_PurityScores <- df_xcell[scores_xCell, ]

    # Store untransformed and transformed scores.
    x@cell_deconv[['xCell']][[i]][['scores']] <- df_xcell_NoPurityScores
    df_xcell_sqrt <- sqrt(df_xcell_NoPurityScores[, -1])
    df_xcell_sqrt <- df_xcell_sqrt %>%
      tibble::add_column(df_xcell_NoPurityScores[, 1], .before=1)
    x@cell_deconv[['xCell']][[i]][['sqrt_scores']] <- df_xcell_sqrt

    # Store p-values.
    x@cell_deconv[['xCell']][[i]][['pvals']] <- xcell_pval

    # Store transformed tumor/stroma scores.
#    df_xcell_purity_sqrt <- sqrt(df_xcell_PurityScores[, -1])
#    df_xcell_purity_sqrt <- df_xcell_purity_sqrt %>%
#      tibble::add_column(df_xcell_PurityScores[, 1], .before=1)
#    x@cell_deconv[[i]]$transf_tumorstroma <- df_xcell_purity_sqrt

    # Calculate cell means and standard deviations, and store in object.
    cell_stdevs <- apply(as.data.frame(
      x@cell_deconv$xCell[[i]]$sqrt_scores[, -1]), 1, sd, na.rm=T)
    # cell_means <- rowMeans(
    #   x@cell_deconv$xCell[[i]]$sqrt_scores[, -1], na.rm=T)
    # cell_stdevs_df <- tibble::tibble(cell_means, cell_stdevs) %>%
    #   tibble::add_column(df_xcell_NoPurityScores[, 1], .before=1)
    cell_stdevs_df <- tibble::tibble(cell_stdevs) %>%
      tibble::add_column(df_xcell_NoPurityScores[, 1], .before=1)
    names(cell_stdevs_df)[1] <- 'cell'
    x@cell_deconv[['xCell']][[i]][['cell_stdev']] <- cell_stdevs_df

  }

  return(x)

}
