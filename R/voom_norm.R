##
# This function takes a STList and normalize the count data matricies within it
# in two steps. In the first step, (edgeR) normalization factors are used to scale
# the counts from each library/spot. In the second step, limma-voom normalization
# is applied. The resulting normalized count matrix is stored in the voom_counts
# slot of the STList. The function also calculates gene-wise mean and standard
# deviation from the normalized counts and stores them in the gene_stdev slot.
#
# NOTE: NEED TO IMPLEMENT GENE LENGTH-BASED NORMALIZATION (TPM?)
# NOTE: ADD MEAN-VAR PLOT TO THIS FUNCTION
#
# @param x, a STList with raw count matrices.
# @return x, an updated STList with normalized counts.
#
#
# Load 'tidyverse' for tibble manipulation.
# Load 'edgeR' and 'limma' for normalization.
# require('tidyverse')
# require('edgeR')
# require('limma')
require('magrittr')

voom_norm <- function(x=NULL) {

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Create list to store normalized counts.
  # counts_df_list <- list()

  # Loop through count matrices in STList
  for(i in 1:length(x@counts)){

    # Show progress.
    cat(paste0("Normalizing spatial array #", i, "...\n"))

    # Calculate edgeR normalization factors.
    norm_factors <- edgeR::calcNormFactors(x@counts[[i]][-1], method='TMM',
                                           lib.size=NULL)

    # Create new object to store the size-corrected counts.
    df <- c()

    # Divide each count value by their respective column (sample) normalization
    # factor.
    for(raw_col in names(x@counts[[i]][-1])){
      df <- dplyr::bind_cols(df, tibble::as_tibble(
        x@counts[[i]][raw_col] / norm_factors[raw_col]))
    }

    # Replace gene names.
    # NOTE: This step may not be necessary as gene names are replaced also at the end.
    df <- df %>% tibble::add_column(x@counts[[i]][1], .before=1)

    # Apply voom transformation to count data.
    df_voom <- limma::voom(df[-1], design=NULL,lib.size=colSums(df[-1]),
                           normalize.method='none', plot=F)

    # Estimate gene-wise means and variance and store in object.
    gene_stdevs <- apply(df_voom$E, 1, sd, na.rm=T)
    gene_means <- rowMeans(df_voom$E, na.rm=T)
    gene_stdevs_df <- tibble::tibble(gene_means, gene_stdevs) %>%
      tibble::add_column(df[1], .before=1)
    x@gene_stdev[[i]] <- gene_stdevs_df

    # Put back gene names to matrix and store in object.
    df_voom <- tibble::as_tibble(df_voom$E) %>%
      tibble::add_column(df[1], .before=1)
    x@voom_counts[[i]] <- df_voom
  }

  return(x)

}
