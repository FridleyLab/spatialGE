##
# This function takes a data frame with gene names in the first column and
# counts from sampled localities in subsequent columns. It also takes a second
# data frame first column containing the same identifiers of the sampled locations
# in the first data frame, and x, y positions of the samples. The coordinates data
# should not have column names. The function returns an STList object to be used
# in spatial transcriptomics analysis
#
# @param counts, the file path or data frame containing gene names and counts.
# @param coords, the file path ot data frame containing the x, y coordinates.
# @return, the STList object containing the counts and coordinates.
#
#
# Load 'tidyverse' for tibble manipulation.
# Load 'magrittr' in case piping was not enabled by 'tidyverse'
require('tidyverse')
require('magrittr')

STList <- function(counts=NULL, coords=NULL, voom_counts=NULL) {
  # Creates object class, expecting tibbles for counts and coordinates.
  setClass("STList", slots=list(counts="tbl",
                                coords="tbl",
                                voom_counts="tbl"))

  # Test whether counts input is a file path (string) or a data frame.
  # Read file path or assign data frame accordingly.
  if(is.data.frame(counts)){
    counts_df <- as_tibble(counts)
  }else{
    counts_df <- readr::read_delim(counts, delim="\t", col_types=cols())
  }

  # Test whether coords input is a file path (string) or a data frame.
  # Read file path or assign data frame accordingly.
  if(is.data.frame(coords)){
    coords_df <- as_tibble(coords)
  }else{
    coords_df <- readr::read_delim(coords, delim="\t", col_types=cols(),
                                   col_names=F)
  }

  # Column names of the count data are simplified using the clean_names()
  # function.
  # NOTE: May need to reconsider use later.
  counts_df <- counts_df %>% janitor::clean_names()

  # Clean sample names in coordinates data using the clean_names(), so that it
  # mirrors the column names in the count data frame.
  # NOTE: May need to reconsider use later.
  coords_df[, 1] <- janitor::make_clean_names(unlist(coords_df[, 1]))

  # Sort coordinate data according to third column in the coordinate data frame.
  coords_df <- coords_df[order(coords_df[, 3]), ]

  # Order column names in count data frame according to sorted coordinate data.
  counts_df <- counts_df[, c(names(counts_df[, 1]), unlist(coords_df[, 1]))]

  # Test if sample names are the same in both count and coordinate data frames.
  if(length(setdiff(names(counts_df)[-1], unlist(coords_df[, 1]))) != 0){
    stop("There are differences in sample names between the count and coordinate
         data frames.")
  }

  # Test for duplicated gene names in count data.
  gene_names <- counts_df[, 1]
  dup_genes_mask <- duplicated(gene_names)
  if(sum(dup_genes_mask) > 0){
    stop("There are duplicated feature/gene names in the data.")
  }

  # Creates STList object from both count and coordinates data.
  STList_obj <- new("STList", counts=counts_df, coords=coords_df,
                    voom_counts=tibble())

  return(STList_obj)

}


setMethod("show", signature="STList",
          function(object){
            cat("Spatial Transcriptomics List (STList)\n")
            cat(dim(object@counts)[2], "sampled positions.\n")
            cat(dim(object@counts)[1], "features/genes.")
          })
