##
# This function takes a data frame with gene names in the first column and
# counts from sampled localities in subsequent columns. It also takes a second
# data frame first column containing the same identifiers of the sampled locations
# in the first data frame, and x, y positions of the samples. The coordinate data
# should not have column names. The function returns an STList object to be used
# in spatial transcriptomics analysis
#
# @param counts, the file path or data frame containing gene names and counts.
# @param coords, the file path or data frame containing the x, y coordinates.
# @param clinical, the file path or data frame containing clinical data associated
# with the count data.
# NOTE: This clinical data frame needs more thinking! What to do with it? format?
# @return, the STList object containing the counts and coordinates.
#
#
# Load 'tidyverse' for tibble manipulation.
# Load 'magrittr' in case piping was not enabled by 'tidyverse'
# require('tidyverse')
require('magrittr')

STList <- function(countfiles=NULL, coordfiles=NULL, clinical=NULL) {

  # Creates object class, expecting tibbles for counts and coordinates.
  setClass("STList", slots=list(counts="list",
                                coords="list",
                                clinical="tbl",
                                voom_counts="tbl",
                                gene_stdev="tbl",
                                gene_het="list",
                                gene_krige="list",
                                cell_deconv="list",
                                cell_stdev="tbl",
                                cell_het="list",
                                cell_krige="list"
  ),
  )

  # Test if counts and coords are NULL. If so, print error.
  if(is.null(countfiles) | is.null(coordfiles)){
    stop("Must have filepaths for BOTH counts and coordinates data.")
  }

  # Get count and coord filepaths and test there is an equal number of both.
  count_fpaths <- scan(countfiles, what=character(), quiet=T)
  coord_fpaths <- scan(coordfiles, what=character(), quiet=T)
  if(length(count_fpaths) != length(coord_fpaths)){
    stop("The number of count files is different to the number of coordinate files.")
  }

  # Create lists to store count matrices from each sample.
  counts_df_list <- list()
  coords_df_list <- list()

  # Set progress bar, forcing to show with pb$tick(0)
  pb <- progress::progress_bar$new(total=length(count_fpaths), clear=F)
  pb$tick(0)
  pb$message("Creating STList...")

  # Process each of the file sets.
  for(i in 1:length(count_fpaths)){

    # Update progress bar.
    pb$tick()
    Sys.sleep(1 / length(count_fpaths))

    counts <- count_fpaths[i]
    coords <- coord_fpaths[i]

    # Test that a pair of count and coord files exist.
    if(!(file.exists(counts) & file.exists(coords))){
      stop("Either one of the count files or coordinate files do not exists.")
    }

    # Read filepaths and create.
    counts_df <- readr::read_delim(counts, delim="\t", col_types=readr::cols(),
                                   progress=F)
    coords_df <- readr::read_delim(coords, delim="\t", col_types=readr::cols(),
                                   progress=F, col_names=F)

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

    counts_df_list[[i]] <- counts_df
    coords_df_list[[i]] <- coords_df

  }

  # Creates STList object from both count and coordinates data.
  STList_obj <- new("STList",
                    counts=counts_df_list,
                    coords=coords_df_list,
                    clinical=tibble::tibble(),
                    voom_counts=tibble::tibble(),
                    gene_stdev=tibble::tibble(),
                    gene_het=list(),
                    gene_krige=list(),
                    cell_deconv=list(),
                    cell_stdev=tibble::tibble(),
                    cell_het=list(),
                    cell_krige=list()
  )

  return(STList_obj)

}
