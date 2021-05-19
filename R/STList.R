##
#' @title STList: Creation of STList objects
#' @description Creates an STList object.
#' @details
#' This function takes two files containing file paths of count matrices and
#' coordinates, one per line. The files containing the counts must have gene names
#' in the first column and counts from sampled localities in subsequent columns.
#' The coordinate files must have a first column containing the same identifiers
#' of the sampled locations as in in the corresponding count file, and x, y positions
#' of the samples in the next two columns. The coordinate data must not have column
#' names. The function returns an STList object to be used in later analyses with
#' spatialGE. Additionally, the function can also take a comma separated file with
#' clinical or phenotype data from the samples. This clinical file has data for each
#' array (in the same order as the count and coordinate file paths), and the first
#' row has names of clinical variables. The first column contains the user-defined
#' array IDs.
#'
#' @param countfiles, the path of a file containing file paths (one per line) of the
#' count data files. Count files must be tab-separated.
#' @param coords, the path of a file containing file paths (one per line) of the
#' coordinate data files. Coordinate files must be tab-separated.
#' @param clinical, the file path of a file containing a table with clinical/phenotype
#' variables associated with the spatial arrays. The order of rows must match that of
#' the count and coordinate file paths. the first column is the user-defined ID of the
#' arrays.
#' @return x, the STList object containing the counts and coordinates, and optionally
#' the clinical data.
#' @export
#
#
STList <- function(countfiles=NULL, coordfiles=NULL, clinical=NULL) {

  # Test if counts and coords are NULL. If so, print error.
  if(is.null(countfiles) | is.null(coordfiles)){
    stop("Must have filepaths for BOTH counts and coordinates data.")
  }

  # Get count and coord filepaths and test there is an equal number of both.
  count_fpaths <- readLines(countfiles) #, what=character(), quiet=T)
  coord_fpaths <- readLines(coordfiles) #, what=character(), quiet=T)
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

    # Get file paths for a specific array.
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

    # Column names of the count data are simplified using the clean_names() function.
    # NOTE: May need to reconsider use later.
    counts_df <- janitor::clean_names(counts_df)

    # Clean sample names in coordinates data using the clean_names(), so that they
    # mirror the column names in the count data frame.
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

    # Store count and coordinate matrices in lists.
    counts_df_list[[i]] <- counts_df
    coords_df_list[[i]] <- coords_df

  }

  # Test if clinical data is available.
  if(!is.null(clinical)){
    clinical_df <- readr::read_delim(clinical, delim=",", col_types=readr::cols())
    if(nrow(clinical_df) != length(count_fpaths)){
      stop("The number of rows in the clinical data is not the same as the
             number of spatial arrays.")
    }
  }else{
    clinical_df <- tibble::tibble()
  }

  # Creates STList object from both count and coordinates data.
  STList_obj <- new("STList",
                    counts=counts_df_list,
                    coords=coords_df_list,
                    clinical=clinical_df,
                    voom_counts=list(),
                    gene_stdev=list(),
                    gene_het=list(),
                    gene_krige=list(),
                    cell_deconv=list(),
                    cell_het=list(),
                    cell_krige=list(),
                    prediction_border=list()
  )

  return(STList_obj)

}
