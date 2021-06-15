##
#' @title STList: Creation of STList objects
#' @description Creates an STList object.
#' @details
#' This function takes raw RNA counts and coordinates of those RNA counts from an
#' spatial transcriptomics array. The function can take these inputs in several
#' ways:
#'
#' A. One data frame for counts + one data frame for coordinates.
#' B. One path to a file containing counts + one path to a file containing coordinates.
#' C. One path to a file with file paths to N count files + one path to a file containing N coordinate files.
#' D. One list containing count data frames + one list containing coordinate data frames.
#'
#' The files containing the counts must have gene names in the first column and
#' counts from sampled localities in subsequent columns. The coordinate files must
#' have a first column containing the same identifiers of the sampled locations as in
#' the corresponding count file, and y, x (rows, columns) positions of the spots in
#' the next two columns. The coordinate data must have headers. The function returns an
#' STList object to be used in later analyses with spatialGE. Additionally, the
#' function can also take a comma separated file with clinical or phenotype
#' data from the samples. This clinical file has data for each array (in the same
#' order as the count and coordinate file paths), and the first row has names of
#' clinical variables. The first column contains the user-defined array IDs.
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

  # Check input type.
  input_check <- detect_input(countfiles=countfiles, coordfiles=coordfiles)

  # Create lists to store count matrices from each sample.
  counts_df_list <- list()
  coords_df_list <- list()

  # Store if data frame or list, or a vector.
  if(input_check == 'is_list'){

    for(i in 1:length(countfiles)){
      counts_df_list[[i]] <- tibble::as_tibble(countfiles[[i]])
      coords_df_list[[i]] <- tibble::as_tibble(coordfiles[[i]])
    }

  } else if(input_check == 'is_df'){
    counts_df_list[[1]] <- tibble::as_tibble(countfiles)
    coords_df_list[[1]] <- tibble::as_tibble(coordfiles)
  } else{

    if(input_check == 'is_csv' | input_check == 'is_vector_csv' | input_check == 'is_several_csv'){
      del <- ','
    }

    if(input_check == 'is_tsv' | input_check == 'is_vector_tsv' | input_check == 'is_several_tsv'){
      del <- '\t'
    }

    if(input_check == 'is_csv' | input_check == 'is_tsv'){
      counts_df_list[[1]] <- readr::read_delim(countfiles, delim=del, col_types=readr::cols(),
                                     progress=F)
      coords_df_list[[1]] <- readr::read_delim(coordfiles, delim=del, col_types=readr::cols(),
                                     progress=F)
    }

    if(input_check == 'is_vector_csv' | input_check == 'is_vector_tsv'){

      # Set progress bar, forcing to show with pb$tick(0)
      pb <- progress::progress_bar$new(total=length(countfiles), clear=F)
      pb$tick(0)
      pb$message("Creating STList...")

      for(i in 1:length(countfiles)){

        pb$tick()
        Sys.sleep(1 / length(countfiles))

        # Read filepaths and create.
        counts_df_list[[i]] <- readr::read_delim(countfiles[i], delim=del, col_types=readr::cols(),
                                       progress=F)
        coords_df_list[[i]] <- readr::read_delim(coordfiles[i], delim=del, col_types=readr::cols(),
                                       progress=F)
      }
    }

    if(input_check == 'is_several_csv' | input_check == 'is_several_tsv'){

      count_fpaths <- readLines(countfiles)
      coord_fpaths <- readLines(coordfiles)

      # Set progress bar, forcing to show with pb$tick(0)
      pb <- progress::progress_bar$new(total=length(count_fpaths), clear=F)
      pb$tick(0)
      pb$message("Creating STList...")

      for(i in 1:length(count_fpaths)){

        pb$tick()
        Sys.sleep(1 / length(count_fpaths))

        # Read filepaths and create.
        counts_df_list[[i]] <- readr::read_delim(count_fpaths[i], delim=del, col_types=readr::cols(),
                                                 progress=F)
        coords_df_list[[i]] <- readr::read_delim(coord_fpaths[i], delim=del, col_types=readr::cols(),
                                                 progress=F)
      }

    }

  }

  for(i in 1:length(counts_df_list)){

    # Column names of the count data are simplified using the clean_names() function.
    # NOTE: May need to reconsider use later.
    counts_df_list[[i]] <- janitor::clean_names(counts_df_list[[i]])

    # Clean sample names in coordinates data using the clean_names(), so that they
    # mirror the column names in the count data frame.
    # NOTE: May need to reconsider use later.
    coords_df_list[[i]][, 1] <- janitor::make_clean_names(unlist(coords_df_list[[i]][, 1]))

    # Sort coordinate data according to third column in the coordinate data frame.
    coords_df_list[[i]] <- coords_df_list[[i]][order(coords_df_list[[i]][, 3]), ]

    # Order column names in count data frame according to sorted coordinate data.
    counts_df_list[[i]] <- counts_df_list[[i]][, c(names(counts_df_list[[i]][1]), unlist( coords_df_list[[i]][, 1]))]

    # Test if sample names are the same in both count and coordinate data frames.
    if(length(setdiff(names(counts_df_list[[i]])[-1], unlist(coords_df_list[[i]][, 1]))) != 0){
      stop("There are differences in sample names between the count and coordinate
         data frames.")
    }

    # Test for duplicated gene names in count data.
    gene_names <- counts_df_list[[i]][, 1]
    dup_genes_mask <- duplicated(gene_names)
    if(sum(dup_genes_mask) > 0){
      stop("There are duplicated feature/gene names in the data.")
    }

    # Put column names to coordinate data.
    colnames(coords_df_list[[i]] ) <- c('libname', 'ypos', 'xpos')

  }

  # Test if clinical data is available.
  if(!is.null(clinical)){
    clinical_df <- readr::read_delim(clinical, delim=",", col_types=readr::cols())
    if(nrow(clinical_df) != length(counts_df_list)){
      stop("The number of rows in the clinical data is not the same as the number of spatial arrays.")
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
                    prediction_border=list(),
                    st_clusters=list()
  )

  return(STList_obj)

}
