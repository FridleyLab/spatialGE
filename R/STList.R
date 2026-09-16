##
#' @title STlist: Creation of STlist objects for spatial transcriptomics analysis
#' @description Creates an STlist object from one or multiple spatial transcriptomic samples.
#' @details
#' Objects of the S4 class STlist are the starting point of analyses in **`spatialGE`**.
#' The STlist contains data from one or multiple samples (i.e., tissue slices), and
#' results from most `spatialGE`'s functions are stored within the object.
#' \itemize{
#' \item Raw gene counts and spatial coordinates. Gene count data have genes in rows and
#' sampling units (e.g., cells, spots) in columns. Spatial coordinates have
#' sampling units in rows and three columns: sample unit IDs, Y position, and X position.
#' \item Visium outputs from *Space Ranger*. The Visium directory must have the directory
#' structure resulting from `spaceranger count`, with either a count matrix represented in
#' MEX files or a h5 file. The directory should also contain a `spatial` sub-directory,
#' with the spatial coordinates (`tissue_positions_list.csv`), and
#' optionally the high resolution tissue image and scaling factor file `scalefactors_json.json`.
#' \item Xenium outputs from *Xenium Ranger*. The Xenium directory must have the directory
#' structure resulting from the `xeniumranger` pipeline, with either a cell-feature matrix
#' represented in MEX files or a h5 file. The directory should also contain a parquet file,
#' with the spatial coordinates (`cells.parquet`).
#' \item CosMx-SMI outputs. Two files are required to process SMI outputs: The `exprMat` and
#' `metadata` files. Both files must contain the "fov" and "cell_ID" columns. In addition,
#' the `metadata` files must contain the "CenterX_local_px" and "CenterY_local_px" columns.
#'  \item Seurat object (V4). A Seurat V4 object produced via `Seurat::Load10X_Spatial`.
#' }
#' Optionally, the user can input a path to a file containing a table of sample-level
#' metadata (e.g., clinical outcomes, tissue type, age). This sample metadata file
#' should contain sample IDs in the first column partially matching the file names of
#' the count/coordinate file paths or Visium directories. _Note:_ The sample ID of a
#' given sample cannot be a substring of the sample ID of another sample. For example,
#' instead of using "tissue1" and "tissue12", use "tissue01" and "tissue12".
#'
#' The function uses parallelization if run in a Unix system. Windows users
#' will experience longer times depending on the number of samples.
#'
#' @param rnacounts the count data which can be provided in one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing raw gene counts, one file
#' for each spatial sample. The first column contains gene names and subsequent columns
#' contain the expression data for each cell/spot. Duplicate gene names will be
#' modified using `make.unique`. Requires `spotcoords` and `samples`.
#' \item File paths to Visium output directories (one per spatial sample). The directory
#' must follow the structure resulting from `spaceranger count`. The directory contains
#' the `.h5` and `spatial` sub-directory. If no `.h5` file is available, sparse
#' matrices (MEX) from `spaceranger count`. In that case a second sub-directory
#' called `filtered_feature_bc_matrix` should contain contain the `barcodes.tsv.gz`,
#' `features.tsv.gz`, and `matrix.mtx.gz` files. The `spatial` sub-directory minimally
#' contains the coordinates (`tissue_positions_list.csv`), and optionally the high
#' resolution PNG image and accompanying scaling factors (`scalefactors_json.json`).
#' Requires `samples`.
#' #' \item File paths to Xenium output directories (one per spatial sample). The directory
#' must follow the structure resulting from the `xeniumranger` pipeline. The directory
#' contains the `.h5` or sparse matrices (MEX). In that case a second sub-directory
#' called `cell_feature_matrix` should contain contain the `barcodes.tsv.gz`,
#' `features.tsv.gz`, and `matrix.mtx.gz` files. The coordinates must be available
#' in the `cells.parquet`. Requires `samples`.
#' \item The `exprMat` file for each slide of a CosMx-SMI output. The file must contain
#' the "fov" and "cell_ID" columns. The `STlist` function will separate data from each
#' FOV, since analysis in spatialGE is conducted at the FOV level. Requires `samples` and
#' `spotcoords`.
#' \item Seurat object (V4). A Seurat V4 object produced via `Seurat::Load10X_Spatial`.
#' Multiple samples are allowed as long as they are stored as "slices" in the Seurat object.
#' Does not require `samples` as sample names are taken from `names(seurat_obj@images)`
#' \item A named list of data frames with raw gene counts (one data frame per spatial
#' sample). Requires `spotcoords`. Argument `samples` only needed when a file path to
#' sample metadata is provided.
#' }
#' @param spotcoords the cell/spot coordinates. Not required if inputs are Visium or
#' Xenium (spaceranger or xeniumranger outputs).
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing cell/spot coordinates, one
#' for each spatial sample. The files must contain three columns: cell/spot IDs, Y positions, and
#' X positions. The cell/spot IDs must match the column names for each cells/spots (columns) in
#' the gene count files. Requires `samples` and `rnacounts`.
#' \item The `metadata` file for each slide of a CosMx-SMI output. The file must contain
#' the "fov", "cell_ID", "CenterX_local_px", and "CenterY_local_px" columns. The `STlist`
#' function will separate data from each FOV, since analysis in spatialGE is conducted at
#' the FOV level. Requires `samples` and `rnacounts`.
#' \item A named list of data frames with cell/spot coordinates. The list names must
#' match list names of the gene counts list
#' }
#' @param samples the sample names/IDs and (optionally) metadata associated with
#' each spatial sample.
#' The following options are available for `samples`:
#' \itemize{
#' \item A vector with sample names, which will be used to match gene the counts and
#' cell/spot coordinates file paths. A sample name must not match file
#' paths for two different samples. For example, instead of using "tissue1" and
#' "tissue12", use "tissue01" and "tissue12".
#' \item A path to a file containing a table with metadata. This file is a comma- or
#' tab-separated table with one sample per row and sample names/IDs in the first
#' column. Subsequent columns may contain variables associated with each spatial sample
#' }
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems
#' @param verbose logical, whether to print text to console
#' @return an STlist object containing the counts and coordinates, and optionally
#' the sample metadata, which can be used for downstream analysis with `spatialGE`
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export STlist
#'
#' @import Matrix
#' @importFrom magrittr %>%
#
STlist = function(rnacounts=NULL, spotcoords=NULL, samples=NULL, cores=NULL, verbose=TRUE){
  # Check input type.
  input_check = detect_input(rnacounts=rnacounts, spotcoords=spotcoords, samples=samples)
  input_rnacounts = detect_input_rnacounts(rnacounts)
  input_spotcoords = detect_input_spotcoords(spotcoords)
  input_samples = detect_input_samples(samples)

  # Define number of available cores to use.
  cores = resolve_cores(length(rnacounts), cores, coerce='as_integer')
  if(is.na(cores)){
    raise_err(err_code='error0032')
  }

  # Output error if input_check is empty (likely input format not recognized).
  if(rlang::is_empty(input_rnacounts)){
    raise_err(err_code='error0033')
  }

  # Parse sample names
  if(input_rnacounts == 'list_dfs'){
    sample_names = names(rnacounts)
  } else if(input_samples == 'sample_names'){
    sample_names = samples
  } else if(input_rnacounts == 'seurat'){
    sample_names = names(rnacounts@images)
  } else {
    if(input_samples == 'samplesfile_excel'){
      samples_df = readxl::read_excel(samples)
    } else if(input_samples == 'samplesfile_tab'){
      samples_df = readr::read_delim(samples, name_repair='unique', delim='\t', show_col_types=F)
    } else if(input_samples == 'samplesfile_comma'){
      samples_df = readr::read_delim(samples, name_repair='unique', delim=',', show_col_types=F)
    }
    sample_names = samples_df[[1]]
  }

  # Test the validity of sample names
  if(length(sample_names) > 0){
    # Sample names SHOULD NOT begin with a number
    test_number = any(sapply(sample_names, function(i){grepl("^[0-9]", i)}))
    if(test_number){
      raise_err(err_code='error0028')
    }
    # Test that sample names contain only alpha-numerics, spaces, dash, and underscores
    test_chars = any(!sapply(sample_names, function(i){grepl("^[ //-_//A-Za-z0-9]+$", i)}))
    if(test_chars){
      raise_err(err_code='error0028')
    }
    rm(test_number, test_chars) # Clean env
  } else{
    raise_err(err_code='error0004')
  }


  #### CASES BEGIN HERE ########################################################


  # CASE: NAMED LIST OF DATAFRAMES WITH COUNTS AND NAMED LIST OF DATA FRAMES WITH COORDINATES.
  # METADATA INFO OPTIONAL.
  if(input_rnacounts == 'list_dfs' && input_spotcoords == 'list_dfs'){
    if(verbose){
      cat("Found list of dataframes.\n")
    }
    pre_lists = read_list_dfs(rnacounts, spotcoords)
  }


  # CASE: SEURAT OBJECT
  if(input_rnacounts == 'seurat'){
    requireNamespace('SeuratObject')
    if(verbose){
      cat("Found Seurat object.\n")
    }
    pre_lists = read_seurat(rnacounts)
    img_obj = pre_lists[['images']]
    platform = 'visium' # FUTURE DEVELOPMENT: CANT ASSUME IT'S VISIUM... MAYBE SHOULD USE INFO FROM SEURAT OBJECT TO IDENTIFY TECH
  }


  # CASE: SAMPLE FILE PLUS FILE PATH(S) TO COUNT COORDINATE MATRICES OR VISIUM/XENIUM DIRS
  visium = c('visium_filtered_h5', 'visium_raw_h5', 'visium_filtered_mex', 'visium_raw_mex')
  xenium = c('xenium_h5', 'xenium_mex')
  cosmx = c('tab_delim_cosmx', 'comma_delim_cosmx')
  delim = c('tab_delim', 'comma_delim')
  if(input_rnacounts %in% c(visium, xenium, cosmx, delim)){
    # Get list of filepaths
    filepaths = process_sample_names(rnacounts, spotcoords, sample_names, input_rnacounts, input_spotcoords)

    # Check if input is Visium, Xenium, CosMx, or count+coord matrices
    if(input_rnacounts %in% visium){
      if(verbose){
        cat("Found Visium data\n")
      }
      pre_lists = read_visium_outs(filepaths, input_rnacounts, cores=cores, verbose=verbose)
      img_obj = pre_lists[['images']]
      image_scale = pre_lists[['json_scale']]
      platform = 'visium'
    } else if(input_rnacounts %in% xenium){
      if(verbose){
        cat("Found Xenium data\n")
      }
      pre_lists = read_xenium_outs(filepaths, input_rnacounts, cores=cores, verbose=verbose)
      platform = 'xenium'
    } else if(input_check$rna[1] == 'cosmx'){
      if(verbose){
        cat("Found CosMx-SMI data\n")
      }
      pre_lists = read_cosmx_input(filepaths, input_check, cores=cores, verbose=verbose)
      img_obj = pre_lists[['images']]
      platform = 'cosmx'
    } else{
      if(verbose){
        cat("Found matrix data\n")
      }
      pre_lists = read_matrices_fps(filepaths, input_check, cores=cores)
    }
  }


  #### CASES FINISH HERE ########################################################


  # No input provided
  if(is.null(input_check$rna) && is.null(input_check$coords) && is.null(input_check$samples)){
    stop('No input provided. Please refer to documentation.')
  }

  # Process count and coordinate lists before placing within STlist
  if(verbose){
    cat(paste("Matching gene expression and coordinate data...\n"))
  }
  procLists = process_lists(counts_df_list=pre_lists[['counts']], coords_df_list=pre_lists[['coords']])

  # Process metadata if provided or make an empty tibble
  samples_df = tibble::tibble()
  if(input_check$samples[1] == 'samplesfile' || input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex') || input_check$samples[1] == 'samplesfile_matrices'){
    samples_df = process_meta(samples=samples, input_check=input_check, counts_df_list=procLists[['counts']])
  }else{
    samples_df = tibble::tibble(sample_name=names(procLists[['counts']]))
  }
  # Make sure sample IDs are character
  samples_df[[1]] = as.character(samples_df[[1]])

  if(!is.null(input_check$rna[1])){
    if(!(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex'))){
      if(verbose){
        cat(paste("Converting counts to sparse matrices\n"))
      }
      procLists[['counts']] = sge_lapply(procLists[['counts']], function(x){
        makeSparse(x)
      }, cores=cores)
    }
  } else if(!(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex'))){
    if(verbose){
      cat(paste("Converting counts to sparse matrices\n"))
    }
    procLists[['counts']] = sge_lapply(procLists[['counts']], function(x){
      makeSparse(x)
    }, cores=cores)
  }

  # Detect if image from Visium out is available
  if(!exists('img_obj')){
    img_obj = NULL
  }

  # Detect if image scale info is available
  if(!exists('image_scale')){
    image_scale = NULL
  }

  # If no specific platform was found, then make generic
  if(!exists('platform')){
    platform = 'generic'
  }

  # Creates STlist object from processed data
  STlist_obj = new("STlist",
                   counts=procLists[['counts']],
                   spatial_meta=procLists[['coords']],
                   gene_meta=list(),
                   sample_meta=samples_df,
                   tr_counts=list(),
                   #gene_het=list(),
                   gene_krige=list(),
                   #cell_deconv=list(),
                   misc=list(platform=platform,
                             sp_images=img_obj,
                             image_scaling=image_scale,
                             sthet=list())
  )

  if(verbose){
    cat("Completed STlist!\n")
  }

  return(STlist_obj)
}
