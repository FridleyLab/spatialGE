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
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
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
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores) && !is.null(rnacounts)){
    cores = count_cores(length(rnacounts))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      raise_err(err_code='error0032')
    }
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
  #}else if(input_check$samples[1] == 'samplesfile_geomx'){
  #  samples_df = process_meta_geomx(samples, input_check, procLists[['counts']], gmx_slide_col, gmx_meta_cols)
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
      procLists[['counts']] = parallel::mclapply(procLists[['counts']], function(x){
        makeSparse(x)
      })
    }
  } else if(!(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex'))){
    if(verbose){
      cat(paste("Converting counts to sparse matrices\n"))
    }
    procLists[['counts']] = parallel::mclapply(procLists[['counts']], function(x){
      makeSparse(x)
    })
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


# Helpers ----------------------------------------------------------------------

##
# makeSparse: takes a dataframe input and makes it sparse saving space
# wrapped from Matrix package
# @param dataframe, a dataframe or tibble with gene names in the first column and
# different spots in the following columns
# @return a sparsed data matrix
#
makeSparse = function(dataframe){

  # To prevent NOTES in R CMD check
  . = NULL

  if(any(class(dataframe) == 'dgCMatrix')){
    numdat = dataframe
  } else if(!is.matrix(dataframe)){
    # Force data frame if tibble (tibbles do not support rownames)
    if(tibble::is_tibble(dataframe)){
      dataframe = as.data.frame(dataframe)
    }
    genecol = colnames(dataframe)[1]

    numdat = dataframe
    rownames(numdat) = dataframe[[genecol]]
    numdat = numdat[, -which(colnames(numdat) == genecol)]
    numdat = as.matrix(numdat)
    numdat = as(numdat, "sparseMatrix")
  } else {
    numdat = dataframe %>% as(., "sparseMatrix")
  }
  return(numdat)
}

##
# expandSparse: takes a sparsed matrix and returns a dataframe
# @param sparsedMatrix, a sparse matrix made from the Matrix
# package
# @return a data frame
#
expandSparse = function(sparsedMatrix){
  NonSparse = data.frame(as.matrix(sparsedMatrix), check.names=F)
  return(NonSparse)
}

##
# read_list_dfs: Takes two named lists of counts and coordinates and performs sorting
# to get count and coordinate data frames in the same order
# @param rnacounts, a named list with counts
# @param spotcoords, a named list with coordinates
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_list_dfs = function(rnacounts, spotcoords){
  # Check that RNA and coord lists have the same number of elements.
  if(length(rnacounts) != length(spotcoords)){
    stop('The RNA count and coordinate do not have the same number of elements.')
  }
  # Test that names within the lists are the same.
  if(length(setdiff(names(rnacounts), names(spotcoords))) != 0){
    stop('The RNA count and coordinate lists do not have the same names.')
  }
  # Sort lists' names to order both.
  sorted_names = sort(names(rnacounts))

  # Process duplicate gene names
  # NOTE: Duplicate genes coming from matrices are modified with `make.unique`, but
  # are not comparable across matrices because unique identifiers for the probes are not used.
  for(i in 1:length(rnacounts)){
    if(any(duplicated(rnacounts[[i]][[1]]))){
      rnacounts[[i]][[1]] = make.unique(rnacounts[[i]][[1]])
    }
  }

  # Create list to be returned
  return_lists = list()
  return_lists[['counts']] = rnacounts[sorted_names]
  return_lists[['coords']] = spotcoords[sorted_names]
  return(return_lists)
}

##
# read_seurat: Takes a Seurat object and converts it to a STlist.
# @param rnacounts, a seurat object with spatial slot
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_seurat = function(rnacounts){
  # Create list to be returned
  return_lists = list()
  return_lists[['counts']] = list()
  return_lists[['coords']] = list()
  # Find assays within Seurat object
  slices = names(rnacounts@images)
  counts_mtx = as.data.frame(rnacounts@assays$Spatial@counts)
  # Extract counts and coordinates from each assay
  for(i in slices){
    spots = rownames(rnacounts@images[[i]]@coordinates)
    return_lists[['coords']][[i]] = rnacounts@images[[i]]@coordinates[, c('imagerow', 'imagecol')]
    return_lists[['coords']][[i]] = tibble::rownames_to_column(return_lists[['coords']][[i]], var='libname')
    return_lists[['coords']][[i]] = tibble::as_tibble(return_lists[['coords']][[i]])

    return_lists[['counts']][[i]] = counts_mtx[, colnames(counts_mtx) %in% spots]
    return_lists[['counts']][[i]] = tibble::rownames_to_column(return_lists[['counts']][[i]], var='gene')
    return_lists[['counts']][[i]] = tibble::as_tibble(return_lists[['counts']][[i]])

    return_lists[['images']][[i]] = rnacounts@images[[i]]@image
  }
  return(return_lists)
}

##
# read_matrices_fps: Takes a list with file paths to count and coordinate matrices
# and sample names and returns a list with count and coordinates data frames per sample
# @param filepaths a list with file paths to Visium outputs and sample IDs
# @param input_check The result of detect_input
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_matrices_fps = function(filepaths, input_check, cores=NULL){
  # Get delimiter of RNA counts and coordinates.
  # NOTE: This only test the first matrix.
  # NOTE2: This is being commented out because presumably, the detect_input function
  # already provides de delimiter to use
  if(!is.null(input_check$rna) & !is.null(input_check$coords)){
    delrna = input_check$rna[2]
    delcoords = input_check$coords[2]
    if(!(delrna %in% c(',', '\t')) | !(delcoords %in% c(',', '\t'))){
      stop('Counts or coordinate files are not comma or tab- delimited') # Print error if files were not detected to be comma- or tab-delimited
    }
  } else{ # If delimiter is not available (for example because only sample file provided), find delimiter
    count_testline = readLines(filepaths[['count_found']][1], n=1)
    if(grepl('\t', count_testline)){
      delrna = '\t'
    } else if(grepl(',', count_testline)){
      delrna = ','
    } else{
      stop('RNA counts are not comma or tab- delimited')
    }
    coord_testline = readLines(filepaths[['coord_found']][1], n=1)
    if(grepl('\t', coord_testline)){
      delcoords = '\t'
    } else if(grepl(',', coord_testline)){
      delcoords = ','
    } else{
      stop('Sport coordinates are not comma or tab- delimited')
    }
  }

  # Use parallelization to read count data if possible.
  counts_df_list = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    # Read filepaths.
    #counts_df = readr::read_delim(filepaths[['count_found']][i], delim=delrna, col_types=readr::cols(), progress=F)
    counts_df = data.table::fread(filepaths[['count_found']][i], sep=delrna)
    counts_df = as.data.frame(counts_df)
    return(counts_df)
  }, mc.cores=cores, mc.preschedule=F)
  # Name list elements.
  names(counts_df_list) = filepaths[['sampleids']]

  # Use parallelization to read coordinate data if possible.
  coords_df_list = parallel::mclapply(seq_along(filepaths[['coord_found']]), function(i){
    # Read filepaths.
    coords_df = readr::read_delim(filepaths[['coord_found']][i], delim=delcoords, col_types=readr::cols(), progress=F)
    return(coords_df)
  }, mc.cores=cores, mc.preschedule=F)
  # Name list elements.
  names(coords_df_list) = filepaths[['sampleids']]

  rm(delrna, delcoords) # Clean environment

  # Process duplicate gene names
  # NOTE: Duplicate genes coming from matrices are modified with `make.unique`, but
  # are not comparable across matrices because unique identifiers for the probes are not used.
  for(i in 1:length(counts_df_list)){
    if(any(duplicated(counts_df_list[[i]][[1]]))){
      counts_df_list[[i]][[1]] = make.unique(counts_df_list[[i]][[1]])
    }
  }

  return_lists = list()
  return_lists[['counts']] = counts_df_list
  return_lists[['coords']] = coords_df_list
  return(return_lists)
}

##
# read_visium_outs: Takes a list with Visium output file paths and sample names and
# returns a list with count and coordinates data frames per sample
# @param filepaths a list with file paths to Visium outputs and sample IDs
# @param input_rnacounts The result of detect_input_rnacounts
# @return a list with two lists within (one with counts, one with coordinates)
#
read_visium_outs = function(filepaths, input_rnacounts, cores=NULL, verbose=TRUE){

  # To prevent NOTES in R CMD check
  . = NULL

  # Find necessary files from visium input
  missingSamples = 0
  fp_list = list()
  for(i in 1:length(filepaths[['count_found']])){
    fp_list[[i]] = list()
    # Get all system paths within output folder (h5/MEX, coordinates, image)
    temp_fps = list.files(filepaths[['count_found']][i], recursive=T, include.dirs=T, full.names=T)

    vimage = grep('spatial\\/', temp_fps, value=T) %>%
      grep('tissue_hires_image.png|tissue_lowres_image.png', ., value=T)
    vjson = grep('spatial\\/', temp_fps, value=T) %>%
      grep('scalefactors_json.json$', ., value=T)
    vcoords = grep('spatial\\/', temp_fps, value=T) %>%
      grep('tissue_positions_list.csv|tissue_positions.csv', ., value=T)
    # Filter out 'SPATIAL_RNA_COUNTER' folders (intermediate files from Space Ranger?).
    vcoords = vcoords[!grepl('SPATIAL_RNA_COUNTER', vcoords)]

    if(input_rnacounts[1] %in% c('visium_filtered_mex', 'visium_raw_mex')){
      if(input_rnacounts[1] == 'visium_filtered_mex'){
        # NOTE: For MEX, assume files are inside a directory named 'filtered_feature_bc_matrix'
        vfeatures = grep('filtered_feature_bc_matrix\\/features.tsv.gz',  temp_fps, value=T)
        vbarcodes = grep('filtered_feature_bc_matrix\\/barcodes.tsv.gz', temp_fps, value=T)
        vcounts = grep('filtered_feature_bc_matrix\\/matrix.mtx.gz', temp_fps, value=T)
      } else if(input_rnacounts[1] == 'visium_raw_mex'){
        vfeatures = grep('raw_feature_bc_matrix\\/features.tsv.gz',  temp_fps, value=T)
        vbarcodes = grep('raw_feature_bc_matrix\\/barcodes.tsv.gz', temp_fps, value=T)
        vcounts = grep('raw_feature_bc_matrix\\/matrix.mtx.gz', temp_fps, value=T)
      }

      # Test that all files have been found.
      needed_mex_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('csv', vcoords))
      if(any(needed_mex_test)){
        raise_err(err_code='error0035')
      }

      fp_list[[i]]$features = vfeatures
      fp_list[[i]]$barcodes = vbarcodes
      fp_list[[i]]$counts = vcounts

      if(rlang::is_empty(vfeatures)) message(paste("Features for", filepaths$sampleids[i], "not able to be found..."))
      if(rlang::is_empty(vbarcodes)) message(paste("Barcodes for", filepaths$sampleids[i], "not able to be found..."))
      if(rlang::is_empty(vcounts)) message(paste("Counts for", filepaths$sampleids[i], "not able to be found..."))

      if(rlang::is_empty(vfeatures) | rlang::is_empty(vbarcodes) | rlang::is_empty(vcounts) | rlang::is_empty(vcoords)){
        fp_list[[i]] = list()
        missingSamples  = missingSamples + 1
      }

      rm(vfeatures, vbarcodes, vcounts) # Clean environment
    } else if(input_rnacounts[1] %in% c('visium_filtered_h5', 'visium_raw_h5')){
      if(input_rnacounts[1] == 'visium_filtered_h5'){
        h5counts = grep('filtered_feature_bc_matrix.h5', temp_fps, value=T)
      } else if(input_rnacounts[1] == 'visium_raw_h5'){
        h5counts = grep('raw_feature_bc_matrix.h5', temp_fps, value=T)
      }

      # Test that all files have been found.
      needed_h5_test = c(!grepl('\\.h5$', h5counts), !grepl('csv', vcoords))
      if(any(needed_h5_test)){
        raise_err(err_code='error0035')
      }

      fp_list[[i]]$counts = h5counts

      if(rlang::is_empty(vcoords)) message(paste("Coordinates for", filepaths$sampleids[i], "not able to be found...\n"))
      if(rlang::is_empty(h5counts) | rlang::is_empty(vcoords)){
        fp_list[[i]] = list()
        missingSamples  = missingSamples + 1
      }
      rm(h5counts) # Clean environment
    }

    fp_list[[i]]$coords = vcoords
    fp_list[[i]]$image = vimage
    fp_list[[i]]$json = vjson
    fp_list[[i]]$runname = filepaths[['sampleids']][i]

    rm(temp_fps, vcoords, vimage) # Clean environment
  }

  if(verbose){
    cat(paste("\tFound", length(filepaths$sampleids)-missingSamples, "Visium samples\n"))
  }

  # Define number of available cores to use.
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(filepaths[['count_found']]))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      raise_err(err_code='error0036')
    }
  }

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    if(length(fp_list[[i]]$counts) == 0){
      return(list())
    }
    if(verbose){
      system(sprintf('echo "%s"', paste0("\t\tProcessing Sample ", i, "....")))
    }

    # Process Visium outputs
    if(input_rnacounts %in% c('visium_filtered_mex', 'visium_raw_mex')){
      visium_processed = import_visium(features_fp=fp_list[[i]][['features']],
                                       barcodes_fp=fp_list[[i]][['barcodes']],
                                       counts_fp=fp_list[[i]][['counts']],
                                       coords_fp=fp_list[[i]][['coords']])
    } else{
      visium_processed = import_visium_h5(counts_fp=fp_list[[i]][['counts']],
                                          coords_fp=fp_list[[i]][['coords']])
    }

    if(verbose){
      system(sprintf('echo "%s"', paste0("\t\tFinished data read Sample ", i)))
    }

    return(visium_processed)
  }, mc.cores=cores, mc.preschedule=F)
  if(verbose){
    cat(paste("\tData read completed\n"))
  }

  # Organize the paralellized output into corresponding lists.
  return_lists = list()
  return_lists[['counts']] = list()
  return_lists[['coords']] = list()
  return_lists[['images']] = list()
  return_lists[['json_scale']] = list()
  for(i in 1:length(output_temp)){
    return_lists[['counts']][[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
    return_lists[['coords']][[fp_list[[i]]$runname]] = output_temp[[i]]$coords
    return_lists[['images']][[fp_list[[i]]$runname]] = NULL
    return_lists[['json_scale']][[fp_list[[i]]$runname]] = NULL
    # Test if image file is present
    if(!rlang::is_empty(fp_list[[i]]$image)){
      # Try reading image
      image_read = tryCatch({png::readPNG(fp_list[[i]]$image[1])},
                            error=function(e){return(e)})
      # Try reading image again using the second element if available
      if(any(class(image_read) == 'simpleError') & length(fp_list[[i]][['image']]) > 1){
        image_read = tryCatch({png::readPNG(fp_list[[i]]$image[2])},
                              error=function(e){return(e)})
      }
      # If no image reading was possible, do not output image
      if(any(class(image_read) == 'simpleError')){
        image_read = NULL
        warning('Tissue image could not be read. Is the image a PNG derived from Space Ranger?')
      }
      return_lists[['images']][[fp_list[[i]]$runname]] = image_read
    }
    # Test if JSON scaling factor file is present
    if(!rlang::is_empty(fp_list[[i]]$json)){
      # image_read = tryCatch({jsonlite::read_json(fp_list[[i]]$json[1])},
      #                       error=function(e){png::readPNG(fp_list[[i]]$image[2])})
      json_read = jsonlite::read_json(fp_list[[i]]$json[1])
      return_lists[['json_scale']][[fp_list[[i]]$runname]] = json_read
    } else{
      return_lists[['json_scale']][[fp_list[[i]]$runname]] = 'no_scaling_available'
    }
  }

  rm(fp_list) # Clean environment
  return(return_lists)
} # CLOSE read_visium_outs

##
# read_xenium_outs: Takes a list with Xenium output file paths and sample names and
# returns a list with count and coordinates data frames per sample
# @param filepaths a list with file paths to Xemium outputs and sample IDs
# @param input_rnacounts The result of detect_input_rnacounts
# @return a list with two lists within (one with counts, one with coordinates)
#
read_xenium_outs = function(filepaths, input_rnacounts, cores=NULL, verbose=TRUE){
  # Find necessary files from Xenium input
  missingSamples = 0
  fp_list = list()
  for(i in 1:length(filepaths[['count_found']])){
    fp_list[[i]] = list()
    # Get all system paths within output folder (h5/MEX, coordinates, image)
    temp_fps = list.files(filepaths[['count_found']][i], recursive=T, include.dirs=T, full.names=T)

    vcoords = grep('cells.parquet$', temp_fps, value=T)
    # Filter out 'SPATIAL_RNA_COUNTER' folders (intermediate files from Space Ranger? Also in Xenium outputs?).
    #vcoords = vcoords[!grepl('SPATIAL_RNA_COUNTER', vcoords)]

    if(input_rnacounts[1] == 'xenium_mex'){
      vfeatures = grep('cell_feature_matrix\\/features.tsv.gz',  temp_fps, value=T)
      vbarcodes = grep('cell_feature_matrix\\/barcodes.tsv.gz', temp_fps, value=T)
      vcounts = grep('cell_feature_matrix\\/matrix.mtx.gz', temp_fps, value=T)

      # Test that all files have been found.
      needed_mex_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('csv', vcoords))
      if(any(needed_mex_test)){
        raise_err(err_code='error0037')
      }

      fp_list[[i]]$features = vfeatures
      fp_list[[i]]$barcodes = vbarcodes
      fp_list[[i]]$counts = vcounts

      if(rlang::is_empty(vfeatures)) message(paste("Features for", filepaths$sampleids[i], "not able to be found..."))
      if(rlang::is_empty(vbarcodes)) message(paste("Barcodes for", filepaths$sampleids[i], "not able to be found..."))
      if(rlang::is_empty(vcounts)) message(paste("Counts for", filepaths$sampleids[i], "not able to be found..."))

      if(rlang::is_empty(vfeatures) | rlang::is_empty(vbarcodes) | rlang::is_empty(vcounts) | rlang::is_empty(vcoords)){
        fp_list[[i]] = list()
        missingSamples  = missingSamples + 1
      }

      rm(vfeatures, vbarcodes, vcounts) # Clean environment
    } else if(input_rnacounts[1] == 'xenium_h5'){
      h5counts = grep('cell_feature_matrix.h5', temp_fps, value=T)

      # Test that all files have been found.
      needed_h5_test = c(!grepl('\\.h5$', h5counts), !grepl('\\.parquet|\\.csv', vcoords))
      if(any(needed_h5_test)){
        raise_err(err_code='error0037')
      }

      fp_list[[i]]$counts = h5counts

      if(rlang::is_empty(vcoords)) message(paste("Coordinates for", filepaths$sampleids[i], "not able to be found...\n"))
      if(rlang::is_empty(h5counts) | rlang::is_empty(vcoords)){
        fp_list[[i]] = list()
        missingSamples  = missingSamples + 1
      }
      rm(h5counts) # Clean environment
    }

    fp_list[[i]]$coords = vcoords
    fp_list[[i]]$runname = filepaths[['sampleids']][i]

    rm(temp_fps, vcoords) # Clean environment
  }

  if(verbose){
    cat(paste("\tFound", length(filepaths$sampleids)-missingSamples, "Xenium samples\n"))
  }

  # Define number of available cores to use.
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(filepaths[['count_found']]))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      raise_err(err_code='error0038')
    }
  }

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    if(length(fp_list[[i]]$counts) == 0){
      return(list())
    }

    if(verbose){
      system(sprintf('echo "%s"', paste0("\t\tProcessing Sample ", i, "....")))
    }

    # Process Visium outputs
    if(input_rnacounts[1] == 'xenium_mex'){
      xenium_processed = import_xenium_mex(features_fp=fp_list[[i]][['features']],
                                       barcodes_fp=fp_list[[i]][['barcodes']],
                                       counts_fp=fp_list[[i]][['counts']],
                                       coords_fp=fp_list[[i]][['coords']])
    } else{
      xenium_processed = import_xenium_h5(counts_fp=fp_list[[i]][['counts']],
                                          coords_fp=fp_list[[i]][['coords']])
    }

    if(verbose){
      system(sprintf('echo "%s"', paste0("\t\tFinished data read Sample ", i)))
    }

    return(xenium_processed)
  }, mc.cores=cores, mc.preschedule=F)
  if(verbose){
    cat(paste("\tData read completed\n"))
  }

  # Organize the paralellized output into corresponding lists.
  return_lists = list()
  return_lists[['counts']] = list()
  return_lists[['coords']] = list()
  for(i in 1:length(output_temp)){
    return_lists[['counts']][[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
    return_lists[['coords']][[fp_list[[i]]$runname]] = output_temp[[i]]$coords
  }

  rm(fp_list) # Clean environment
  return(return_lists)
} # CLOSE read_xenium_outs

##
# read_cosmx_input: Takes a list with CosMx-SMI output file paths and sample names and
# returns a list with count and coordinates data frames per sample
# @param filepaths a list with file paths to CosMx-SMI outputs and sample IDs
# @param input_check The result of detect_input
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_cosmx_input = function(filepaths, input_check, cores=NULL, verbose=TRUE){
  # Verify that the files contain the necessary columns: fov, cell_ID, CenterX_local_px, and CenterY_local_px
  missingSamples = 0
  fp_list = list()
  for(i in 1:length(filepaths[['count_found']])){
    # Does RNA count file have the required columns?
    firstline = readLines(filepaths[['count_found']][[i]], n=1)
    rna_cols = (grepl('fov', firstline) & grepl('cell_ID', firstline))
    rm(firstline)

    # Does coordinate file have the required columns?
    firstline = readLines(filepaths[['coord_found']][[i]], n=1)
    coord_cols = (grepl('fov', firstline) & grepl('cell_ID', firstline) & grepl('CenterX_local_px', firstline) & grepl('CenterY_local_px', firstline))
    rm(firstline)

    # If columns are not present in files, count as missing
    fp_list[[i]] = list()
    if(!rna_cols | !coord_cols){
      missingSamples  = missingSamples + 1
    } else{
      fp_list[[i]][['counts']] = filepaths[['count_found']][[i]]
      fp_list[[i]][['coords']] = filepaths[['coord_found']][[i]]
      fp_list[[i]][['runname']] = filepaths[['sampleids']][[i]]
    }
  }

  if(verbose){
    cat(paste("Found", length(filepaths$sampleids)-missingSamples, "CosMx-SMI samples\n"))
  }

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(1:length(fp_list)), function(i){
    if(length(fp_list[[i]][['counts']]) == 0){
      return(list())
    }

    if(verbose){
      system(sprintf('echo "%s"', paste0("\tProcessing sample ", i, "....")))
    }

    # Process CosMx outputs.
    cosmx_processed = import_smi(counts_fp=fp_list[[i]][['counts']],
                                 coords_fp=fp_list[[i]][['coords']],
                                 slidename=fp_list[[i]][['runname']])

    if(verbose){
      system(sprintf('echo "%s"', paste0("\tFinished data read sample ", i)))
    }

    return(cosmx_processed)
  }, mc.cores=cores, mc.preschedule=F)

  if(verbose){
    cat("\tData read completed\n")
  }

  # Organize the paralellized output into corresponding lists.
  return_lists = list()
  return_lists[['counts']] = list()
  return_lists[['coords']] = list()
  for(i in 1:length(output_temp)){
    return_lists[['counts']] = append(return_lists[['counts']], output_temp[[i]][['rawcounts']])
    return_lists[['coords']] = append(return_lists[['coords']], output_temp[[i]][['coords']])
  }

  rm(fp_list) # Clean environment
  return(return_lists)
}


##
# process_lists: Takes two named lists of counts and coordinates and process gene and
# spot names before placing them within an STlist
# @param rnacounts, a named list with counts
# @param spotcoords, a named list with coordinates
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
process_lists = function(counts_df_list, coords_df_list){
  # Keep track of FOVs with zero counts
  rm_fov = c()
  # Process the count and coordinate list.
  for(i in 1:length(names(counts_df_list))){

    # Get loop current name to be processed.
    name_i = names(counts_df_list)[i]

    # One of the test data sets has column names begin with a number. To avoid this to
    # become a problem for data sets with the same issue, will add an 'x', similar to what R would do
    if(any(grepl("^[0-9]", colnames(counts_df_list[[name_i]])))){
      colnames(counts_df_list[[name_i]])[-1] = paste0('x', colnames(counts_df_list[[name_i]][, -1]))
      # Make the same addition in coordinate data to match column names in count data
      coords_df_list[[name_i]][, 1] = paste0('x', coords_df_list[[name_i]][[1]])
    }

    # Test that spot names are the same in both count and coordinate data frames.
    if(length(setdiff(colnames(counts_df_list[[name_i]])[-1], unlist(coords_df_list[[name_i]][, 1]))) != 0){
      raise_err(err_code='error0002', samplename=name_i)
    }

    array_col = names(coords_df_list[[name_i]])[3]
    # Sort coordinate data according to third column in the coordinate data frame.
    coords_df_list[[name_i]] = coords_df_list[[name_i]] %>% dplyr::arrange(!!array_col)

    # Convert to sparse matrix if not one
    if(!any(class(counts_df_list[[name_i]]) == 'dgCMatrix')){
      counts_df_list[[name_i]] = as.data.frame(counts_df_list[[name_i]])
      counts_df_list[[name_i]] = makeSparse(counts_df_list[[name_i]])
    }

    # Order column names in count data frame according to sorted coordinate data.
    if(!all(as.vector(coords_df_list[[name_i]][[1]]) %in% colnames(counts_df_list[[name_i]]))){
      warning('Not all spots/cells in coordinates are present in count data. Removing spots/cells not in counts data') # WARNING IF COUNTS ARE SUBSET OF COORDINATES
      coords_df_list[[name_i]] = coords_df_list[[name_i]][as.vector(coords_df_list[[name_i]][[1]]) %in% colnames(counts_df_list[[name_i]]), ]
    }
    counts_df_list[[name_i]] = counts_df_list[[name_i]][, as.vector(coords_df_list[[name_i]][[1]])]

    # Put column names to coordinate data (if not already there)
    if(sum(grepl('libname|ypos|xpos', colnames(coords_df_list[[name_i]]))) != 3){
      colnames(coords_df_list[[name_i]]) = c('libname', 'ypos', 'xpos')
    }
    # Force numeric to 2 and 3 column of coordinates table to ensure coordinates can be used
    coords_df_list[[name_i]][[2]] = as.numeric(coords_df_list[[name_i]][[2]])
    coords_df_list[[name_i]][[3]] = as.numeric(coords_df_list[[name_i]][[3]])

    # Get total gene counts and genes with no-zero counts
#    if(class(counts_df_list[[name_i]])[1] == 'dgCMatrix'){
      coords_df_list[[name_i]][['total_counts']] = as.vector(colSums(as.matrix(counts_df_list[[name_i]])))
      coords_df_list[[name_i]][['total_genes']] = as.vector(colSums(as.matrix(counts_df_list[[name_i]]) != 0))
#    } else{
#      coords_df_list[[name_i]][['total_counts']] = colSums(counts_df_list[[name_i]][, -1])
#      coords_df_list[[name_i]][['total_genes']] = colSums(counts_df_list[[name_i]][, -1] != 0)
#    }

    # If no counts in the entire FOV, mark for removal
    if(sum(coords_df_list[[name_i]][['total_counts']]) < 1){
      rm_fov = append(rm_fov, name_i)
      warning(paste0('No counts present in FOV ', name_i, '. Removing from data set.\n'))
    }
  }

  # Remove FOVs without counts
  if(length(rm_fov) > 0){
    counts_df_list = counts_df_list[ !(names(counts_df_list) %in% rm_fov) ]
    coords_df_list = coords_df_list[ !(names(coords_df_list) %in% rm_fov) ]
  }

  proc_return_lists = list()
  proc_return_lists[['counts']] = counts_df_list
  proc_return_lists[['coords']] = coords_df_list
  return(proc_return_lists)
}

##
# process_sample_filepaths: Gets file paths to Visium outputs or count/coord matrices
# from a sample metadata file
# @param samples, the file path to metadata file
# @param input_check, a list resulting from detect_input
# @return filepaths, a list with found file paths per each sample, and sample IDs
#
process_sample_filepaths = function(samples, input_check){
  # Check that all system paths provided exist.
  if(!file.exists(samples)){
    stop('Could not find samples file.')
  }
  # Read sample file and get sample names.
  delsample = input_check$samples[2]
  sample_file_df = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F, show_col_types=F)
  sample_names = sample_file_df[[1]]
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex')){
    for(i in sample_names){
      sample_i = as.vector(unlist(sample_file_df[, 2][grep(i, sample_file_df[[1]]), ]))
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_i[1])
    }
    filepaths[['sampleids']] = sample_names
  } else{
    for(i in sample_names){
      sample_i = as.vector(unlist(sample_file_df[, c(2, 3)][grep(i, sample_file_df[[1]]), ]))
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_i[1])
      filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_i[2])
    }
    filepaths[['sampleids']] = sample_names
  }
  return(filepaths)
}

##
# process_sample_names: Match file paths to Visium outputs or count/coord matrices
# to sample names for processing
# @param rnacounts, vector of file paths to Visium outputs or count matrices
# @param spotcoords, vector of file paths to coordinate matrices
# @param samples, vector of sample names
# @param input_rnacounts, a string resulting from detect_input_rnacounts
# @param input_spotcoords, a string resulting from detect_input_spotcoords
# @return a list with found file paths per each sample, and sample IDs
#
process_sample_names = function(rnacounts, spotcoords, samples, input_rnacounts, input_spotcoords){
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store paths to files
  filepaths = list()
  vis_xen = c('visium_filtered_h5', 'visium_raw_h5', 'visium_filtered_mex', 'visium_raw_mex', 'xenium_h5', 'xenium_mex')
  delim = c('tab_delim', 'comma_delim', 'tab_delim_cosmx', 'comma_delim_cosmx')
  if(input_rnacounts %in% vis_xen){
    for(i in samples){
      sample_i = grep(i, rnacounts, value=T, fixed=T)
      if(length(sample_i) == 1){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      } else if(length(sample_i) > 1){
        raise_err(err_code='error0001')
      } else{
        warning(paste0('Sample ', i, ' was not found among the Visium or Xenium outputs.'))
      }
    }
  } else if(input_rnacounts %in% delim && input_spotcoords %in% delim){
    for(i in samples){
      sample_count = grep(i, rnacounts, value=T, fixed=T)
      sample_coord = grep(i, spotcoords, value=T, fixed=T)
      if(length(sample_count) != 0 & length(sample_coord) != 0){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_count)
        filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_coord)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      } else{
        warning(paste0('Sample ', i, ' was not found among the count/coordinate files.'))
      }
    }
  } else{
    raise_err(err_code='error0034')
  }
  return(filepaths)
}

##
# process_sample_names_from_file: Match file paths to Visium outputs or count/coord matrices
# to sample names from a meta data file for later processing
# @params rnacounts, vector of file paths to Visium outputs or count matrices
# @params spotcoords, vector of file paths to coordinate matrices
# @param samples, the file path to metadata file
# @param input_check, a list resulting from detect_input
# @return filepaths, a list with found file paths per each sample, and sample IDs
#
process_sample_names_from_file = function(rnacounts, spotcoords, samples, input_check){
  # Check that all system paths provided exist.
  if(!file.exists(samples)){
    stop('Could not find samples file.')
  }
  # Read sample file and get sample names.
  delsample = input_check$samples[2]
  sample_file_df = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F, show_col_types=F, name_repair='minimal')
  sample_names = sample_file_df[[1]]
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex')){
      for(i in sample_names){
        sample_i = grep(i, rnacounts, value=T)
        if(length(sample_i) != 0){
          filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
          filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
        }
      }
    }
  } else if(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex')){
    for(i in sample_names){
      visium_fps = sample_file_df[[2]] # If only samples file provided, then visium directories must be in second column of file
      sample_i = grep(i, visium_fps, value=T)
      if(length(sample_i) != 0){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      }
    }
  } else{
    for(i in sample_names){
      # Get file paths to counts and coords from samples file
      rna_fps = sample_file_df[[2]]
      coords_fp = sample_file_df[[3]]
      sample_count = grep(i, rna_fps, value=T)
      sample_coord = grep(i, coords_fp, value=T)
      if(length(sample_count) != 0 & length(sample_coord) != 0){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_count)
        filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_coord)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      }
    }
  }
  return(filepaths)
}

##
# process_meta: Takes the input file path of a metadata table and return the processed
# data frame for the STlist
# @param samples, the file path to metadata file
# @param input_check, a list resulting from detect_input
# @param counts_df_list, a list with counts resulting from process_lists
# @return samples_df, a data frame to be placed within the clinical slot of the STlist
#
process_meta = function(samples, input_check, counts_df_list){

  # To prevent NOTES in R CMD check
  . = NULL

  # Get delimiter of file from input_check
  del = input_check[['samples']][2]
  # Read file.
  samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F, show_col_types=F)
  # Remove leading spaces in columan names if any
  colnames(samples_df) = gsub('^ ', '', colnames(samples_df))

  # Get sorted list names and fetch corresponding rows from the samplefile.
  sorted_names = sort(names(counts_df_list))

  # Check if all sample names in metadata table are present
  # If CosMx, expand metadata to FOV level
  if(!is.null(input_check[['rna']])){
    if(input_check[['rna']][[1]] == 'cosmx'){
      samples_df_fov = tibble::tibble(fovid=sorted_names)
      samples_df_fov[[colnames(samples_df)[[1]]]] = stringr::str_extract(sorted_names, paste0('^', samples_df[[1]], collapse='|'))
      samples_df = samples_df %>% dplyr::left_join(samples_df_fov, ., by=colnames(samples_df)[[1]])
      colnames(samples_df)[c(1,2)] = c('sample_name', 'slide_name')
      rm(samples_df_fov) # Clean env
    } else{
      samples_df = samples_df[samples_df[[1]] %in% sorted_names, ]
    }
  } else{
    samples_df = samples_df[samples_df[[1]] %in% sorted_names, ]
  }

  # Remove file paths from sample file if existent
  if(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex')){
    samples_df = samples_df[, -2]
  } else if(input_check$samples[1] == 'samplesfile_matrices'){
    samples_df = samples_df[, -c(2:3)]
  }
  # Test that number of rows in sample file correspond to number of elements in list.
  if(nrow(samples_df) != length(counts_df_list)){
    stop("The number of rows in the sample data is not the same as the number of spatial arrays.")
  }
  return(samples_df)
}

# process_meta_geomx = function(samples, input_check, counts_df_list, gmx_slide_col, gmx_meta_cols){
#   if(input_check$samples[2] == 'xls'){
#     samples_df = readxl::read_excel(samples)
#   } else {
#     # Get delimiter of file from input_check
#     del = input_check$samples[2]
#     samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F, show_col_types=F)
#   }
#   samples_df = samples_df[samples_df[[gmx_slide_col]] %in% names(counts_df_list), ]
#   samples_df = samples_df[, c(gmx_slide_col, gmx_meta_cols)]
#
#   samples_df_summ = tibble::tibble()
#   for(slide_id in unique(samples_df[[gmx_slide_col]])){
#     meta_row = c(sampleID=slide_id)
#     for(metacol in gmx_meta_cols){
#       if(is.numeric(samples_df[[metacol]])){
#         tmp_val = mean(samples_df[[metacol]][samples_df[[gmx_slide_col]] == slide_id], na.rm=T)
#       } else{
#         tmp_val = unique(samples_df[[metacol]][samples_df[[gmx_slide_col]] == slide_id])
#         tmp_val = paste(tmp_val, collapse='|')
#       }
#       names(tmp_val) = metacol
#       meta_row = append(meta_row, tmp_val)
#     }
#     samples_df_summ = dplyr::bind_rows(samples_df_summ, meta_row)
#   }
#   return(samples_df_summ)
# }

