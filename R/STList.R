##
#' @title STlist: Creation of STlist objects for spatial transcriptomics analysis
#' @description Creates an STlist object from one or multiple spatial transcriptomic samples.
#' @details
#' Objects of the S4 class STlist are the starting point of analyses in **`spatialGE`**.
#' The STlist contains data from one or multiple samples (i.e., tissue slices), and
#' results from most `spatialGE`'s functions are stored within the object.
#' \itemize{
#' \item Raw gene counts and spatial coordinates. Gene count data have genes in rows and
#' sampling units (e.g., cells, spots, ROIs) in columns. Spatial coordinates have
#' sampling units in rows and three columns: sample unit IDs, Y position, and X position.
#' \item Visium outputs from *space ranger*. The Visium directory should generally have
#' the file structure resulting from `spaceranger count`, with either a count matrix
#' represented in MEX files or a h5 file. The directory should also contain a `spatial`
#' sub-directory, with the spatial coordinates (`tissue_positions_list.csv`), and
#' optionally the high resolution tissue image and scaling factor file `scalefactors_json.json`.
#' \item CosMx-SMI outputs. Two files are required to process SMI outputs: The `exprMat` and
#' `metadata` files. Both files must contain the "fov" and "cell_ID" columns. In addition,
#' the `metadata` files must contain the "CenterX_local_px" and "CenterY_local_px" columns.
#' }
#' Optionally, the user can input a path to a file containing a table of sample-level
#' metadata (e.g., clinical outcomes, tissue type, age). This sample metadata file
#' should contain sample IDs in the first column partially matching the file names of
#' the count/coordinate file paths or Visium directories. _Note:_ The sample ID of a
#' given sample cannot be a substring of the sample ID of another sample. For example,
#' instead of using "tissue1" and "tissue12", use "tissue01" and "tissue12".
#'
#' The function uses parallelization if run in a unix system. Windows users
#' will experience longer times depending on the number of samples.
#'
#' @param rnacounts the count data which can be provided in one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing raw gene counts, one file
#' for each spatial sample. The first column contains gene names and subsequent columns
#' contain the expression data for each cell/spot. Duplicate gene names will be
#' modified using `make.unique`. Requires `spotcoords` and `samples`
#' \item File paths to Visium output directories (one per spatial sample). The directory
#' should follow the structure resulting from `spaceranger count`. The directory contains
#' the `.h5` and `spatial` sub-directory. If no `.h5` file is available, sparse
#' matrices (MEX) from `spaceranger count`. In that case a second sub-directory
#' called `filtered_feature_bc_matrix` should contain contain the `barcodes.tsv.gz`,
#' `features.tsv.gz`, and `matrix.mtx.gz` files. The `spatial` sub-directory minimally
#' contains the coordinates (`tissue_positions_list.csv`), and optionally the high
#' resolution PNG image and accompanying scaling factors (`scalefactors_json.json`).
#' Requires `samples`
#' \item The `exprMat` file for each slide of a CosMx-SMI output. The file must contain
#' the "fov" and "cell_ID" columns. The `STlist` function will separate data from each
#' FOV, since analysis in spatialGE is conducted at the FOV level. Requires `samples` and
#' `spotcoords`
#' \item A named list of data frames with raw gene counts (one data frame per spatial
#' sample). Requires `spotcoords`. Argument `samples` only needed when a file path to
#' sample metadata is the input.
# \item File path to `.dcc` files from GeoMx output. Requires `samples`
#' }
#' @param spotcoords the cell/spot coordinates. Not required if inputs are Visium
#' space ranger outputs
# or GeoMx DCC files
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing cell/spot coordinates, one
#' for each spatial sample. The files must contain three columns: cell/spot IDs, Y positions, and
#' X positions. The cell/spot IDs must match the column names for each cells/spots (columns) in
#' the gene count files
#' \item The `metadata` file for each slide of a CosMx-SMI output. The file must contain
#' the "fov", "cell_ID", "CenterX_local_px", and "CenterY_local_px" columns. The `STlist`
#' function will separate data from each FOV, since analysis in spatialGE is conducted at
#' the FOV level. Requires `samples` and `rnacounts`
#' \item A named list of data frames with cell/spot coordinates. The list names must
#' match list names of the gene counts list
#' }
#' @param samples the sample names/IDs and (optionally) metadata associated with
#' each spatial sample. This file can also include files paths to gene counts and
#' cell/spot coordinate files, bypassing the need to specify `rnacounts` and `spotcoords`.
#' The following options are available for `samples`:
#' \itemize{
#' \item A vector with sample names, which will be used to partially match gene
#' counts and cell/spot coordinates file paths. A sample name must not match file
#' paths for two different samples. For example, instead of using "tissue1" and
#' "tissue12", use "tissue01" and "tissue12".
#' \item A path to a file containing a table with metadata. This file is a comma- or
#' tab-separated table with one sample per row and sample names/IDs in the first
#' column. Paths to gene counts and coordinate files can be placed in the second and
#' third columns respectively (while leaving the `rnacounts` and `spotcoords` arguments
#' empty). If Visium directories are provided, only the second column with paths to
#' `spaceranger count` directories are expected. Subsequent columns can contain
#' variables associated with each spatial sample
# Note: For GeoMx, the metadata file contains one row per ROI. This information is
# automatically summarized to one row per tissue slice. This metadata file usually
# also contains the X and Y positions, which can be identified via `gmx_y_col` and
# `gmx_x_col` arguments
#' }
# @param gmx_pkc the file path to the `.pkc` probe set file (for GeoMx input)
# @param gmx_slide_col the name of the column in the metadata table containing
# the slide names (for GeoMx input)
# @param gmx_roi_col the name of the column in the metadata table containing the
# ROI IDs, matching IDs in the DCC files (for GeoMx input)
# @param gmx_y_col the name of the column in the metadata table containing the
# Y positions (for GeoMx input)
# @param gmx_x_col the name of the column in the metadata table containing the
# X positions (for GeoMx input)
# @param gmx_meta_cols a vector with column names in the metadata table containing
# clinical data (for GeoMx input)
#' @param cores integer indicating the number of cores to use during parallelization.
#' If NULL, the function uses half of the available cores at a maximum. The parallelization
#' uses `parallel::mclapply` and works only in Unix systems.
#' @return x an STlist object containing the counts and coordinates, and optionally
#' the sample metadata, which can be used for downstream analysis with `spatialGE`
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
#' count_files <- grep("counts", data_files, value=T)
#' coord_files <- grep("mapping", data_files, value=T)
#' clin_file <- grep("thrane_clinical", data_files, value=T)
#' melanoma <- STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples
#' melanoma
#
#' @export STlist
#'
#' @import Matrix
#' @importFrom magrittr %>%
#'
STlist = function(rnacounts=NULL, spotcoords=NULL, samples=NULL,
#                   gmx_pkc=NULL, gmx_slide_col=NULL, gmx_roi_col=NULL, gmx_y_col=NULL, gmx_x_col=NULL, gmx_meta_cols=NULL, ### DISABLED UNTIL SUPPORT FOR GEOMX-DCC IS SORTED OUT
                  cores=NULL){
  # Check input type.
  input_check = detect_input(rnacounts=rnacounts, spotcoords=spotcoords, samples=samples)

  # Define number of available cores to use.
  if(is.null(cores)){
    if(!is.null(rnacounts)){
      cores = count_cores(length(rnacounts))
    } else{ # In case samplefile including filepaths is provided
      cores_tmp = length(readLines(samples))
      cores = count_cores(cores_tmp)
      rm(cores_tmp) # Clean env
    }
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Output error if input_check is empty (likely input format not recognized).
  if(rlang::is_empty(input_check)){
    stop('Input not recognized. Please refer to documentation.')
  }

  # CASE: NAMED LIST OF DATAFRAMES WITH COUNTS AND NAMED LIST OF DATA FRAMES WITH COORDINATES.
  # METADATA INFO OPTIONAL.
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'list_dfs' && input_check$coords[1] == 'list_dfs'){
      cat(crayon::green(paste("Found list of dataframes.\n")))
      pre_lists = read_list_dfs(rnacounts, spotcoords)
    }
  }

  # CASE: SEURAT OBJECT
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'seurat'){
      cat(crayon::green(paste("Found Seurat object.\n")))
      pre_lists = read_seurat(rnacounts)
      img_obj = pre_lists[['images']]
      platform = 'visium'
    }
  }

  # CASE: GEOMX INPUT
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'geomx_dcc'){
      cat(crayon::green(paste("Found GeoMx DCC output.\n")))
      pre_lists = import_Geomx(dcc=rnacounts, pkc=gmx_pkc, annots=samples, slide_col=gmx_slide_col, id_col=gmx_roi_col, x_col=gmx_x_col, y_col=gmx_y_col)
      platform = 'geomx'
    }
  }

  # CASE: SAMPLE FILE ONLY CONTAINING FILE PATH(S) TO COUNT COORDINATE MATRICES OR VISIUM DIRS
  if(is.null(rnacounts) && is.null(input_check$coords) && !is.null(input_check$samples)){
    # Get list of filepaths
    filepaths = process_sample_names_from_file(samples=samples, input_check=input_check)
    # Check if input is Visium or count/coord matrices
    if(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex')){
      cat(crayon::green(paste("Found Visium data\n")))
      pre_lists = read_visium_outs(filepaths, input_check=input_check, cores=cores)
      img_obj = pre_lists[['images']]
      image_scale = pre_lists[['json_scale']]
      platform = 'visium'
    } else{
      cat(crayon::green(paste("Found matrix data\n")))
      pre_lists = read_matrices_fps(filepaths, input_check, cores=cores)
    }
  }

  # CASE: SAMPLE FILE PLUS FILE PATH(S) TO COUNT COORDINATE MATRICES OR VISIUM DIRS
  if(!is.null(rnacounts) && (input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex') || input_check$samples[1] == 'samplesfile')){
    if(input_check$rna[1] != 'list_dfs'){
      # Read sample nanes
      sample_names = readr::read_delim(samples, show_col_types=F)
      sample_names = as.character(sample_names[[1]])
      # Get list of filepaths
      filepaths = process_sample_names(rnacounts, spotcoords, sample_names, input_check)

      # Check if input is Visium, CosMx, or count/coord matrices
      if(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex')){
        cat(crayon::green(paste("Found Visium data\n")))
        pre_lists = read_visium_outs(filepaths, input_check, cores=cores)
        img_obj = pre_lists[['images']]
        image_scale = pre_lists[['json_scale']]
        platform = 'visium'
      } else if(input_check$rna[1] == 'cosmx'){
        cat(crayon::green(paste("Found CosMx-SMI data\n")))
        pre_lists = read_cosmx_input(filepaths, input_check, cores=cores)
        img_obj = pre_lists[['images']]
        platform = 'cosmx'
      } else{
        cat(crayon::green(paste("Found matrix data\n")))
        pre_lists = read_matrices_fps(filepaths, input_check, cores=cores)
      }
    }
  }

  # CASE: FILE PATH(S) TO COUNT/COORDINATES MATRICES OR VISIUM DIRS, AND SAMPLE NAMES VECTOR.
  if(any(input_check$samples == 'sample_names') && !is.null(rnacounts)){
    if(input_check$rna[1] != 'list_dfs'){
      # Get list of filepaths
      filepaths = process_sample_names(rnacounts, spotcoords, as.character(samples), input_check)

      # Check if input is Visium, CosMx, or count/coord matrices
      if(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex')){
        cat(crayon::green(paste("Found Visium data\n")))
        pre_lists = read_visium_outs(filepaths, input_check, cores=cores)
        img_obj = pre_lists[['images']]
        image_scale = pre_lists[['json_scale']]
        platform = 'visium'
      } else if(input_check$rna[1] %in% c('cosmx')){
        cat(crayon::green(paste("Found CosMx-SMI data\n")))
        pre_lists = read_cosmx_input(filepaths, input_check, cores=cores)
        img_obj = pre_lists[['images']]
        platform = 'cosmx'
      } else{
        cat(crayon::green(paste("Found matrix data\n")))
        pre_lists = read_matrices_fps(filepaths, input_check, cores=cores)
      }
    }
  }

  # No input provided
  if(is.null(input_check$rna) && is.null(input_check$coords) && is.null(input_check$samples)){
    stop('No input provided. Please refer to documentation.')
  }

  # Process count and coordinate lists before placing within STlist
  cat(crayon::yellow(paste("Matching gene expression and coordinate data...\n")))
  procLists = process_lists(counts_df_list=pre_lists[['counts']], coords_df_list=pre_lists[['coords']])

  # Process metadata if provided or make an empty tibble
  samples_df = tibble::tibble()
  if(input_check$samples[1] == 'samplesfile' || input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex') || input_check$samples[1] == 'samplesfile_matrices'){
    samples_df = process_meta(samples=samples, input_check=input_check, counts_df_list=procLists[['counts']])
  }else if(input_check$samples[1] == 'samplesfile_geomx'){
    samples_df = process_meta_geomx(samples, input_check, procLists[['counts']], gmx_slide_col, gmx_meta_cols)
  }else{
    samples_df = tibble::tibble(sample_name=names(procLists[['counts']]))
  }
  # Make sure sample IDs are character
  samples_df[[1]] = as.character(samples_df[[1]])

  if(!is.null(input_check$rna[1])){
    if(!(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex'))){
      cat(crayon::yellow(paste("Converting counts to sparse matrices\n")))
      procLists[['counts']] = parallel::mclapply(procLists[['counts']], function(x){
        makeSparse(x)
      })
    }
  } else if(!(input_check$samples[1] %in% c('samplesfile_visium_h5', 'samplesfile_visium_mex'))){
    cat(crayon::yellow(paste("Converting counts to sparse matrices\n")))
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
  cat(crayon::green$bold(paste("Completed STlist!\n")))
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
  #suppressMessages({library(Matrix)})
  if(!is.matrix(dataframe)){
    genecol = colnames(dataframe)[1]
    # numdat = dataframe %>%
    #   tibble::column_to_rownames(genecol) %>%
    #   as.matrix() %>% as(., "sparseMatrix")
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

  # Define number of available cores to use.
  if(is.null(cores)){
    cores = count_cores(length(filepaths[['count_found']]))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Use parallelization to read count data if possible.
  counts_df_list = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    # Read filepaths.
    counts_df = readr::read_delim(filepaths[['count_found']][i], delim=delrna, col_types=readr::cols(), progress=F)
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
# @param input_check The result of detect_input
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_visium_outs = function(filepaths, input_check, cores=NULL){
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
      grep('tissue_positions_list.csv', ., value=T)
    # Filter out 'SPATIAL_RNA_COUNTER' folders (intermediate files from Space Ranger?).
    vcoords = vcoords[!grepl('SPATIAL_RNA_COUNTER', vcoords)]

    # Check if data is MEX of HDF5
    if(!is.null(input_check$rna)){
      if(input_check$rna[1] == 'visium_out_mex'){
        mtx_type = 'mex'
      } else if(input_check$rna[1] == 'visium_out_h5'){
        mtx_type = 'h5'
      }
    } else{
      if(input_check$samples[1] == 'samplesfile_visium_mex'){
        mtx_type = 'mex'
      } else if(input_check$samples[1] == 'samplesfile_visium_h5'){
        mtx_type = 'h5'
      }
    }

    if(mtx_type == 'mex'){
      # NOTE: For MEX, assume files are inside a directory named 'filtered_feature_bc_matrix'
      vfeatures = grep('[raw|filtered]_feature_bc_matrix\\/features.tsv.gz',  temp_fps, value=T)
      vbarcodes = grep('[raw|filtered]_feature_bc_matrix\\/barcodes.tsv.gz', temp_fps, value=T)
      vcounts = grep('[raw|filtered]_feature_bc_matrix\\/matrix.mtx.gz', temp_fps, value=T)

      # In case there are both raw and filtered data within output directory, choose filtered
      if(length(vfeatures) > 1){
        vfeatures = grep("filtered_feature_bc_matrix", vfeatures, value=T)
        vbarcodes = grep("filtered_feature_bc_matrix", vbarcodes, value=T)
        vcounts = grep("filtered_feature_bc_matrix", vcounts, value=T)
      }

      # Test that all files have been found.
      needed_mex_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('csv', vcoords))
      if(any(needed_mex_test)){
        stop(paste0('Visium output folder (', filepaths[['count_found']][i], ') does not have all necessary files.'))
      }

      fp_list[[i]]$features = vfeatures
      fp_list[[i]]$barcodes = vbarcodes
      fp_list[[i]]$counts = vcounts

      if(rlang::is_empty(vfeatures)) cat(crayon::red(paste("Features for", filepaths$sampleids[i], "not able to be found...")))
      if(rlang::is_empty(vbarcodes)) cat(crayon::red(paste("Barcodes for", filepaths$sampleids[i], "not able to be found...")))
      if(rlang::is_empty(vcounts)) cat(crayon::red(paste("Counts for", filepaths$sampleids[i], "not able to be found...")))

      if(rlang::is_empty(vfeatures) | rlang::is_empty(vbarcodes) | rlang::is_empty(vcounts) | rlang::is_empty(vcoords)){
        fp_list[[i]] = list()
        missingSamples  = missingSamples + 1
      }

      rm(vfeatures, vbarcodes, vcounts) # Clean environment
    } else if(mtx_type == 'h5'){
      h5counts = grep('[raw|filtered]_feature_bc_matrix.h5', temp_fps, value=T)

      # In case there are both raw and filtered data within output directory, choose filtered
      if(length(h5counts) > 1){
        h5counts = grep("filtered_feature_bc_matrix", h5counts, value=T)
      }

      # Test that all files have been found.
      needed_h5_test = c(!grepl('\\.h5$', h5counts), !grepl('csv', vcoords))
      if(any(needed_h5_test)){
        stop(paste0('Visium output folder (', filepaths[['count_found']][i], ') does not have all necessary files.'))
      }

      fp_list[[i]]$counts = h5counts

      if(rlang::is_empty(vcoords)) cat(crayon::red(paste("Coordinates for", filepaths$sampleids[i], "not able to be found...\n")))
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

  cat(crayon::green$bold(paste("Found", length(filepaths$sampleids)-missingSamples, "Visium samples\n")))

  # Define number of available cores to use.
  if(is.null(cores)){
    cores = count_cores(length(filepaths[['count_found']]))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    if(length(fp_list[[i]]$counts) == 0){
      return(list())
    }
    system(sprintf('echo "%s"', crayon::yellow(paste0("\tProcessing Sample ", i, "...."))))

    # Process Visium outputs.
    if(mtx_type == 'mex'){
      visium_processed = import_visium(features_fp=fp_list[[i]][['features']],
                                       barcodes_fp=fp_list[[i]][['barcodes']],
                                       counts_fp=fp_list[[i]][['counts']],
                                       coords_fp=fp_list[[i]][['coords']])
    } else{
      visium_processed = import_visium_h5(counts_fp=fp_list[[i]][['counts']],
                                          coords_fp=fp_list[[i]][['coords']])
    }
    system(sprintf('echo "%s"', crayon::green(paste0("\tFinished data read Sample ", i))))
    return(visium_processed)
  }, mc.cores=cores, mc.preschedule=F)
  cat(crayon::green$bold(paste("\tData read completed\n")))

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
}


##
# read_cosmx_input: Takes a list with CosMx-SMI output file paths and sample names and
# returns a list with count and coordinates data frames per sample
# @param filepaths a list with file paths to CosMx-SMI outputs and sample IDs
# @param input_check The result of detect_input
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_cosmx_input = function(filepaths, input_check, cores=NULL){
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

  cat(crayon::green$bold(paste("Found", length(filepaths$sampleids)-missingSamples, "CosMx-SMI samples\n")))

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(1:length(fp_list)), function(i){
    if(length(fp_list[[i]][['counts']]) == 0){
      return(list())
    }
    system(sprintf('echo "%s"', crayon::yellow(paste0("\tProcessing sample ", i, "...."))))

    # Process CosMx outputs.
    cosmx_processed = import_smi(counts_fp=fp_list[[i]][['counts']],
                                 coords_fp=fp_list[[i]][['coords']],
                                 slidename=fp_list[[i]][['runname']])

    system(sprintf('echo "%s"', crayon::green(paste0("\tFinished data read sample ", i))))
    return(cosmx_processed)
  }, mc.cores=cores, mc.preschedule=F)
  cat(crayon::green$bold(paste("\tData read completed\n")))

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

    # Order column names in count data frame according to sorted coordinate data.
    if(class(counts_df_list[[name_i]])[1] == 'dgCMatrix'){
      counts_df_list[[name_i]] = counts_df_list[[name_i]][, as.vector(unlist(coords_df_list[[name_i]][, 1]))]
    } else{
      counts_df_list[[name_i]] = counts_df_list[[name_i]][, c(colnames(counts_df_list[[name_i]])[1], coords_df_list[[name_i]][[1]])]
    }

    # Put column names to coordinate data (if not already there)
    if(sum(grepl('libname|ypos|xpos', colnames(coords_df_list[[name_i]]))) != 3){
      colnames(coords_df_list[[name_i]]) = c('libname', 'ypos', 'xpos')
    }
    # Force numeric to 2 and 3 column of coordinates table to ensure coordinates can be used
    coords_df_list[[name_i]][[2]] = as.numeric(coords_df_list[[name_i]][[2]])
    coords_df_list[[name_i]][[3]] = as.numeric(coords_df_list[[name_i]][[3]])

    # Get total gene counts and genes with no-zero counts
    if(class(counts_df_list[[name_i]])[1] == 'dgCMatrix'){
      coords_df_list[[name_i]][['total_counts']] = colSums(as.matrix(counts_df_list[[name_i]]))
      coords_df_list[[name_i]][['total_genes']] = colSums(as.matrix(counts_df_list[[name_i]]) != 0)
    } else{
      coords_df_list[[name_i]][['total_counts']] = colSums(counts_df_list[[name_i]][, -1])
      coords_df_list[[name_i]][['total_genes']] = colSums(counts_df_list[[name_i]][, -1] != 0)
    }
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
# @param sample_names, vector of sample names
# @param input_check, a list resulting from detect_input
# @return filepaths, a list with found file paths per each sample, and sample IDs
#
process_sample_names = function(rnacounts, spotcoords, samples, input_check){
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(input_check$rna[1] %in% c('visium_out_h5', 'visium_out_mex')){
    for(i in samples){
      sample_i = grep(i, rnacounts, value=T)
      if(length(sample_i) == 1){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      } else if(length(sample_i) > 1){
        raise_err(err_code='error0001')
      } else{
        warning(paste0('Sample ', i, ' was not found among the Visium outputs.'))
      }
    }
  } else{
    for(i in samples){
      sample_count = grep(i, rnacounts, value=T)
      sample_coord = grep(i, spotcoords, value=T)
      if(length(sample_count) != 0 & length(sample_coord) != 0){
        filepaths[['count_found']] = append(filepaths[['count_found']], sample_count)
        filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_coord)
        filepaths[['sampleids']] = append(filepaths[['sampleids']], i)
      } else{
        warning(paste0('Sample ', i, ' was not found among the count/coordinate files.'))
      }
    }
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

process_meta_geomx = function(samples, input_check, counts_df_list, gmx_slide_col, gmx_meta_cols){
  if(input_check$samples[2] == 'xls'){
    samples_df = readxl::read_excel(samples)
  } else {
    # Get delimiter of file from input_check
    del = input_check$samples[2]
    samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F, show_col_types=F)
  }
  samples_df = samples_df[samples_df[[gmx_slide_col]] %in% names(counts_df_list), ]
  samples_df = samples_df[, c(gmx_slide_col, gmx_meta_cols)]

  samples_df_summ = tibble::tibble()
  for(slide_id in unique(samples_df[[gmx_slide_col]])){
    meta_row = c(sampleID=slide_id)
    for(metacol in gmx_meta_cols){
      if(is.numeric(samples_df[[metacol]])){
        tmp_val = mean(samples_df[[metacol]][samples_df[[gmx_slide_col]] == slide_id], na.rm=T)
      } else{
        tmp_val = unique(samples_df[[metacol]][samples_df[[gmx_slide_col]] == slide_id])
        tmp_val = paste(tmp_val, collapse='|')
      }
      names(tmp_val) = metacol
      meta_row = append(meta_row, tmp_val)
    }
    samples_df_summ = dplyr::bind_rows(samples_df_summ, meta_row)
  }
  return(samples_df_summ)
}

