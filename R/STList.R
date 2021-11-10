##
#' @title STList: Creation of STList objects
#' @description Creates an STList object from one or several spatial transcriptomic samples.
#' @details
#' Objects of the S4 class STList can be created from two sources:
#' \itemize{
#' \item Raw RNA counts (genes in rows, spots in columns) and spot coordinates (spots
#' in rows, with three columns: spot IDs, y position, and x position).
#' \item Visium output folders from Space Ranger. The folder should have the
#' structure resulting from `spaceranger count`.
#'}
#' Optionally, the user can input a system path to a table containing data associated with
#' each spatial array (e.g., clinical variables). This sample metadata file should contain
#' sample IDs in the first column matching partially the file names of the
#' count/coordinate file paths or Visium directories.
#'
#' The function will read data in parallel if unix system is available (Windows users
#' will experience longer times depending on the number of samples).
#'
#' @param rnacounts The count data which can be provided in one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing raw RNA counts, one
#' for each spatial array. The first column contains gene names. Subsequent columns
#' contain data for each spot. Duplicate gene names will be appended a suffix ('_d**').
#' Requires `spotcoords` and `samples`.
#' \item File paths to Visium output folders (one per spatial array). The folder should have the
#' structure resulting from `spaceranger count` and contain the `.h5` files and `spatial` folder.
#' Requires `samples`.
#' \item File path to `.dcc` files from GeoMx output. Requires `samples`.
#' \item One named list of data frames with raw RNA counts (one data frame per spatial array).
#' Requires `spotcoords` and `samples`.
#' }
#' @param spotcoords The spot coordinates, which are required if inputs are not Visium outputs.
#' The user can provide one of these formats:
#' \itemize{
#' \item File paths to comma- or tab-delimited files containing spot coordinates, one
#' for each spatial array. The files must contain three columns: spot IDs, y positions, and
#' x positions. The spot IDs must match the column names for each spot (column) in the RNA count files.
#' \item One named list of data frames with spot coordinates. The list names must match list
#' names of the RNA counts list.
#' }
#' @param samples The metadata associated to each spatial array. This file can also include
#' system paths to RNA count and spot coordinate files, bypassing the need to specify `rnacounts`
#' and `spotcoords`. One of the following options should be entered to create an STList:
#' \itemize{
#' \item A vector with sample names, which will be used to partial match RNA count and
#' spot coordinates file paths. A sample name must not match two file paths.
#' \item A path to  the file containing metadata for each spatial array. This file is a comma-
#' or tab-separated file with one sample per row and sample names in the first column.
#' Paths to RNA count and coordinate files can be in the second and third column respectively,
#' and omitting the `rnacounts` and `spotcoords` arguments. If Visium, only the second
#' column with paths to output folders is expected. Subsequent columns can contain
#' variables associated to each spatial array.
#' Note: For GeoMx, the metadata file contains one row per ROI. This information is later summarized
#' to one row per slide by this function.
#' }
#' @param gmx_pkc, the file path to the `.pkc` (for GeoMx input)
#' @param gmx_slide_col, the name of the column in the metadata table containing the slide names (for GeoMx input)
#' @param gmx_roi_col, the name of the column in the metadata table containing the ROI IDs, matching IDs in the DCC files (for GeoMx input)
#' @param gmx_x_col, the name of the column in the metadata table containing the x coordinates (for GeoMx input)
#' @param gmx_y_col, the name of the column in the metadata table containing the y coordinates (for GeoMx input)
#' @param gmx_meta_cols, a vector with column names in the metadata table containing clinical data (for GeoMx input)
#' @return x, the STList object containing the counts and coordinates, and optionally
#' the sample metadata.
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' # data_files <- system.file("extdata", package="spatialGEdev")
#' # count_files <- grep("genes", data_files, value=T)
#' # coord_files <- grep("mapping", data_files, value=T)
#' # clin_file <- grep("clinical", data_files, value=T)
#' # melanoma <- STList(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
#' # melanoma
#
#' @export STList
<<<<<<< HEAD
STList = function(rnacounts=NULL, spotcoords=NULL, samples=NULL) {
  cat(crayon::yellow(paste("Requested", length(samples), "Samples\n")))
=======
#'
STList = function(rnacounts=NULL, spotcoords=NULL, samples=NULL, gmx_pkc=NULL, gmx_slide_col=NULL, gmx_roi_col=NULL, gmx_x_col=NULL, gmx_y_col=NULL, gmx_meta_cols=NULL) {
  require('magrittr')
>>>>>>> oscar_dev
  # Check input type.
  input_check = detect_input(rnacounts=rnacounts, spotcoords=spotcoords, samples=samples)

  # Output error if input_check is empty (likely input format not recognized).
  if(rlang::is_empty(input_check)){
    stop('Input not recognized. Please refer to documentation.')
  }

  # CASE: NAMED LIST OF DATAFRAMES WITH COUNTS AND NAMED LIST OF DATA FRAMES WITH COORDINATES.
  # METADATA INFO OPTIONAL.
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'list_dfs' && input_check$coords[1] == 'list_dfs'){
      cat(crayon::yellow(paste("Found List of Dataframes...\n")))
      pre_lists = read_list_dfs(rnacounts, spotcoords)
    }
  }

  # CASE: SEURAT OBJECT
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'seurat'){
<<<<<<< HEAD
      cat(crayon::yellow(paste("Found Seurat Data.\n")))
=======
      cat(crayon::yellow(paste("Found Seurat object.\n")))
>>>>>>> oscar_dev
      pre_lists = read_seurat(rnacounts)
      img_obj = pre_lists[['images']]
      platform = 'visium'
    }
  }

  # CASE: GEOMX INPUT
  if(!is.null(input_check$rna)){
    if(input_check$rna[1] == 'geomx_dcc'){
      cat(crayon::yellow(paste("Found GeoMx DCC output.\n")))
      pre_lists = import_Geomx(dcc=rnacounts, pkc=gmx_pkc, annots=samples, slide_col=gmx_slide_col, id_col=gmx_roi_col, x_col=gmx_x_col, y_col=gmx_y_col)
      platform = 'geomx'
    }
  }

  # CASE: SAMPLE FIlE ONLY CONTAINING FILE PATH(S) TO COUNT COORDINATE MATRICES OR VISIUM DIRS
  if(is.null(rnacounts) && is.null(input_check$coords) && !is.null(input_check$samples)){
    # Get list of filepaths
    filepaths = process_sample_filepaths(samples, input_check)
    # Check if input is Visium or count/coord matrices
    if(input_check$samples[1] == 'samplesfile_visium'){
      cat(crayon::yellow(paste("Found Visium Data.\n")))
      pre_lists = read_visium_outs(filepaths)
      img_obj = pre_lists[['images']]
      platform = 'visium'
    } else{
      cat(crayon::yellow(paste("Found Matrix Data.\n")))
      pre_lists = read_matrices_fps(filepaths)
    }
  }

  # CASE: SAMPLE FIlE PLUS FILE PATH(S) TO COUNT COORDINATE MATRICES OR VISIUM DIRS
  if(!is.null(rnacounts) && (input_check$samples[1] == 'samplesfile_visium' || input_check$samples[1] == 'samplesfile')){
    if(input_check$rna[1] != 'list_dfs'){
      # Get list of filepaths
      filepaths = process_sample_names_from_file(rnacounts, spotcoords, samples, input_check)
      # Check if input is Visium or count/coord matrices
      if(input_check$rna[1] == 'visium_out'){
        cat(crayon::yellow(paste("Found Visium Data.\n")))
        pre_lists = read_visium_outs(filepaths)
        img_obj = pre_lists[['images']]
        platform = 'visium'
      } else{
        cat(crayon::yellow(paste("Found Matrix Data.\n")))
        pre_lists = read_matrices_fps(filepaths)
      }
    }
  }

  # CASE: FILE PATH(S) TO COUNT/COORDINATES MATRICES OR VISIUM DIRS, AND SAMPLE NAMES VECTOR.
  if(input_check$samples == 'sample_names' && !is.null(rnacounts)){
    if(input_check$rna[1] != 'list_dfs'){
      # Get list of filepaths
      filepaths = process_sample_names(rnacounts, spotcoords, samples, input_check)
      # Check if input is Visium or count/coord matrices
      if(input_check$rna[1] == 'visium_out'){
        cat(crayon::yellow(paste("Found Visium Data.\n")))
        pre_lists = read_visium_outs(filepaths)
        img_obj = pre_lists[['images']]
        platform = 'visium'
      } else{
        cat(crayon::yellow(paste("Found Matrix Data.\n")))
        pre_lists = read_matrices_fps(filepaths)
      }

    }
  }

  # No input provided
  if(is.null(input_check$rna) && is.null(input_check$coords) && is.null(input_check$samples)){
    stop('No input provided. Please refer to documentation.')
  }

  cat(crayon::yellow(paste("Requested", length(pre_lists[['counts']]), "Samples\n")))

  # Process count and coordinate lists before placing within STList
  cat(crayon::yellow(paste("Cleaning Count and Coordinate Data Gene Names.\n")))
  procLists = process_lists(pre_lists[['counts']], pre_lists[['coords']])

  # Process metadata if provided or make an empty tibble
  samples_df = tibble::tibble()
  if(input_check$samples[1] == 'samplesfile' || input_check$samples[1] == 'samplesfile_visium' || input_check$samples[1] == 'samplesfile_matrices'){
    samples_df = process_meta(samples, input_check, procLists[['counts']])
  }else if(input_check$samples[1] == 'samplesfile_geomx'){
    samples_df = process_meta_geomx(samples, input_check, procLists[['counts']], gmx_slide_col, gmx_meta_cols)
  }else{
    samples_df = tibble::tibble()
  }

  cat(crayon::yellow(paste("Converting Counts to Sparse Matrices\n")))
<<<<<<< HEAD
  procLists[['counts']] = mclapply(procLists[['counts']], function(x){
    makeSparse(x)
  })

=======
  procLists[['counts']] = parallel::mclapply(procLists[['counts']], function(x){
    makeSparse(x)
  })

  # Detect if image from Visium out is available
  if(!exists('img_obj')){
    img_obj = NULL
  }

  # If no specific platform was found, then make generic
  if(!exists('platform')){
    platform = 'generic'
  }

>>>>>>> oscar_dev
  # Creates STList object from both count and coordinates data.
  STList_obj = new("STList",
                   counts=procLists[['counts']],
                   coords=procLists[['coords']],
                   clinical=samples_df,
                   tr_counts=list(),
                   #voom_counts=list(),
                   #log_counts=list(),
                   gene_var=list(),
                   #log_stdev=list(),
                   gene_het=list(),
                   gene_krige=list(),
                   cell_deconv=list(),
                   cell_het=list(),
                   cell_krige=list(),
                   gene_krige_data=list(),
                   deconv_krige_data=list(),
                   st_clusters=list(),
                   pheno_plots=list(),
                   misc=list(sp_images=img_obj, platform=platform)
  )
  cat(crayon::green$bold(paste("Completed STList!\n")))
  return(STList_obj)
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
  suppressMessages({library(Matrix)})
  numdat = dataframe %>%
<<<<<<< HEAD
    column_to_rownames("gene") %>%
=======
    tibble::column_to_rownames("gene") %>%
>>>>>>> oscar_dev
    as.matrix() %>% as(., "sparseMatrix")
}

##
# expandSparse: takes a sparsed matrix and returns a dataframe
# @param sparsedMatrix, a sparse matrix made from the Matrix
# package
# @return a data frame
#
expandSparse = function(sparsedMatrix){
<<<<<<< HEAD
  sparsedList %>% data.frame(check.names=F)
=======
  NonSparse = data.frame(as.matrix(sparsedMatrix), check.names=F)
  return(NonSparse)
>>>>>>> oscar_dev
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
  # NOTE: Duplicate genes coming from matrices are only appended a '_d[0-9]+', but
  # are not comparable across matrices because unique identifiers for the probes are not used.
  for(i in 1:length(rnacounts)){
    dup_genes_Mask = duplicated(rnacounts[[i]][[1]])
    dup_genes = rnacounts[[i]][[1]][dup_genes_Mask]
    if(length(dup_genes) != 0){
      for(gene in 1:length(dup_genes)){
        geneToChange = rnacounts[[i]][[1]][rnacounts[[i]][[1]] == dup_genes[gene]]
        if(!is.null(geneToChange) & !is.na(geneToChange)){
          rnacounts[[i]][[1]][rnacounts[[i]][[1]] == dup_genes[gene]] = paste0(geneToChange, '_d', gene)
        }
      }
    }
  }

  # Create list to be returned
  return_lists = list()
  return_lists[['counts']] = rnacounts[sorted_names]
  return_lists[['coords']] = spotcoords[sorted_names]
  return(return_lists)
}

##
# read_seurat: Takes a Seurat object and converts it to a STList.
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
    return_lists[['coords']][[i]] = rnacounts@images[[i]]@coordinates[, c('row', 'col')]
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
# @param filepaths, a list with file paths to Visium outputs and sample IDs
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_matrices_fps = function(filepaths){
  # Get delimiter of RNA counts and coordinates.
  # NOTE: This only test the first matrix.
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

  # Define number of available cores to use.
  cores = count_cores(length(filepaths[['count_found']]))

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

  # Process duplicate gene names
  # NOTE: Duplicate genes coming from matrices are only appended a '_d[0-9]+', but
  # are not comparable across matrices because unique identifiers for the probes are not used.
  for(i in 1:length(counts_df_list)){
    dup_genes_idx = which(duplicated(counts_df_list[[i]][[1]]))
    dup_genes_names = counts_df_list[[i]][[1]][dup_genes_idx]
    if(length(dup_genes_idx) != 0){
      for(gene in 1:length(dup_genes_idx)){
        #geneToChange = counts_df_list[[i]][[1]][counts_df_list[[i]][[1]] == dup_genes[gene]]
        #if(!is.null(geneToChange) & !is.na(geneToChange)){
        #counts_df_list[[i]][[1]][counts_df_list[[i]][[1]] == dup_genes[gene]] = paste0(geneToChange, '_d', gene)
        counts_df_list[[i]][[1]][dup_genes_idx[gene]] = paste0(counts_df_list[[i]][[1]][dup_genes_idx[gene]], '_d', gene)
        #}
      }
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
# @param filepaths, a list with file paths to Visium outputs and sample IDs
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
read_visium_outs = function(filepaths){
  # Find necessary files from visium input
  missingSamples = 0
  fp_list = list()
  for(i in 1:length(filepaths[['count_found']])){
    # Get all system paths within output folder.
    temp_fps = list.files(filepaths[['count_found']][i], recursive=T, include.dirs=T, full.names=T)
    # Find matches of relevant files wihtin output folders.
    vfeatures = grep('filtered_feature_bc_matrix/features.tsv.gz',  temp_fps, value=T)
    vbarcodes = grep('filtered_feature_bc_matrix/barcodes.tsv.gz', temp_fps, value=T)
    vcounts = grep('filtered_feature_bc_matrix/matrix.mtx.gz', temp_fps, value=T)
    vcoords = grep('spatial/tissue_positions_list.csv', temp_fps, value=T)
    vimage = grep('spatial/tissue_lowres_image.png', temp_fps, value=T)

    # Filter out 'SPATIAL_RNA_COUNTER' folders (intermediate files from Space Ranger?).
    vcoords = vcoords[!grepl('SPATIAL_RNA_COUNTER', vcoords)]

    # Test that all files have been found.
    needed_fps_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('csv', vcoords))
    if(any(needed_fps_test)){
      stop(paste0('Visium output folder (', filepaths[['count_found']][i], ') does not have all necessary files.'))
    }

    fp_list[[i]] = list()
    fp_list[[i]]$features = vfeatures
    fp_list[[i]]$barcodes = vbarcodes
    fp_list[[i]]$counts = vcounts
    fp_list[[i]]$coords = vcoords
    fp_list[[i]]$image = vimage
    fp_list[[i]]$runname = filepaths[['sampleids']][i]
<<<<<<< HEAD
    if(is_empty(vfeatures)) cat(crayon::red(paste("Features for", filepaths$sampleids[i], "not able to be found...")))
    if(is_empty(vbarcodes)) cat(crayon::red(paste("Barcodes for", filepaths$sampleids[i], "not able to be found...")))
    if(is_empty(vcounts)) cat(crayon::red(paste("Counts for", filepaths$sampleids[i], "not able to be found...")))
    if(is_empty(vcoords)) cat(crayon::red(paste("Cooredinates for", filepaths$sampleids[i], "not able to be found...")))
    if(is_empty(vfeatures) | is_empty(vbarcodes) | is_empty(vcounts) | is_empty(vcoords)){
=======
    if(rlang::is_empty(vfeatures)) cat(crayon::red(paste("Features for", filepaths$sampleids[i], "not able to be found...")))
    if(rlang::is_empty(vbarcodes)) cat(crayon::red(paste("Barcodes for", filepaths$sampleids[i], "not able to be found...")))
    if(rlang::is_empty(vcounts)) cat(crayon::red(paste("Counts for", filepaths$sampleids[i], "not able to be found...")))
    if(rlang::is_empty(vcoords)) cat(crayon::red(paste("Coordinates for", filepaths$sampleids[i], "not able to be found...")))
    if(rlang::is_empty(vfeatures) | rlang::is_empty(vbarcodes) | rlang::is_empty(vcounts) | rlang::is_empty(vcoords)){
>>>>>>> oscar_dev
      fp_list[[i]] = list()
      missingSamples  = missingSamples + 1
    }
  }

  cat(crayon::green$bold(paste("Found", length(filepaths$sampleids)-missingSamples, "Visium Samples\n")))

  # Define number of available cores to use.
<<<<<<< HEAD
  cores = parallel::detectCores()
  if(cores > 1){
    cores = cores - 1
  }
  if(cores > length(filepaths[['count_found']])){
    cores = length(filepaths[['count_found']])
  }

  # Use parallelization to read count data if possible.
  cat(crayon::yellow(paste("Importing Samples....\n")))
  output_temp = mclapply(seq_along(filepaths[['count_found']]), function(i){
    if(length(fp_list[[i]]$counts) == 0){
      return(list())
    }
=======
  cores = count_cores(length(filepaths[['count_found']]))
  # cores = parallel::detectCores()
  # if(cores > 1){
  #   cores = cores - 1
  # }

  # Use parallelization to read count data if possible.
  output_temp = parallel::mclapply(seq_along(filepaths[['count_found']]), function(i){
    if(length(fp_list[[i]]$counts) == 0){
      return(list())
    }

>>>>>>> oscar_dev
    system(sprintf('echo "%s"', crayon::yellow(paste0("\tProcessing Sample ", i, "...."))))
    #cat(crayon::yellow(paste("\tSample", i, "\n")))

    # Process Visium folder.
    visium_processed = import_Visium(features_fp=fp_list[[i]][['features']],
                                     barcodes_fp=fp_list[[i]][['barcodes']],
                                     counts_fp=fp_list[[i]][['counts']],
                                     coords_fp=fp_list[[i]][['coords']])
    system(sprintf('echo "%s"', crayon::green(paste0("\tFinished Processing Sample ", i, "...."))))
    return(visium_processed)
<<<<<<< HEAD
  }, mc.cores=cores, mc.preschedule=T)
=======
  }, mc.cores=cores, mc.preschedule=F)
>>>>>>> oscar_dev
  cat(crayon::green$bold(paste("\tCompleted!\n")))

  # Organize the paralellized output into corresponding lists.
  return_lists = list()
  return_lists[['counts']] = list()
  return_lists[['coords']] = list()
  return_lists[['images']] = list()
  for(i in 1:length(output_temp)){
<<<<<<< HEAD
    if(is_empty(output_temp[[i]])){
      return_lists[['counts']][[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
      return_lists[['coords']][[fp_list[[i]]$runname]]  = output_temp[[i]]$coords
=======
    if(rlang::is_empty(output_temp[[i]])){
      return_lists[['counts']][[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
      return_lists[['coords']][[fp_list[[i]]$runname]] = output_temp[[i]]$coords
      return_lists[['images']][[fp_list[[i]]$runname]] = NULL
>>>>>>> oscar_dev
    }
    return_lists[['counts']][[fp_list[[i]]$runname]] = output_temp[[i]]$rawcounts
    return_lists[['coords']][[fp_list[[i]]$runname]] = output_temp[[i]]$coords
    return_lists[['images']][[fp_list[[i]]$runname]] = png::readPNG(fp_list[[i]]$image)
  }

  # Process duplicates genes if output comes from Visium
  dup_emsb = c()
  for(i in 1:length(names(return_lists[['counts']]))){
    dup_genes_Mask = duplicated(return_lists[['counts']][[i]]$gene)
    dup_emsb = append(dup_emsb, return_lists[['counts']][[i]]$emsb[dup_genes_Mask])
  }
  dup_emsb = unique(dup_emsb)
  for(i in 1:length(names(return_lists[['counts']]))){
    for(emsb in 1:length(dup_emsb)){
      geneToChange = return_lists[['counts']][[i]]$gene[return_lists[['counts']][[i]]$emsb == dup_emsb[emsb]]
      if(!is.null(geneToChange)){
        return_lists[['counts']][[i]]$gene[return_lists[['counts']][[i]]$emsb == dup_emsb[emsb]] = paste0(geneToChange, '_d', emsb)
      }
    }
    return_lists[['counts']][[i]] = return_lists[['counts']][[i]][, -1]
  }
  return(return_lists)
}

##
# process_lists: Takes two named lists of counts and coordinates and process gene and
# spot names before placing them within an STList
# @param rnacounts, a named list with counts
# @param spotcoords, a named list with coordinates
# @return return_lists a list with two lists within (one with counts, one with coordinates)
#
process_lists = function(counts_df_list, coords_df_list){
  require('dplyr')
  # Process the count and coordinate list.
  for(i in 1:length(names(counts_df_list))){

    # Get loop current name to be processed.
    name_i = names(counts_df_list)[i]

    # Column names of the count data are simplified using the clean_names() function.
    # NOTE: May need to reconsider use later.
    counts_df_list[[name_i]] = janitor::clean_names(counts_df_list[[name_i]])

    # Clean sample names in coordinates data using the clean_names(), so that they
    # mirror the column names in the count data frame.
    # NOTE: May need to reconsider use later.
    coords_df_list[[name_i]][, 1] = janitor::make_clean_names(unlist(coords_df_list[[name_i]][, 1]))

    # Test that spot names are the same in both count and coordinate data frames.
    if(length(setdiff(colnames(counts_df_list[[name_i]])[-1], unlist(coords_df_list[[name_i]][, 1]))) != 0){
      stop(paste0('The spots in the count  data (columns) and coordinate data (rows) do not match in spatial array (\"', name_i, '\").'))
    }

    array_col = names(coords_df_list[[name_i]])[3]
    # Sort coordinate data according to third column in the coordinate data frame.
    coords_df_list[[name_i]] = coords_df_list[[name_i]] %>%
      arrange(array_col)#[order(coords_df_list[[name_i]][, 3]), ]

    # Order column names in count data frame according to sorted coordinate data.
    counts_df_list[[name_i]] = counts_df_list[[name_i]][, c(colnames(counts_df_list[[name_i]][1]), unlist(coords_df_list[[name_i]][, 1]))]

    # Put column names to coordinate data.
    colnames(coords_df_list[[name_i]] ) = c('libname', 'ypos', 'xpos')
  }

<<<<<<< HEAD

=======
>>>>>>> oscar_dev
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
  sample_file_df = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F)
  sample_names = sample_file_df[[1]]
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(input_check$samples[1] == 'samplesfile_visium'){
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
# @params rnacounts, vector of file paths to Visium outputs or count matrices
# @params spotcoords, vector of file paths to coordinate matrices
# @param samples, vector of sample names
# @param input_check, a list resulting from detect_input
# @return filepaths, a list with found file paths per each sample, and sample IDs
#
process_sample_names = function(rnacounts, spotcoords, sample_names, input_check){
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(input_check$rna[1] == 'visium_out'){
    for(i in sample_names){
      sample_i = grep(i, rnacounts, value=T)
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
    }
    filepaths[['sampleids']] = sample_names
  } else{
    for(i in sample_names){
      sample_count = grep(i, rnacounts, value=T)
      sample_coord = grep(i, spotcoords, value=T)
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_count)
      filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_coord)
    }
    filepaths[['sampleids']] = sample_names
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
  sample_file_df = readr::read_delim(samples, delim=delsample, col_types=readr::cols(), progress=F)
  sample_names = sample_file_df[[1]]
  # Create objects to store sample names found. Will be used to assess if expected samples exist.
  # Also store file paths to open
  filepaths = list()
  if(input_check$rna[1] == 'visium_out'){
    for(i in sample_names){
      sample_i = grep(i, rnacounts, value=T)
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_i)
    }
    filepaths[['sampleids']] = sample_names
  } else{
    for(i in sample_names){
      sample_count = grep(i, rnacounts, value=T)
      sample_coord = grep(i, spotcoords, value=T)
      filepaths[['count_found']] = append(filepaths[['count_found']], sample_count)
      filepaths[['coord_found']] = append(filepaths[['coord_found']], sample_coord)
    }
    filepaths[['sampleids']] = sample_names
  }
  return(filepaths)
}

##
# process_meta: Takes the input file path of a metadata table and return the processed
# data frame for the STList
# @param samples, the file path to metadata file
# @param input_check, a list resulting from detect_input
# @param counts_df_list, a list with counts resulting from process_lists
# @return samples_df, a data frame to be placed within the clinical slot of the STList
#
process_meta = function(samples, input_check, counts_df_list){
  # Get delimiter of file from input_check
  del = input_check$samples[2]
  # Read file.
  samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F)
  # Get sorted list names and fetch corresponding rows from the samplefile.
  sorted_names = sort(names(counts_df_list))
  samples_df = samples_df[samples_df[[1]] %in% sorted_names, ]
  # Remove file paths from sample file if existent
  if(input_check$samples[1] == 'samplesfile_visium'){
    samples_df = samples_df[, -2]
  } else if(input_check$samples[1] == 'samplesfile_matrices'){
    samples_df = samples_df[, -c(2:3)]
  }
  # Test that number of rows in sample file correspond to number of elements in list.
  if(nrow(samples_df) != length(counts_df_list)){
    stop("The number of rows in the clinical data is not the same as the number of spatial arrays.")
  }
  return(samples_df)
}

process_meta_geomx = function(samples, input_check, counts_df_list, gmx_slide_col, gmx_meta_cols){
  if(input_check$samples[2] == 'xls'){
    samples_df = readxl::read_excel(samples)
  } else {
    # Get delimiter of file from input_check
    del = input_check$samples[2]
    samples_df = readr::read_delim(samples, delim=del, col_types=readr::cols(), progress=F)
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




