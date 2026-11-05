# Readers: import raw counts/coordinates for each supported input format into
# per-sample lists, dispatched by STlist().

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
read_seurat = function(rnacounts, assay=NULL){
  if(!requireNamespace('SeuratObject', quietly=TRUE)){
    raise_err(err_code='error0052')
  }

  # Assay name. Previously hardcoded to 'Spatial', so an object whose spatial
  # assay is named anything else -- 'SCT' after SCTransform, 'RNA' for several
  # pipelines, a custom name -- silently produced NULL counts.
  if(is.null(assay)){
    assay = SeuratObject::DefaultAssay(rnacounts)
  }
  if(!(assay %in% names(rnacounts@assays))){
    raise_err(err_code='error0053')
  }

  # Counts, sparse to sparse.
  #
  # v5 stores an Assay5 whose data lives in @layers, so the v4 access
  # rnacounts@assays$Spatial@counts returns NULL there -- the whole reason this
  # function did not work on v5. LayerData() is the version-agnostic accessor and
  # is available in SeuratObject >= 5.
  #
  # Also note what is NOT done here any more: the previous implementation called
  # as.data.frame() on the full matrix up front, densifying every gene x every
  # cell across ALL slices before subsetting to one. For a multi-slice Visium
  # object that is the largest allocation in the import path, and it happened
  # before a single slice was looked at.
  counts_mtx = tryCatch(
    SeuratObject::LayerData(rnacounts, assay=assay, layer='counts'),
    error=function(e) NULL)
  if(is.null(counts_mtx) || prod(dim(counts_mtx)) == 0){
    # v4 objects, and v5 objects whose counts layer is absent but data is present.
    counts_mtx = tryCatch(
      SeuratObject::GetAssayData(rnacounts, assay=assay, slot='counts'),
      error=function(e) NULL)
  }
  if(is.null(counts_mtx) || prod(dim(counts_mtx)) == 0){
    raise_err(err_code='error0054')
  }
  if(!inherits(counts_mtx, 'dgCMatrix')){
    counts_mtx = methods::as(methods::as(methods::as(counts_mtx, 'dMatrix'),
                                        'generalMatrix'), 'CsparseMatrix')
  }

  return_lists = list(counts=list(), coords=list(), images=list())

  slices = names(rnacounts@images)
  if(length(slices) == 0){
    raise_err(err_code='error0055')
  }

  for(i in slices){
    img = rnacounts@images[[i]]
    coords_df = seurat_image_coords(img)
    spots = intersect(coords_df[['libname']], colnames(counts_mtx))
    if(length(spots) == 0) next

    coords_df = coords_df[match(spots, coords_df[['libname']]), , drop=FALSE]
    return_lists[['coords']][[i]] = tibble::as_tibble(coords_df)

    # Subset first, densify once, and only for this slice's spots.
    return_lists[['counts']][[i]] = tibble::as_tibble(
      tibble::rownames_to_column(
        expandSparse(counts_mtx[, spots, drop=FALSE]), var='gene'))

    if(.hasSlot(img, 'image')) return_lists[['images']][[i]] = img@image
  }

  if(length(return_lists[['counts']]) == 0){
    raise_err(err_code='error0056')
  }
  return(return_lists)
}


##
# seurat_image_coords: Extracts spot/cell coordinates from one Seurat image
# object, across the layouts v4 and v5 use.
#
# v4 VisiumV1 carries a @coordinates data.frame with imagerow/imagecol. v5
# VisiumV2 and FOV instead hold a @boundaries list of Centroids/Segmentation
# objects, for which SeuratObject::GetTissueCoordinates() is the accessor. The
# returned column names differ between those paths ('x'/'y' vs
# 'imagerow'/'imagecol'), so they are normalized here rather than at the call
# site.
#
# @param img one element of a Seurat object's @images.
# @return a data.frame with libname, ypos, xpos.
#
seurat_image_coords = function(img){
  if(.hasSlot(img, 'coordinates') && is.data.frame(img@coordinates) &&
     all(c('imagerow', 'imagecol') %in% colnames(img@coordinates))){
    d = img@coordinates[, c('imagerow', 'imagecol'), drop=FALSE]
    return(data.frame(libname=rownames(d), ypos=d[['imagerow']], xpos=d[['imagecol']],
                      stringsAsFactors=FALSE))
  }

  tc = tryCatch(SeuratObject::GetTissueCoordinates(img), error=function(e) NULL)
  if(is.null(tc)){
    raise_err(err_code='error0057')
  }
  tc = as.data.frame(tc)

  # GetTissueCoordinates() returns x/y, and may carry the barcode either in a
  # 'cell' column (v5 Centroids) or in the row names.
  libname = if('cell' %in% colnames(tc)) as.character(tc[['cell']]) else rownames(tc)
  ycol = if('y' %in% colnames(tc)) 'y' else 'imagerow'
  xcol = if('x' %in% colnames(tc)) 'x' else 'imagecol'
  if(!all(c(ycol, xcol) %in% colnames(tc))){
    raise_err(err_code='error0057')
  }
  data.frame(libname=libname, ypos=tc[[ycol]], xpos=tc[[xcol]], stringsAsFactors=FALSE)
}

##
# read_one_counts_file: Reads a single delimited gene count file
# @param filepath the count file's path
# @param delim the file's delimiter
# @return a data frame with the file's raw gene counts
#
read_one_counts_file = function(filepath, delim){
  counts_df = data.table::fread(filepath, sep=delim)
  counts_df = as.data.frame(counts_df)
  return(counts_df)
}

##
# read_one_coords_file: Reads a single delimited spot/cell coordinate file
# @param filepath the coordinate file's path
# @param delim the file's delimiter
# @return a data frame with the file's spot/cell coordinates
#
read_one_coords_file = function(filepath, delim){
  coords_df = readr::read_delim(filepath, delim=delim, col_types=readr::cols(), progress=F)
  return(coords_df)
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
  counts_df_list = sge_lapply(seq_along(filepaths[['count_found']]), function(i){
    read_one_counts_file(filepaths[['count_found']][i], delrna)
  }, cores=cores, preschedule=F)
  # Name list elements.
  names(counts_df_list) = filepaths[['sampleids']]

  # Use parallelization to read coordinate data if possible.
  coords_df_list = sge_lapply(seq_along(filepaths[['coord_found']]), function(i){
    read_one_coords_file(filepaths[['coord_found']][i], delcoords)
  }, cores=cores, preschedule=F)
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
# read_one_visium_sample: Reads and imports a single Visium sample's counts
# and coordinates (MEX or h5)
# @param fp_entry the sample's entry from `fp_list` (features/barcodes/counts
# or h5 counts, plus coords)
# @param sample_index the sample's position, used only for progress messages
# @param input_rnacounts The result of detect_input_rnacounts
# @param verbose logical, whether to print text to console
# @return a list with the sample's raw counts and coordinates, or an empty
# list if the sample's counts were not found
#
read_one_visium_sample = function(fp_entry, sample_index, input_rnacounts, verbose){
  if(length(fp_entry$counts) == 0){
    return(list())
  }
  if(verbose){
    system(sprintf('echo "%s"', paste0("\t\tProcessing Sample ", sample_index, "....")))
  }

  # Process Visium outputs
  if(input_rnacounts %in% c('visium_filtered_mex', 'visium_raw_mex')){
    visium_processed = import_visium(features_fp=fp_entry[['features']],
                                     barcodes_fp=fp_entry[['barcodes']],
                                     counts_fp=fp_entry[['counts']],
                                     coords_fp=fp_entry[['coords']])
  } else{
    visium_processed = import_visium_h5(counts_fp=fp_entry[['counts']],
                                        coords_fp=fp_entry[['coords']])
  }

  if(verbose){
    system(sprintf('echo "%s"', paste0("\t\tFinished data read Sample ", sample_index)))
  }

  return(visium_processed)
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
  cores = resolve_cores(length(filepaths[['count_found']]), cores, coerce='as_integer')
  if(is.na(cores)){
    raise_err(err_code='error0036')
  }

  # Use parallelization to read count data if possible.
  output_temp = sge_lapply(seq_along(filepaths[['count_found']]), function(i){
    read_one_visium_sample(fp_list[[i]], i, input_rnacounts, verbose)
  }, cores=cores, preschedule=F)
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
# read_one_xenium_sample: Reads and imports a single Xenium sample's counts
# and coordinates (MEX or h5)
# @param fp_entry the sample's entry from `fp_list` (features/barcodes/counts
# or h5 counts, plus coords)
# @param sample_index the sample's position, used only for progress messages
# @param input_rnacounts The result of detect_input_rnacounts
# @param verbose logical, whether to print text to console
# @return a list with the sample's raw counts and coordinates, or an empty
# list if the sample's counts were not found
#
read_one_xenium_sample = function(fp_entry, sample_index, input_rnacounts, verbose){
  if(length(fp_entry$counts) == 0){
    return(list())
  }

  if(verbose){
    system(sprintf('echo "%s"', paste0("\t\tProcessing Sample ", sample_index, "....")))
  }

  # Process Visium outputs
  if(input_rnacounts[1] == 'xenium_mex'){
    xenium_processed = import_xenium_mex(features_fp=fp_entry[['features']],
                                     barcodes_fp=fp_entry[['barcodes']],
                                     counts_fp=fp_entry[['counts']],
                                     coords_fp=fp_entry[['coords']])
  } else{
    xenium_processed = import_xenium_h5(counts_fp=fp_entry[['counts']],
                                        coords_fp=fp_entry[['coords']])
  }

  if(verbose){
    system(sprintf('echo "%s"', paste0("\t\tFinished data read Sample ", sample_index)))
  }

  return(xenium_processed)
}

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

      # Test that all files have been found. vcoords is always a cells.parquet
      # path (matched above), never csv, so this must check for 'parquet' --
      # checking 'csv' made this fire unconditionally, masking any real result.
      needed_mex_test = c(!grepl('gz', vfeatures), !grepl('gz', vbarcodes), !grepl('gz', vcounts), !grepl('parquet', vcoords))
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
  cores = resolve_cores(length(filepaths[['count_found']]), cores, coerce='as_integer')
  if(is.na(cores)){
    raise_err(err_code='error0038')
  }

  # Use parallelization to read count data if possible.
  output_temp = sge_lapply(seq_along(filepaths[['count_found']]), function(i){
    read_one_xenium_sample(fp_list[[i]], i, input_rnacounts, verbose)
  }, cores=cores, preschedule=F)
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
# read_one_cosmx_slide: Reads and imports a single CosMx-SMI slide's counts
# and coordinates
# @param fp_entry the slide's entry from `fp_list` (counts, coords, runname)
# @param sample_index the slide's position, used only for progress messages
# @param verbose logical, whether to print text to console
# @return a list with the slide's raw counts and coordinates (one entry per
# FOV), or an empty list if the slide's counts were not found
#
read_one_cosmx_slide = function(fp_entry, sample_index, verbose){
  if(length(fp_entry[['counts']]) == 0){
    return(list())
  }

  if(verbose){
    system(sprintf('echo "%s"', paste0("\tProcessing sample ", sample_index, "....")))
  }

  # Process CosMx outputs.
  cosmx_processed = import_smi(counts_fp=fp_entry[['counts']],
                               coords_fp=fp_entry[['coords']],
                               slidename=fp_entry[['runname']])

  if(verbose){
    system(sprintf('echo "%s"', paste0("\tFinished data read sample ", sample_index)))
  }

  return(cosmx_processed)
}

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
  output_temp = sge_lapply(seq_along(1:length(fp_list)), function(i){
    read_one_cosmx_slide(fp_list[[i]], i, verbose)
  }, cores=cores, preschedule=F)

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
# seurat_platform: Infers the spatialGE platform label from a Seurat object's
# image classes.
#
# STlist() previously hardcoded 'visium' for any Seurat input, with a comment
# admitting it could not assume that. @misc$platform is not cosmetic: it selects
# STgradient's default nb_dist_thr (c(0.75, 1.25) for visium against c(0.25, 3)
# otherwise) and branches in plot_image() and STplot_interpolation(), so a wrong
# label changes numbers, not just metadata.
#
# @param obj a Seurat object.
# @return one of 'visium' or 'xenium'.
#
seurat_platform = function(obj){
  classes = vapply(obj@images, function(im) class(im)[1], character(1))
  if(length(classes) == 0) return('visium')
  # VisiumV1 (Seurat v4) and VisiumV2 (v5) are spot-based arrays. FOV is v5's
  # imaging-based container, used for Xenium and CosMx.
  if(any(grepl('^Visium', classes))) return('visium')
  if(any(classes == 'FOV' | grepl('Xenium|CosMx|Fov', classes))) return('xenium')
  'visium'
}


##
#' @title as_Seurat: Convert an STlist to a Seurat object
#' @description Converts an STlist into a Seurat (v5) object, one image/slice per
#' sample, so spatialGE results can be taken into Seurat-based workflows.
#' @details
#' Each sample in the STlist becomes one FOV in the returned object, named after
#' the sample. Raw counts go to the `counts` layer and, when `transform_data()`
#' has been run, transformed counts go to the `data` layer. Spot-level columns
#' from `@spatial_meta` -- including `STclust` cluster assignments -- are carried
#' into `@meta.data`, prefixed by sample where a merge would otherwise collide.
#'
#' Cell names are prefixed with the sample name (`<sample>_<libname>`) because
#' spot barcodes are only unique within a sample in the general case, and a merged
#' Seurat object requires globally unique cell names.
#'
#' @param x an STlist
#' @param assay the assay name to create. Default `'Spatial'`
#' @param samples the samples to include. If NULL (default), all samples
#' @return a Seurat object
#'
#' @export
#
as_Seurat = function(x=NULL, assay='Spatial', samples=NULL){
  assert_stlist(x)
  if(!requireNamespace('SeuratObject', quietly=TRUE)){
    raise_err(err_code='error0052')
  }
  samples = resolve_samples(x, samples, slot='counts', numeric_handling='na_omit',
                            char_transform='exact', post_check='any_in',
                            message='None of the requested samples are present in the STlist.')

  has_tr = !rlang::is_empty(x@tr_counts)

  # Prefix cell names with the sample. Spot barcodes are unique only within a
  # sample -- the fixture's own barcodes are, but Visium barcodes repeat across
  # sections -- and a merged Seurat object requires globally unique names.
  prefix_cols = function(m, s){
    colnames(m) = paste0(s, '_', colnames(m))
    m
  }

  counts_list = lapply(samples, function(s) prefix_cols(sample_counts(x, s), s))
  names(counts_list) = samples

  # Genes present in every requested sample. Seurat needs one matrix, so genes
  # absent from a sample cannot simply be dropped for that sample only; taking the
  # intersection is reported rather than silently applied.
  gene_sets = lapply(samples, function(s) sample_genes(x, s, 'counts'))
  common_genes = Reduce(intersect, gene_sets)
  if(length(common_genes) == 0){
    raise_err(err_code='error0058')
  }
  dropped = length(Reduce(union, gene_sets)) - length(common_genes)
  if(dropped > 0){
    message(paste0('as_Seurat: ', dropped,
                    ' gene(s) not shared by every requested sample were dropped.\n'))
  }

  merged_counts = do.call(cbind, lapply(counts_list, function(m) m[common_genes, , drop=FALSE]))

  obj = SeuratObject::CreateSeuratObject(counts=merged_counts, assay=assay)

  if(has_tr){
    tr_list = lapply(samples, function(s){
      m = prefix_cols(sample_tr_counts(x, s), s)
      # Transformed counts can hold fewer genes than raw -- sct drops low-expression
      # genes -- so align to the merged gene set, filling absent rows with zero.
      out = Matrix::Matrix(0, nrow=length(common_genes), ncol=ncol(m), sparse=TRUE,
                           dimnames=list(common_genes, colnames(m)))
      shared = intersect(common_genes, rownames(m))
      out[shared, ] = m[shared, , drop=FALSE]
      out
    })
    merged_tr = do.call(cbind, tr_list)
    SeuratObject::LayerData(obj, assay=assay, layer='data') = merged_tr
  }

  # Spot-level metadata.
  meta = do.call(rbind, lapply(samples, function(s){
    d = sample_spatial_meta(x, s)
    d = tibble::remove_rownames(as.data.frame(d))
    d[['sample_name']] = s
    rownames(d) = paste0(s, '_', d[['libname']])
    d
  }))
  meta = meta[colnames(obj), , drop=FALSE]
  obj = SeuratObject::AddMetaData(obj, metadata=meta)

  # One FOV per sample, so each keeps its own coordinate frame rather than being
  # overlaid in a single space.
  for(s in samples){
    d = sample_spatial_meta(x, s)
    cents = SeuratObject::CreateCentroids(data.frame(
      x = as.numeric(d[['xpos']]),
      y = as.numeric(d[['ypos']]),
      cell = paste0(s, '_', d[['libname']]),
      stringsAsFactors = FALSE))
    obj[[s]] = SeuratObject::CreateFOV(cents, type='centroids', assay=assay)
  }

  return(obj)
}
