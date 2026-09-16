# Post-read processing: matches counts to coordinates, resolves sample names
# against file paths, and reads/validates sample metadata for STlist().

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
