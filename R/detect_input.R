##
# @title detect_input_rnacounts
# @description Detects the type of gene expression data being provided to the function STList.
# @details
# This function detects the input provided to the `rnacounts` argument of the STList()
# function. It also detects the delimiter of the file when relevant. NOTE that the
# function does minimum checking on the contents of files, limited mostly to detect
# csv/tsv, data frames, Visium/Xeniumss files, or Seurat objects.
#
# @param rnacounts the object, file, or directory with counts provided to STList.
# @return inputtype a string indicated the input type.
#
#
detect_input_rnacounts = function(rnacounts=NULL){
  # Define output/return variable.
  # Return NULL if no input was provided to rnacounts
  inputtype_rnacounts = c()
  if(is.null(rnacounts) || length(rnacounts) == 0){
    inputtype_rnacounts = NULL

    return(inputtype_rnacounts)
  }


  # CASE: SEURAT OBJECT
  # FUTURE DEV: May be good to be able to input a list of Seurat objects, instead of a single object at a time
  # ALSO CONSIDER: Currently multiple samples from Seurat are allowed as long as they are stored as "slices" (names(seurat_oj@images))
  # FUTURE DEV: Seurat V5 -_-
  if(is(rnacounts, 'Seurat')){
    inputtype_rnacounts = 'seurat'

    return(inputtype_rnacounts)
  }

  # CASE: NAMED LIST OF DATAFRAMES
  if(inherits(rnacounts, 'list')){
    # Test that lists are named.
    if(!is.null(names(rnacounts))){
      # Test that all elements in the list are data frames
      tests_df = all(sapply(rnacounts, is.data.frame))
      if(tests_df){
        inputtype_rnacounts = 'list_dfs'
      } else{
        raise_err(err_code='error0022')
      }
      rm(tests_df) # Clean env
    } else {
      raise_err(err_code='error0023')
    }

    return(inputtype_rnacounts)
  }


  # CASE: FILE PATH(S) TO COUNT MATRICES/TABLES (INCLUDES COSMX-SMI)
  # Test that all elements of the input vector exist as files and not as directories.
  test_files = all(sapply(rnacounts, file.exists))
  test_dirs = all(!sapply(rnacounts, dir.exists))
  if(test_files && test_dirs){
    # Test that files are CSV or TSV
    tests_text = lapply(rnacounts, readLines, n=2)
    tests_tab = unlist(lapply(tests_text, function(i){grepl("\t", i[2])}))
    tests_comma = unlist(lapply(tests_text, function(i){grepl(",", i[2])}))
    if(all(tests_tab)){
      inputtype_rnacounts = 'tab_delim'
    } else if(all(tests_comma)){
      inputtype_rnacounts = 'comma_delim'
    } else{
      raise_err(err_code='error0024')
    }
    rm(tests_tab, tests_comma) # Clean env

    # Check if inputs were COSMX-SMI
    tests_cosmx = unlist(lapply(tests_text, function(i){ ( grepl('fov', i[1]) && grepl('cell_ID|cell_id', i[1]) )}))
    if(all(tests_cosmx)){
      inputtype_rnacounts = paste0(inputtype_rnacounts, '_cosmx')
    }
    rm(tests_text, tests_cosmx) # Clean env

    return(inputtype_rnacounts)
  }
  rm(test_files, test_dirs) # Clean env


  # CASE: FILE PATHS TO VISIUM DIRECTORIES
  # Test that elements in `rnacounts` are directories.
  test_dirs = all(sapply(rnacounts, dir.exists))
  if(test_dirs){
    # Check that directories contains an element with name matching:
    # 'filtered_feature_bc' (visium), 'raw_feature_bc' (visium), or 'cell_feature_matrix' (xenium)
    visium_h5_filt = lapply(rnacounts, function(i){list.files(i, pattern='filtered_feature_bc[_a-zA-Z0-9]*\\.h5$', include.dirs=T, full.names=T)})
    visium_h5_raw = lapply(rnacounts, function(i){list.files(i, pattern='raw_feature_bc[_a-zA-Z0-9]*\\.h5$', include.dirs=T, full.names=T)})
    xenium_h5 = lapply(rnacounts, function(i){list.files(i, pattern='cell_feature_matrix[_a-zA-Z0-9]*\\.h5$', include.dirs=T, full.names=T)})
    visium_mex_filt = lapply(rnacounts, function(i){list.files(i, pattern='filtered_feature_bc[_a-zA-Z0-9]*\\.tar\\.gz$', include.dirs=T, full.names=T)})
    visium_mex_raw = lapply(rnacounts, function(i){list.files(i, pattern='raw_feature_bc[_a-zA-Z0-9]*\\.tar\\.gz$', include.dirs=T, full.names=T)})
    xenium_mex = lapply(rnacounts, function(i){list.files(i, pattern='cell_feature_matrix[_a-zA-Z0-9]*\\.tar\\.gz$', include.dirs=T, full.names=T)})

    # Are h5 or MEX files present?
    if(!(rlang::is_empty(unlist(visium_h5_filt)))){
      tests_files = lapply(visium_h5_filt, function(i){ hdf5r::is_hdf5(i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'visium_filtered_h5' }
    } else if(!(rlang::is_empty(unlist(visium_h5_raw)))){
      tests_files = lapply(visium_h5_raw, function(i){ hdf5r::is_hdf5(i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'visium_raw_h5' }
    } else if(!(rlang::is_empty(unlist(xenium_h5)))){
      tests_files = lapply(xenium_h5, function(i){ hdf5r::is_hdf5(i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'xenium_h5' }
    } else if(!(rlang::is_empty(unlist(visium_mex_filt)))){
      tests_files = lapply(visium_mex_filt, function(i){ grepl('filtered_feature_bc[_a-zA-Z0-9]*\\.tar\\.gz$', i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'visium_filtered_mex' }
    } else if(!(rlang::is_empty(unlist(visium_mex_raw)))){
      tests_files = lapply(visium_mex_raw, function(i){ grepl('raw_feature_bc[_a-zA-Z0-9]*\\.tar\\.gz$', i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'visium_raw_mex' }
    } else if(!(rlang::is_empty(unlist(xenium_mex)))){
      tests_files = lapply(xenium_mex, function(i){ grepl('cell_feature_matrix[_a-zA-Z0-9]*\\.tar\\.gz$', i[1]) })
      if(all(unlist(tests_files))){ inputtype_rnacounts = 'xenium_mex' }
    } else{
      raise_err(err_code='error0026')
    }

    return(inputtype_rnacounts)

  } else{
    raise_err(err_code='error0025')
  }
  rm(test_dirs) # Clean env

} # CLOSE detect_input_rnacounts


##
# @title detect_input_spotcoords
# @description Detects the type of coordinates file being provided to the function STList.
# @details
# This function detects the input provided to the `spotcoords` argument of the STList()
# function. It also detects the delimiter of the file when relevant. NOTE that the
# function does minimum checking on the contents of files, limited mostly to detect
# csv/tsv, data frames, or parquet files.
#
# @param spotcoords the object, file, or directory with counts provided to STList.
# @return inputtype a string indicated the input type.
#
#
detect_input_spotcoords = function(spotcoords=NULL){
  # Define output/return variable.
  # Return NULL if no input was provided to rnacounts
  inputtype_spotcoords = c()
  if(is.null(spotcoords) || length(spotcoords) == 0){
    inputtype_spotcoords = NULL

    return(inputtype_spotcoords)
  }

  # CASE: NAMED LIST OF DATAFRAMES (COORDINATES)
  if(inherits(spotcoords, 'list')){
    # Test that lists are named.
    if(!is.null(names(spotcoords))){
      # Test that all elements in the list are data frames
      tests_df = all(sapply(spotcoords, is.data.frame))
      if(tests_df){
        inputtype_spotcoords = 'list_dfs'
      } else{
        raise_err(err_code='error0030')
      }
      rm(tests_df) # Clean env
    } else {
      raise_err(err_code='error0031')
    }

    return(inputtype_spotcoords)
  }

  # CASE: FILE PATH(S) TO COORDINATES TABLES (INCLUDES COSMX-SMI)
  # Test that all elements of the input vector exist as files and not as directories.
  test_files = all(sapply(spotcoords, file.exists))
  test_dirs = all(!sapply(spotcoords, dir.exists))
  if(test_files && test_dirs){
    # Test that files are CSV or TSV
    tests_text = lapply(spotcoords, readLines, n=2)
    tests_tab = unlist(lapply(tests_text, function(i){grepl("\t", i[2])}))
    tests_comma = unlist(lapply(tests_text, function(i){grepl(",", i[2])}))
    if(all(tests_tab)){
      inputtype_spotcoords = 'tab_delim'
    } else if(all(tests_comma)){
      inputtype_spotcoords = 'comma_delim'
    } else{
      raise_err(err_code='error0029')
    }
    rm(tests_tab, tests_comma) # Clean env

    # Check if spotcoords inputs are COSMX-SMI
    tests_cosmx = unlist(lapply(tests_text, function(i){ ( grepl('fov', i[1]) && grepl('cell_ID|cell_id', i[1]) )}))
    if(all(tests_cosmx)){
      inputtype_spotcoords = paste0(inputtype_spotcoords, '_cosmx')
    }
    rm(tests_text, tests_cosmx) # Clean env

    return(inputtype_spotcoords)
  }
  rm(test_files, test_dirs) # Clean env

} # CLOSE detect_input_spotcoords


##
# @title detect_input_samples
# @description Detects the type of metadata (`samples`) being provided to the function STList.
# @details
# This function detects the input provided to the `rnacounts` argument of the
# STList() function. It also detects the delimiter of the file when relevant. NOTE
# that the function does minimum checking on the contents of files, limited mostly
# to detect csv/tsv, or a vector with strings.
#
# @param samples the object, file, or directory with counts provided to STList.
# @return inputtype a string indicated the input type.
#
#
detect_input_samples = function(samples=NULL){
  # Define output/return variable.
  # Return NULL if no input was provided to samples
  inputtype_samples = c()
  if(is.null(samples) || length(samples) == 0){
    inputtype_samples = NULL

    return(inputtype_samples)
  }

  # CASE: SAMPLES FILE
  # Test that input is a file path to a table (tsv, csv, or xlsx)
  if(length(samples) == 1 && file.exists(samples)){
    # Test if file path ends in 'xlsx'
    tests_xlsx = grepl("\\.xlsx$", samples)
    if(tests_xlsx){
      inputtype_samples = 'samplesfile_excel'
      rm(tests_xlsx) # Clean env
    } else{
      # Test if file is tab- or comma-delimited
      tests_text = readLines(samples, n=2)
      tests_tab = grepl("\t", tests_text)
      tests_comma = grepl(",", tests_text)

      if(all(tests_tab)){
        inputtype_samples = 'samplesfile_tab'
      } else if(all(tests_comma)){
        inputtype_samples = 'samplesfile_comma'
      } else{
        raise_err(err_code='error0027')
      }
      rm(tests_text, tests_tab, tests_comma) # Clean env
    }

    return(inputtype_samples)
  }

  # CASE: VECTOR WITH SAMPLE NAMES
  if(length(samples) >= 1){
    # Test that input is a vector with strings
    #test_vector = is.vector(samples, mode='character')
    inputtype_samples = 'sample_names'

    return(inputtype_samples)
  }

  # Test validity of sample names
  # if(test_vector){
  #   # Test if elements of the begin with a number (SHOULD NOT begin with a number)
  #   test_number = any(sapply(samples, function(i){grepl("^[0-9]", i)}))
  #   if(test_number){
  #     raise_err(err_code='error0028')
  #   }
  #   rm(test_number) # Clean env
  #
  #   # Test if elements contain only alpha-numerics, spaces, dash, and underscores
  #   test_chars = any(!sapply(samples, function(i){grepl("^[ //-_//A-Za-z0-9]+$", i)}))
  #   if(test_chars){
  #     raise_err(err_code='error0028')
  #   }
  #   rm(test_chars) # Clean env
  # }
  # rm(test_vector) # Clean env

} # CLOSE detect_input_samples


##
# @title detect_input
# @description Detects the type of data input being provided to the function STList.
# @details
# This function detects what input is being provided to the STList() function. It
# also detects the delimiter of the file when relevant. NOTE that the function does
# minimum checking on the contents of file, limited mostly to detect csv or tsv,
# Visium files, or Seurat objects. Checks are performed on the first element only,
# and thus other elements could not comply with the format.
#
# Dispatches across a fixed sequence of input-shape cases (see
# `detect_input_case_*` in detect_input_cases.R), threading the accumulated
# `inputtype` list through each. The two Seurat/GeoMx cases can short-circuit
# the sequence (matching the original function's early `return`s); the
# remaining cases only ever add to `inputtype` and always fall through to the
# next case.
#
# @param rnacounts the file/directory with counts provided to STList.
# @param spotcoords the file with coordinates  provided to STList.
# @param samples the metadata or sample names provided to STList.
# @return inputtype a list containing file types of input arguments.
#
#' @importFrom magrittr %>%
#
#
detect_input = function(rnacounts=NULL, spotcoords=NULL, samples=NULL){

  # Define output/return variable.
  # If variable remains NULL, then no valid input was given by the user.
  inputtype = list()
  inputtype$rna = NULL
  inputtype$coords = NULL
  inputtype$samples = NULL

  # CASE SEURAT OBJECT(S) WITH SAMPLE NAMES OR SAMPLE FILE
  case_result = detect_input_case_seurat(rnacounts, samples, inputtype)
  if(case_result$early_return){ return(case_result$inputtype) }
  inputtype = case_result$inputtype

  # CASE DCC FILES FROM GEOMX
  case_result = detect_input_case_geomx_dcc(rnacounts, samples, inputtype)
  if(case_result$early_return){ return(case_result$inputtype) }
  inputtype = case_result$inputtype
  rm(case_result) # Clean env

  # CASE ONLY SAMPLEFILE PROVIDED.
  inputtype = detect_input_case_samplefile_only(rnacounts, samples, inputtype)

  # CASE: NAMED LIST OF DATAFRAMES WITH OR WTHOUTH SAMPLEFILE, OR SAMPLE NAMES.
  inputtype = detect_input_case_list_dfs(rnacounts, spotcoords, samples, inputtype)

  # CASE: FILE PATH(S) TO COUNT AND COORDINATE MATRICES, AND SAMPLE NAMES (FILE OR VECTOR). COSMX-SMI INCLUDED
  inputtype = detect_input_case_matrices_with_samples(rnacounts, spotcoords, samples, inputtype)

  # CASE: FILE PATHS TO VISIUM DIRECTORIES.
  inputtype = detect_input_case_visium_dirs(rnacounts, spotcoords, samples, inputtype)

  return(inputtype)

} # CLOSE detect_input

