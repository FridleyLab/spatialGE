# Case handlers for detect_input()'s dispatch over rnacounts/spotcoords/samples
# combinations. Each function runs one of detect_input()'s original CASE blocks
# verbatim; the two that used to `return(inputtype)` early signal that via
# `early_return` so detect_input() can still short-circuit in the same places.

##
# detect_input_case_seurat
# CASE: Seurat object(s) with sample names or sample file.
# @param rnacounts the object, file, or directory with counts provided to STList.
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return a list with `inputtype` and `early_return` (whether detect_input()
# should return `inputtype` immediately).
#
detect_input_case_seurat = function(rnacounts, samples, inputtype){
  if(!is.null(rnacounts)){
    if(is(rnacounts, 'Seurat')){
      inputtype$rna = 'seurat'
      inputtype$samples = 'samples_from_seurat'
      return(list(inputtype=inputtype, early_return=TRUE))
    }
  }
  return(list(inputtype=inputtype, early_return=FALSE))
}


##
# detect_input_case_geomx_dcc
# CASE: DCC files from GeoMx.
# @param rnacounts the object, file, or directory with counts provided to STList.
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return a list with `inputtype` and `early_return` (whether detect_input()
# should return `inputtype` immediately).
#
detect_input_case_geomx_dcc = function(rnacounts, samples, inputtype){
  # To prevent NOTES in R CMD check
  . = NULL

  if(!is.null(rnacounts) && !is.null(samples)){
    if(is.character(rnacounts)){
      if(dir.exists(rnacounts[1])){
        dcc_files = list.files(rnacounts, full.names=T, pattern='.dcc$', recursive=T)
        if(!is.null(dcc_files)){
          if(length(dcc_files) != 0){
            test_dcc = readLines(dcc_files[1]) %>% grep('<Code_Summary>', .)
            if(length(test_dcc) != 0){
              inputtype$rna = 'geomx_dcc'

              # Read metadata file and get coordinate information
              if(grepl('.xls', samples)){
                inputtype$samples = c('samplesfile_geomx', 'xls')
              } else{
                samples_file = readLines(samples, n=2)
                is_tab_samples = grepl("\t", samples_file[2])
                is_comma_samples = grepl(",", samples_file[2])
                # Determine delimiter of file.
                if(is_tab_samples){
                  del = '\t'
                } else if(is_comma_samples){
                  del = ','
                } else{
                  stop('Samples file is not comma, tab-delimited, or .xls file')
                }
                inputtype$samples = c('samplesfile_geomx', del)
              }
            }
            return(list(inputtype=inputtype, early_return=TRUE))
          }
        }
      }
    }
  }
  return(list(inputtype=inputtype, early_return=FALSE))
}


##
# detect_input_case_samplefile_only
# CASE: only samplefile provided (`rnacounts` not entered, `samples` entered).
# @param rnacounts the object, file, or directory with counts provided to STList.
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return the (possibly updated) `inputtype` list.
#
detect_input_case_samplefile_only = function(rnacounts, samples, inputtype){
  if(is.null(rnacounts) && !is.null(samples)){
    # test that `samples` is a single string and a file path.
    if(length(samples) == 1 && file.exists(samples)){
      # Read two first lines of file and see if it's csv or tsv.
      samples_file = readLines(samples, n=2)
      is_tab_samples = grepl("\t", samples_file[2])
      is_comma_samples = grepl(",", samples_file[2])
      # Determine delimiter of file.
      if(is_tab_samples){
        del = '\t'
      } else if(is_comma_samples){
        del = ','
      } else{
        stop('Samples file is not comma or tab-delimited')
      }
      samples_file_path_test = unlist(strsplit(samples_file[2], del))
      # See if second column is a file path.
      if(file.exists(samples_file_path_test[2]) && file.exists(samples_file_path_test[3])){
        inputtype$samples = c('samplesfile_matrices', del)
      } else if(dir.exists(samples_file_path_test[2]) && !dir.exists(samples_file_path_test[3])){

        if(dir.exists(samples_file_path_test[2])){
          # Check that dirctory contains an element with name matching 'filtered_feature_bc'.
          visium_check = list.files(samples_file_path_test[2], pattern='[raw|filtered]_feature_bc', include.dirs=T, full.names=T)
          if(!(rlang::is_empty(visium_check))){
            h5_test = grep('\\.h5$', visium_check, value=T)
            if(!(rlang::is_empty(h5_test))){
              if(hdf5r::is_hdf5(h5_test)){
                inputtype$samples = c('samplesfile_visium_h5', del)
              } else{
                warning('The .h5 file does not seem to be in HDF5 format')
              }
            } else{
              inputtype$samples = c('samplesfile_visium_mex', del)
            }
          }
        } else{
          stop('If intended input is a Visium output, could not find directory path.')
        }
      } else(
        stop('Samples file does not contain file paths or format is not compatible.')
      )
    }
  }
  return(inputtype)
}


##
# detect_input_case_list_dfs
# CASE: named list of dataframes, with or without samplefile/sample names.
# @param rnacounts the object, file, or directory with counts provided to STList.
# @param spotcoords the file with coordinates provided to STList.
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return the (possibly updated) `inputtype` list.
#
detect_input_case_list_dfs = function(rnacounts, spotcoords, samples, inputtype){
  # Test if inputs were entered for both 'rnacounts' and 'spotcoords'.
  if(!is.null(rnacounts) && !is.null(spotcoords)){
    # Test that `rnacounts` and `spotcoords` are lists.
    if(inherits(rnacounts, 'list') && inherits(spotcoords, 'list')){
      # Test that lists are named.
      if(!is.null(names(rnacounts)) && !is.null(names(spotcoords))){
        inputtype$rna = 'list_dfs'
        inputtype$coords = 'list_dfs'
        # Test if samples file (metadata) was provided.
        if(!is.null(samples)){
          if(length(samples) == 1 && file.exists(samples)){
            # Read samples file and see which delimiter has.
            samples_file = readLines(samples, n=2)
            is_tab_samples = grepl("\t", samples_file[2])
            is_comma_samples = grepl(",", samples_file[2])
            # Determine delimiter of file.
            if(is_tab_samples){
              del = '\t'
            } else if(is_comma_samples){
              del = ','
            } else{
              stop('Samples file is not comma or tab-delimited')
            }
            inputtype$samples = c('samplesfile', del)
          } else{
            inputtype$samples = 'names_from_list_or_df'
          }
        } else{
          inputtype$samples = 'names_from_list_or_df'
        }
      } else(
        raise_err(err_code='error0003')
      )
    }
  }
  return(inputtype)
}


##
# detect_input_case_matrices_with_samples
# CASE: file path(s) to count and coordinate matrices, and sample names
# (file or vector). CosMx-SMI included.
# @param rnacounts the file/directory with counts provided to STList.
# @param spotcoords the file with coordinates provided to STList.
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return the (possibly updated) `inputtype` list.
#
detect_input_case_matrices_with_samples = function(rnacounts, spotcoords, samples, inputtype){
  # Test that there is an input for both `rnacounts` and `spotcoords`.
  if(!is.null(rnacounts) && !is.null(spotcoords) && !is.null(samples)){
    # Test that the first (or only) element of input vector exist, and input is not list or a directory.
    if(!is.list(rnacounts[1]) || !is.list(spotcoords[1])){
      # Determine what was entered as `samples`
      # i.e., string with path to clinical file or sample names vector, and NOT a directory with a name matching the name of the sample (when single sample entered)
      if(length(samples) == 1 && file.exists(samples) && !dir.exists(samples)){
        # Read samples file and see which delimiter has.
        samples_file = readLines(samples, n=2)
        is_tab_samples = grepl("\t", samples_file[2])
        is_comma_samples = grepl(",", samples_file[2])
        # Determine delimiter of file.
        if(is_tab_samples){
          del = '\t'
        } else if(is_comma_samples){
          del = ','
        } else{
          stop('Samples file is not comma or tab-delimited')
        }
        inputtype$samples = c('samplesfile', del)
      } else if(length(samples) == length(rnacounts)){ # Check that sample names were entered instead samplefile.
        inputtype$samples = 'sample_names'
      } else{
        stop('Number of sample names do not match number of RNA counts tables.')
      }


      if(file.exists(rnacounts[1]) && file.exists(spotcoords[1]) && !dir.exists((rnacounts[1]))){
        if(length(rnacounts) == length(spotcoords)){
          # Read first `rnacounts` and `spotcoords` to detect delimiters.
          rna_file = readLines(rnacounts[1], n=2)
          coord_file = readLines(spotcoords[1], n=2)
          is_tab_rna = grepl("\t", rna_file[2])
          is_comma_rna = grepl(",", rna_file[2])
          is_tab_coord = grepl("\t", coord_file[2])
          is_comma_coord = grepl(",", coord_file[2])
          # Determine delimiter of first `rnacounts` file.
          if(is_tab_rna){
            del_rna = '\t'
          } else if(is_comma_rna){
            del_rna = ','
          } else{
            stop('RNA counts file is not comma or tab-delimited')
          }
          # Determine delimiter of first `spotcoords` file.
          if(is_tab_coord){
            del_coords = '\t'
          } else if(is_comma_coord){
            del_coords = ','
          } else{
            stop('Coordinates file is not comma or tab-delimited')
          }

          # Check if COSMX-SMI was input
          if(grepl('fov', rna_file[1]) & grepl('cell_ID|cell_id', rna_file[1])){
            inputtype$rna = c('cosmx', del_rna)
            inputtype$coords = c('cosmx', del_coords)
          } else{
            inputtype$rna = c('rnapath', del_rna)
            inputtype$coords = c('coordpath', del_coords)
          }
        }
      }

      if(is.null(inputtype$rna)){
        stop('Could not open the RNA count file.')
      }

    }
  }
  return(inputtype)
}


##
# detect_input_case_visium_dirs
# CASE: file paths to Visium directories. Needs sample names that partially
# match the file path to be provided.
# @param rnacounts the directory with counts provided to STList.
# @param spotcoords the file with coordinates provided to STList (must be NULL
# for this case to apply).
# @param samples the metadata or sample names provided to STList.
# @param inputtype the inputtype list accumulated so far.
# @return the (possibly updated) `inputtype` list.
#
detect_input_case_visium_dirs = function(rnacounts, spotcoords, samples, inputtype){
  # Test that `rnacounts` were provided and first element is a directory.
  if(!is.null(rnacounts) && is.null(spotcoords) && !is.null(samples)){
    if(dir.exists(rnacounts[1])){
      # Check that dirctory contains an element with name matching 'filtered_feature_bc'.
      visium_check = list.files(rnacounts[1], pattern='[raw|filtered]_feature_bc', include.dirs=T, full.names=T)
      if(!(rlang::is_empty(visium_check))){
        h5_test = grep('\\.h5$', visium_check, value=T)
        if(!(rlang::is_empty(h5_test))){
          if(hdf5r::is_hdf5(h5_test)){
            inputtype$rna = 'visium_out_h5'
          } else{
            warning('The .h5 file does not seem to be in HDF5 format')
          }
        } else{
          inputtype$rna = 'visium_out_mex'
        }
      }
    } else{
      stop('If intended input is a Visium output, could not find directory path.')
    }

    # Determine what was entered as `samples`.
    if(length(samples) == 1 && file.exists(samples) && !dir.exists(samples)){ # Suggested by Mr. Manjarres
      # Read samples file and see which delimiter has.
      samples_file = readLines(samples, n=2)
      is_tab_samples = grepl("\t", samples_file[2])
      is_comma_samples = grepl(",", samples_file[2])
      # Determine delimiter of file.
      if(is_tab_samples){
        del = '\t'
      } else if(is_comma_samples){
        del = ','
      } else{
        stop('Samples file is not comma or tab-delimited')
      }
      inputtype$samples = c('samplesfile', del)
    } else if(length(samples) == length(rnacounts)){
      inputtype$samples = 'sample_names'
    } else if(is.data.frame(samples)){
      raise_err(err_code='error0004')
    } else{
      stop('Number of sample names do not match number of Visium output folders.')
    }
  }
  return(inputtype)
}
