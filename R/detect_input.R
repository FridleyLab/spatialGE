##
# @title detect_input: Determine what is being provided to STList
# @description Detects the type of input being provided to the fucntion STList.
# @details
# This function detects what input is being provided to the STList() function. It
# also detects the delimiter of the file when relevant. NOTE that the function does
# minimum checking on the contents of file, limited mostly to detect csv or tsv,
# Visium files, or Seurat objects. Checks are performed on the first element only,
# and thus other elements could not comply with the format.
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
  if(!is.null(rnacounts)){
    if(class(rnacounts) == 'Seurat'){
      inputtype$rna = 'seurat'
      inputtype$samples = 'samples_from_seurat'
      return(inputtype)
    }
  }

  # CASE DCC FILES FROM GEOMX
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
            return(inputtype)
          }
        }
      }
    }
  }

  # CASE ONLY SAMPLEFILE PROVIDED.
  # Test if `rnacounts` was not entered and `samples` argument was entered.
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

  # CASE: NAMED LIST OF DATAFRAMES WITH OR WTHOUTH SAMPLEFILE, OR SAMPLE NAMES.
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

  # CASE: FILE PATH(S) TO COUNT AND COORDINATE MATRICES, AND SAMPLE NAMES (FILE OR VECTOR). COSMX-SMI INCLUDED
  # Test that there is an input for both `rnacounts` and `spotcoords`.
  if(!is.null(rnacounts) && !is.null(spotcoords) && !is.null(samples)){
    # Test that the first (or only) element of input vector exist, and input is not list or a directory.
    if(!is.list(rnacounts[1]) || !is.list(spotcoords[1])){
      # Determine what was entered as `samples`.
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
            del = '\t'
          } else if(is_comma_rna){
            del = ','
          } else{
            stop('RNA counts file is not comma or tab-delimited')
          }
          # Determine delimiter of first `spotcoords` file.
          if(is_tab_coord){
            del = '\t'
          } else if(is_comma_coord){
            del = ','
          } else{
            stop('Coordinates file is not comma or tab-delimited')
          }

          # Check if COSMX-SMI was input
          if(grepl('fov', rna_file[1]) & grepl('cell_ID|cell_id', rna_file[1])){
            inputtype$rna = c('cosmx', del)
            inputtype$coords = c('cosmx', del)
          } else{
            inputtype$rna = c('rnapath', del)
            inputtype$coords = c('coordpath', del)
          }
        }
      }

      if(is.null(inputtype$rna)){
        stop('Could not open the RNA count file.')
      }

    }
  }

  # CASE: FILE PATHS TO VISIUM DIRECTORIES.
  # Test that `rnacounts` were provided and first element is a directory.
  # Need also sample names that partially match the file path to be provided
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
#    if(length(samples) == 1 && file.exists(samples)){
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
    } else{
      stop('Number of sample names do not match number of Visium output folders.')
    }
  }

  return(inputtype)

} # CLOSE ENTIRE FUNCTION
