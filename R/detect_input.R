##
# @title detect_input: Determine what is being provided to STList
# @description Detects the input type to the fucntion STList
# @details
# This function detects if a list of data frames, a vector with filepaths, or a
# file with file paths is being provided to the list.
#
# @param countfiles, the counts input provided to STList.
# @param coordfiles, the ccoordinates input provided to STList.
# @return inputtype, a string defining what type of input was provided by the user.
#
#
detect_input <- function(countfiles=NULL, coordfiles=NULL){

  inputtype <- NULL

  # Test that count data and coordinate data was provided.
  if(!is.null(countfiles) & !is.null(coordfiles)){

    # Test whether both count and coordinate data are lists of data frames.
    if(inherits(countfiles, 'list') && inherits(coordfiles, 'list')){

      inputtype<- 'is_list'
      return(inputtype)

      # Test if input is a dataframe count and a dataframe for coordinates.
    } else if(is.data.frame(countfiles) && is.data.frame(coordfiles)){

      inputtype<- 'is_df'
      return(inputtype)

      # Test if filepaths was provided.
    } else if(is.vector(countfiles) & is.vector(coordfiles)){

      if(length(countfiles) != 1 && length(coordfiles) != 1){
        for(i in 1:length(countfiles)){
          # Test that a pair of count and coord files exist.
          if(!(file.exists(countfiles[i]) && file.exists(coordfiles[i]))){
            stop("Either one of the count files or coordinate files do not exists.")
          }
        }

        # Test if file path is a single array or several arrays.
        count_test <- readLines(countfiles[1]) #, what=character(), quiet=T)
        coord_test <- readLines(coordfiles[1]) #, what=character(), quiet=T)

        is_tab_count <- grepl("\t", count_test[1])
        is_tab_coord <- grepl("\t", coord_test[1])
        is_comma_count <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", count_test[1])
        is_comma_coord <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", coord_test[1])

        if(is_tab_count && is_tab_coord){
          inputtype <- 'is_vector_tsv'
          return(inputtype)
        }

        if(is_comma_count && is_comma_coord){
          inputtype <- 'is_vector_csv'
          return(inputtype)
        } else{
          stop('Both count and coordinate files must have the same delimiter.')
        }

        return(inputtype)

      } else if(length(countfiles) == 1 && length(coordfiles) == 1){
        if(file.exists(countfiles) && file.exists(coordfiles)){

          # Test if file path is a single array or several arrays.
          count_test <- readLines(countfiles)
          coord_test <- readLines(coordfiles)

          is_tab_count <- grepl("\t", count_test[1])
          is_tab_coord <- grepl("\t", coord_test[1])
          is_comma_count <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", count_test[1])
          is_comma_coord <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", coord_test[1])

          if(is_tab_count && is_tab_coord){
            inputtype <- 'is_tsv'
            return(inputtype)
          } else if (is_comma_count && is_comma_coord){
            inputtype <- 'is_csv'
            return(inputtype)
          } else{
            # Get count and coord filepaths and test there is an equal number of both.
            if(length(count_test) != length(coord_test)){
              stop("The number of count files is different to the number of coordinate files.")
            }
            count_test_single <- readLines(count_test[1])
            coord_test_single <- readLines(coord_test[1])

            is_tab_count <- grepl("\t", count_test_single[1])
            is_tab_coord <- grepl("\t", coord_test_single[1])
            is_comma_count <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", count_test_single[1])
            is_comma_coord <- grepl("^[ A-Za-z0-9]+,[ A-Za-z0-9]+", count_test_single[1])

            if(is_tab_count && is_tab_coord){
              inputtype <- 'is_several_tsv'
              return(inputtype)
            } else if (is_comma_count && is_comma_coord){
              inputtype <- 'is_several_csv'
              return(inputtype)
            } else{
              stop('Both count and coordinate files must have the same delimiter.')
            }

          }

        } else{
          stop('Could not find either the counts file or coordinates file.')
        }

      } else{
        stop('The count and coordinate inputs do not have the same number of elements.')
      }

    } else{
      stop('Please provide one of these valid input options:\n
      A. One data frame for counts + one data frame for coordinates.\n
      B. One path to a file containing counts + one path to a file containing coordinates.\n
      C. One path to a file with file paths to N count files + one path to a file containing N coordinate files.\n
      D. One list containing count data frames + one list containing coordinate data frames.')
    }

  } else{
    stop('Please provide both count and coordinate data.')
  }

}
