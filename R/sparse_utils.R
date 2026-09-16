# Sparse-matrix conversion helpers shared by STlist() and its readers.

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
