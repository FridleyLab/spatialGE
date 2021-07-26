##
# @title show: Prints overview of STList oject.
# @description Prints overview/summary of STList oject.
# @details
# This function takes an STList and prints a the number of spatial arrays in that
# object and other information about the object.
#
# @param object, an STList object to show summary from.
#
#
setMethod("show", signature="STList",
          function(object){
            cat("Spatial Transcriptomics List (STList).\n")
            cat(length(object@counts), "spatial array(s):\n")
            cat(paste0('\t', names(object@counts), '\n'))
            cat('\n')
            if(!rlang::is_empty(object@clinical)){
              cat(paste0((ncol(object@clinical)-1), " variables in sample data:\n"))
              cat('\t', paste0(colnames(object@clinical[, -1]), collapse = ', '))
            }
          }
)
