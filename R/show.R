##
# This function takes an STList and prints a the number of spatial arrays in that
# object.
#
# @param object, an STList object to show summary of.
#
#
# Load 'tidyverse' for tibble manipulation.
# require('tidyverse')

setMethod("show", signature="STList",
          function(object){
            cat("Spatial Transcriptomics List (STList)\n")
            cat(length(object@counts), "spatial arrays.\n")
#            cat((dim(x@counts)[2]-1), "sampled positions.\n")
#            cat(dim(x@counts)[1], "features/genes.")
            if(!is_empty(object@clinical)){
              cat(paste0((ncol(object@clinical)-1),
                " variables in clinical data."))
            }
          }
)
