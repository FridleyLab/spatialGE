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
          function(x){
            cat("Spatial Transcriptomics List (STList)\n")
            cat(length(x@counts), "spatial arrays.\n")
#            cat((dim(x@counts)[2]-1), "sampled positions.\n")
#            cat(dim(x@counts)[1], "features/genes.")
            if(!is.null(x@clinical)){
              cat(paste0((ncol(x@clinical)-1),
                " variables in clinical data."))
            }
          }
)
