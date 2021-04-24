##
#' @title show
#' @description Prints summary of STList
#' @details
#' This function takes an STList and prints a the number of spatial arrays in that
#' object.
#'
#' @param object, an STList object to show summary from.
#' @export
#
#
# Set object class.
# setClass("STList", slots=list(counts="list",
#                               coords="list",
#                               clinical="tbl",
#                               voom_counts="list",
#                               gene_stdev="list",
#                               gene_het="list",
#                               gene_krige="list",
#                               cell_deconv="list",
#                               cell_stdev="list",
#                               cell_het="list",
#                               cell_krige="list",
#                               prediction_grid="list",
#                               prediction_border="list"
# ),
# )

setMethod("show", signature="STList",
          function(object){
            cat("Spatial Transcriptomics List (STList)\n")
            cat(length(object@counts), "spatial arrays.\n")
#            cat((dim(x@counts)[2]-1), "sampled positions.\n")
#            cat(dim(x@counts)[1], "features/genes.")
            if(!rlang::is_empty(object@clinical)){
              cat(paste0((ncol(object@clinical)-1),
                " variables in clinical data."))
            }
          }
)
