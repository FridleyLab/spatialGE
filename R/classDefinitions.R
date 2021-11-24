## Class definitions -----------------------------------------------------------

##
# Definition of an STList object class.
#
#
setClass(Class="STList",
         slots=list(counts="list",
                    coords="list",
                    clinical="tbl",
                    tr_counts="list",
                    gene_var="list",
                    gene_het="list",
                    gene_krige="list",
                    cell_deconv="list",
                    cell_krige="list",
                    st_clusters="list",
                    pheno_plots="list",
                    misc="list"
         )
)


# STList Methods ---------------------------------------------------------------

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

##
# @title summary: Prints overview of STList oject.
# @description Prints overview/summary of STList oject.
# @details
# This function takes an STList and prints a the number of spatial arrays in that
# object and other information about the object.
#
# @param object, an STList object to show summary from.
#
#
setMethod("summary", signature="STList",
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

##
# @title dim: Prints the dimensions of count arrays within an STList object.
# @description Returns the number of genes and spots for each array within an STList object
# @details
# This function takes an STList and prints the number of genes (rows) and spots (columns) of
# each spatial array within that object.
#
# @param object, an STList object to show summary from.
#
#
setMethod(dim, signature(x="STList"),
          function(x){
            dim_res = list()
            for(i in seq(x@counts)){
              dim_res[[i]] = c(base::nrow(x@counts[[i]]), base::ncol(x@counts[[i]]) - 1)
            }
            return(dim_res)
          }
)

