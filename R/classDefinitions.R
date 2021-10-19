## Class definitions -----------------------------------------------------------

##
# Definition of an STList object class.
#
#
setClass(Class="STList",
         slots=list(counts="list",
                    coords="list",
                    clinical="tbl",
                    voom_counts="list",
                    log_counts="list",
                    gene_stdev="list",
                    log_stdev="list",
                    gene_het="list",
                    gene_krige="list",
                    cell_deconv="list",
                    cell_het="list",
                    cell_krige="list",
                    gene_krige_data="list",
                    deconv_krige_data="list",
                    st_clusters="list",
                    pheno_plots="list"
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

# setMethod("dim", signature="STList",
#           function(object){
#             if(length(object@counts == 1)){
#               dim_res = base::dim(object@counts[[1]])
#             } else{
#               dim_res = list()
#               for(i in 1:length(object@counts)){
#                 dim_res[[i]] = base::dim(object@counts[[i]])
#               }
#             }
#             return(dim_res)
#           },
# )
