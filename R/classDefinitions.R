## Class definitions -----------------------------------------------------------

##
#' Definition of an STlist object class.
#'
#' @slot counts per spot RNA counts
#' @slot spatial_meta per spot x,y coordinates
#' @slot gene_meta per gene statistics (e.g., average expression, variance, Moran's I)
#' @slot sample_meta dataframe with metadata per sample
#' @slot tr_counts transfromed per spot counts
#' @slot gene_krige results from kriging on gene expression
#' @slot misc Parameters and images from ST data
#
#
setClass(Class="STlist",
         slots=list(counts="list",
                    spatial_meta="list",
                    gene_meta="list",
                    sample_meta=class(tibble::tibble())[1],
                    tr_counts="list",
                    gene_krige="list",
                    misc="list"
         )
)


# STList Methods ---------------------------------------------------------------

##
#' @title show: Prints overview of STList oject.
#' @description Prints overview/summary of STList oject.
#' @details
#' This function takes an STList and prints a the number of spatial arrays in that
#' object and other information about the object.
#'
#' @param object an STList object to show summary from.
#
#
setMethod("show", signature="STlist",
          function(object){
            cat("Spatial Transcriptomics List (STlist).\n")
            cat(length(object@counts), "spatial array(s):\n")
            if(length(names(object@counts)) >= 11){
              counter = 1
              while(counter < 11){
                i = names(object@counts)[counter]
                cat(paste0('\t', i, " (", ncol(object@counts[[i]]), ' ROIs|spots|cells x ', nrow(object@counts[[i]]), ' genes)\n'))
                counter = counter + 1
              }
            } else{
              for(i in names(object@counts)){
                cat(paste0('\t', i, " (", ncol(object@counts[[i]]), ' ROIs|spots|cells x ', nrow(object@counts[[i]]), ' genes)\n'))
              }
            }
            if(length(names(object@counts)) > 10){
              cat(paste0('\tPlus ', (length(names(object@counts))-10), ' additional spatial array(s)\n'))
            }
            cat('\n')
            if(!rlang::is_empty(object@sample_meta)){
              cat(paste0((ncol(object@sample_meta)-1), " variables in sample-level data:\n"))
              cat('\t', paste0(colnames(object@sample_meta[, -1]), collapse = ', '), '\n')
            }
          }
)

##
#' @title summary: Prints overview of STList oject.
#' @description Prints overview/summary of STList oject.
#' @details
#' This function takes an STList and prints a the number of spatial arrays in that
#' object and other information about the object.
#'
#' @param object an STList object to show summary from.
#
#
setMethod("summary", signature="STlist",
          function(object){
            cat("Spatial Transcriptomics List (STlist).\n")
            cat(length(object@counts), "spatial array(s):\n")
            if(length(names(object@counts)) >= 11){
              counter = 1
              while(counter < 11){
                i = names(object@counts)[counter]
                cat(paste0('\t', i, " (", ncol(object@counts[[i]]), ' ROIs|spots|cells x ', nrow(object@counts[[i]]), ' genes)\n'))
                counter = counter + 1
              }
            } else{
              for(i in names(object@counts)){
                cat(paste0('\t', i, " (", ncol(object@counts[[i]]), ' ROIs|spots|cells x ', nrow(object@counts[[i]]), ' genes)\n'))
              }
            }
            if(length(names(object@counts)) > 10){
              cat(paste0('\tPlus ', (length(names(object@counts))-10), ' additional spatial array(s)\n'))
            }
            cat('\n')
            if(!rlang::is_empty(object@sample_meta)){
              cat(paste0((ncol(object@sample_meta)-1), " variables in sample-level data:\n"))
              cat('\t', paste0(colnames(object@sample_meta[, -1]), collapse = ', '), '\n')
            }
          }
)

##
#' @title dim: Prints the dimensions of count arrays within an STList object.
#' @description Returns the number of genes and spots for each array within an STList object
#' @details
#' This function takes an STList and prints the number of genes (rows) and spots (columns) of
#' each spatial array within that object.
#'
#' @param x an STList object to show summary from.
#
#
setMethod(f="dim", signature="STlist",
          definition=function(x){
            dim_res = list()
            for(i in seq(x@counts)){
              dim_res[[i]] = c(base::nrow(x@counts[[i]]), base::ncol(x@counts[[i]]))
            }
            names(dim_res) = names(x@counts)
            return(dim_res)
          }
)

