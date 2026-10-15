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
#'
#' @return Invisibly returns `NULL`. Called for the side effect of printing an
#' overview: the number of spatial arrays, their dimensions and a sample count.
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
#'
#' @return Invisibly returns `NULL`. Called for the side effect of printing an
#' overview: the number of spatial arrays, their dimensions and a sample count.
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
#' @title dim: Dimensions of the count array in an STList object.
#' @description Returns the number of genes (rows) and spots/cells (columns)
#' observed in an STList object.
#' @details
#' An `STlist` is a container for one or more spatial arrays, and different
#' samples may have different numbers of genes. `dim()` therefore has a single
#' well-defined answer only when the object holds exactly one sample; in that
#' case it returns `c(n_genes, n_spots)` as an integer vector, like any other
#' matrix-like object.
#'
#' For a multi-sample STlist `dim()` returns `NULL` (the base R convention for
#' "no dimensions"), because there is no single gene x spot shape to report. Use
#' [summarize_STlist()] for a per-sample table of dimensions, or [tissue_names()]
#' to iterate samples.
#'
#' Previously this method returned a *named list* of per-sample dimensions, which
#' violated the base `dim()` contract: it made [nrow()]/[ncol()] meaningless, and
#' `R CMD check` reported it as "dim() returns class 'list'".
#'
#' @param x an STList object.
#'
#' @return An integer vector `c(n_genes, n_spots)` when the object holds exactly
#' one sample; `NULL` for a multi-sample object, where there is no single gene x
#' spot shape to report.
#
setMethod(f="dim", signature="STlist",
          definition=function(x){
            # Contract: dim() must return an integer vector, or NULL.
            if(length(x@counts) != 1){
              return(NULL)
            }
            as.integer(sample_dim(x, names(x@counts)[1], 'counts'))
          }
)


##
# Validity checks for the STlist class.
#
# Ensures an STlist is internally consistent before it is used: that every
# sample-keyed slot refers to the same samples as @counts, that each count matrix
# has gene and spot names, and that the number of rows in @spatial_meta[[i]]
# matches the number of spots/cells (columns) in @counts[[i]].
#
setValidity("STlist", function(object){
  counts <- object@counts

  # An empty container is allowed (e.g. mid-construction); nothing to check.
  if(length(counts) == 0){
    return(TRUE)
  }

  sample_names <- names(counts)

  # Every slot keyed by sample must use the same keys as @counts.
  #
  # @gene_krige is deliberately NOT in this list: it is keyed by GENE first and
  # sample second (x@gene_krige[[gene]][[sample]]), so its top-level names are
  # gene symbols. Checking it here rejected every kriged object -- but only when
  # validity actually ran, and slot assignment via `@<-` does not trigger it, so
  # gene_interpolation() appeared to succeed and validObject() failed afterwards.
  for(slot_name in c("spatial_meta", "gene_meta", "tr_counts")){
    slot_value <- methods::slot(object, slot_name)
    if(length(slot_value) == 0){
      next
    }
    extra <- setdiff(names(slot_value), sample_names)
    if(length(extra) > 0){
      return(paste0("slot '", slot_name, "' has entries for samples not in @counts: ",
                    paste(extra, collapse=', '), "."))
    }
  }

  # @gene_krige: check the sample level, one gene at a time.
  if(length(object@gene_krige) > 0){
    for(g in names(object@gene_krige)){
      extra_samples <- setdiff(names(object@gene_krige[[g]]), sample_names)
      if(length(extra_samples) > 0){
        return(paste0("@gene_krige[['", g,
                      "']] has entries for samples not in @counts: ",
                      paste(extra_samples, collapse=', '), "."))
      }
    }
  }

  # Per-sample consistency between the count matrix and its metadata.
  #
  # Only entries that actually hold a matrix are checked: unit tests build
  # deliberately minimal stubs (e.g. counts = list(a = 1, b = 2)) to exercise
  # name-resolution helpers in isolation, and those must stay constructible.
  for(i in sample_names){
    cm <- counts[[i]]
    if(!(is.matrix(cm) || inherits(cm, "Matrix"))){
      next
    }
    if(is.null(rownames(cm)) || is.null(colnames(cm))){
      return(paste0("the count matrix for sample '", i,
                    "' must have both gene (row) and spot/cell (column) names."))
    }

    sm <- object@spatial_meta[[i]]
    if(is.null(sm)){
      return(paste0("no @spatial_meta entry for sample '", i, "'."))
    }
    if(!is.data.frame(sm)){
      return(paste0("@spatial_meta[['", i, "']] must be a data.frame."))
    }
    if(!("libname" %in% colnames(sm))){
      return(paste0("@spatial_meta[['", i,
                    "']] must have a 'libname' column identifying each spot/cell."))
    }
    if(nrow(sm) != ncol(cm)){
      return(paste0("sample '", i, "': @spatial_meta has ", nrow(sm),
                    " rows but the count matrix has ", ncol(cm), " spots/cells."))
    }
  }

  TRUE
})


##
#' @title length: Number of spatial arrays in an STList object.
#' @description Returns the number of samples (spatial arrays) held by an
#' STList object.
#' @details
#' An `STlist` is a container for one or more spatial arrays. `length()` reports
#' how many, which makes the object usable with the usual idioms --
#' `seq_along(x)`, `for(i in seq_len(length(x)))` -- instead of having to know
#' that the count lives in `length(x@counts)`.
#'
#' @param x an STList object.
#
#' @return A single integer: the number of spatial arrays in `x`.
#
setMethod("length", signature="STlist",
          definition=function(x){
            length(x@counts)
          }
)


##
#' @title [: Subset an STList object by sample.
#' @description Returns a new, valid STlist containing only the selected
#' spatial arrays (samples).
#' @details
#' `x[i]` selects samples, not genes or spots: an STlist holds several spatial
#' arrays and each has its own gene and spot sets, so there is no shared second
#' dimension to index. `i` may be a character vector of sample names (as returned
#' by [tissue_names()]), a numeric vector of positions, or a logical vector as
#' long as the number of samples.
#'
#' Every slot keyed by sample is subset alongside `@counts` -- `@spatial_meta`,
#' `@gene_meta`, `@tr_counts` and the per-sample entries inside each
#' `@gene_krige` element -- so the result is internally consistent and passes
#' validity. `@sample_meta` rows are matched by their `sample_name` column rather
#' than by position, so the metadata stays aligned with its sample.
#'
#' An unknown name is an error rather than a silent omission, matching the
#' package-wide rule for `samples=` arguments.
#'
#' @param x an STList object.
#' @param i sample names, positions or a logical vector selecting spatial arrays.
#' @param ... unused; present for compatibility with the `[` generic.
#' @param drop unused. Subsetting an STlist always returns an STlist, never a
#' bare matrix: dropping to the underlying counts would discard the metadata that
#' makes it an STlist.
#
#' @return An STlist containing only the selected spatial arrays. Subsetting to
#' zero samples is an error.
#
setMethod("[", signature(x="STlist", i="ANY", j="missing"),
          definition=function(x, i, ..., drop=TRUE){
            sample_names = names(x@counts)
            if(length(sample_names) == 0){
              raise_err(err_code='error0096')
            }

            if(is.character(i)){
              unknown = setdiff(i, sample_names)
              if(length(unknown) > 0){
                raise_err(err_code='error0095', samplename=NULL,
                          paste(unknown, collapse=', '),
                          paste(sample_names, collapse=', '))
              }
              sel = match(i, sample_names)
            } else if(is.logical(i)){
              # A short logical vector would silently recycle, selecting the
              # wrong samples, so require an exact length.
              if(length(i) != length(sample_names)){
                raise_err(err_code='error0096')
              }
              sel = which(i)
            } else if(is.numeric(i)){
              # Base R semantics, including negative selection. An out-of-range
              # positive index yields NA and is rejected rather than silently
              # returning a shorter object.
              sel = seq_along(sample_names)[i]
              if(anyNA(sel)){
                raise_err(err_code='error0096')
              }
            } else {
              raise_err(err_code='error0096')
            }

            if(length(sel) == 0){
              raise_err(err_code='error0096')
            }

            kept = sample_names[sel]

            # @sample_meta is one row per sample, but the column identifying that
            # sample is not fixed: STlist() names it 'sample_name' for unnamed or
            # list input, while a clinical CSV keeps its own header (the Thrane
            # data uses 'sampleID'). Match on the first column -- which is what
            # show()/summarize_STlist() already treat as the sample identifier --
            # and fall back to position, since STlist() builds the rows in
            # @counts order.
            sm = x@sample_meta
            if(length(sm) > 0 && nrow(sm) == length(sample_names)){
              sm_key = as.character(sm[[1]])
              if(all(kept %in% sm_key)){
                sm = sm[match(kept, sm_key), , drop=FALSE]
              } else {
                sm = sm[sel, , drop=FALSE]
              }
            }

            new("STlist",
                counts       = x@counts[kept],
                spatial_meta = x@spatial_meta[kept],
                gene_meta    = if(length(x@gene_meta) > 0) x@gene_meta[kept] else list(),
                sample_meta  = sm,
                tr_counts    = if(length(x@tr_counts) > 0) x@tr_counts[kept] else list(),
                # @gene_krige is keyed by GENE first, then sample.
                gene_krige   = if(length(x@gene_krige) > 0){
                                 lapply(x@gene_krige, function(g) g[intersect(names(g), kept)])
                               } else {
                                 list()
                               },
                misc         = x@misc)
          }
)


##
#' @title [[: Not supported for STList objects.
#' @description Raises an informative error. Subsetting an STlist is
#' one-dimensional and returns another STlist; use `x[i]`.
#
setMethod("[[", signature(x="STlist"),
          definition=function(x, i, ...){
            raise_err(err_code='error0096')
          }
)


##
#' @title [: Reject two-dimensional subsetting of an STList object.
#' @description Raises an informative error. An STlist holds several spatial
#' arrays, each with its own gene and spot sets, so there is no shared second
#' dimension to index; select samples with `x[i]` and then address genes or spots
#' through that sample's slots.
#
setMethod("[", signature(x="STlist", i="ANY", j="ANY"),
          definition=function(x, i, j, ..., drop=TRUE){
            raise_err(err_code='error0096')
          }
)

