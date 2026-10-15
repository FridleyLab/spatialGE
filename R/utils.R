##
# Utility and helper functions

##
# Get the number of cores to use in parallel as a function of the number of
# data to process. It will not yield a higher number of cores than half of the
# total available cores. Will default to 1 core if not Unix
# @param n an integer representing the number of units to process
# @return the number of cores to use
#
# @importFrom methods as is new
#
#
count_cores = function(n){
  cores = 1
  if(.Platform$OS.type == 'unix'){
    # Use parallelization (if possible) to read data.
    avail_cores = (parallel::detectCores()) / 2
    if(avail_cores <= n){
      cores = avail_cores
    } else{
      cores = ifelse(n > 0, n, 1)
    }
  }
  # mclapply's internal .check_ncores() stops (or warns) when more than 2
  # processes are requested and _R_CHECK_LIMIT_CORES_ is set (R CMD check
  # sets this) -- clamp the auto-detected count so it never trips that check
  # on a many-core machine, independent of how many units of work there are.
  limit_cores = tolower(Sys.getenv('_R_CHECK_LIMIT_CORES_', ''))
  if(nzchar(limit_cores) && limit_cores != 'false' && cores > 2){
    cores = 2
  }
  return(cores)
}


##
#' @title spatial_metadata: Prints the names of the available spot/cell annotations
#' @description returns the names of the annotations in the `x@spatial_meta` slot.
#'
#' @param x an STList object
#' @return a list of character vectors containing the column names of `x@spatial_meta`
#'
#' @export
#
#
spatial_metadata = function(x){

  # Validate, and accept a Seurat object by coercing it (see assert_stlist).
  x = assert_stlist(x)
  meta_res = list()
  for(i in seq(x@spatial_meta)){
    meta_res[[i]] = base::colnames(x@spatial_meta[[i]][, -c(1:3)])
  }
  names(meta_res) = names(x@spatial_meta)

  return(meta_res)
}


##
#' @title tissue_names: Prints the names of the tissue samples in the STlist
#' @description returns a character vector with the names of tissue samples in the
#' STlist.
#'
#' @param x an STList object
#' @return a character vector with the sample names in the STlist object
#'
#' @export
#
#
tissue_names = function(x){

  # Validate, and accept a Seurat object by coercing it (see assert_stlist).
  x = assert_stlist(x)
  sample_res = names(x@counts)

  return(sample_res)
}


##
# @title color_parse: Creates a color palette
# @description Uses Khroma or RColorBrewer to return the colors of a palette name.
# @details
# This function takes a character string and uses either khroma or RColoBrewer to
# create a color palette. The function first looks in khroma, then RColoBrewer. In
# other words, khroma colors have priority.
#
# @param color_pal A name of a Khroma or RColorBrewer. Alternatively, a vector with
# colors of lenght equal or larger than number of categories.
# @param n_cats The number of colors to produce. If NULL, assumes 5 colors.
# @return cat_cols A color palette
#
#' @importFrom methods as is new
#' @importFrom grDevices colorRampPalette
#
# .hasSlot is declared here rather than beside its two uses in
# R/stlist_readers.R because roxygen2 >= 8 drops an @importFrom whose symbol it
# cannot see being called -- and those are conditional calls on a Seurat object,
# which is in Suggests. Without the import, reading a Seurat image slot fails on
# a minimal install. A CI step greps NAMESPACE for it, so a future document()
# run cannot silently remove it again.
#' @importFrom methods .hasSlot
#
#
color_parse = function(color_pal=NULL, n_cats=NULL){
  # Get color palette and number of colors needed.
  # Get names of Khroma colors.
  khroma_cols = khroma::info()
  khroma_cols = khroma_cols$palette

  # Assume 5 categories if n_cats not provided (for kriging/quilts).
  if(is.null(n_cats)){
    n_cats = 5
  }

  # Test if input is a Khroma name or RColorBrewer.
  # If so, create palette.
  if(color_pal[1] %in% khroma_cols){
    cat_cols = as.vector(khroma::colour(color_pal[1], force=T)(n_cats))
  }else if(color_pal[1] %in% rownames(RColorBrewer::brewer.pal.info)){
    cat_cols = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n_cats, color_pal[1]))(n_cats)
  }else{ # Test if user provided a vector of colors.
    if(length(color_pal) >= n_cats){
      cat_cols = color_pal[1:n_cats]
    } else{
      stop('Provide enough colors to plot or an appropriate Khroma/RColorBrewer palette.')
    }
  }

  # Subset colors is more colors in palette than categories
  if(length(cat_cols) > n_cats){
    cat_cols = cat_cols[1:n_cats]
  }
  return(cat_cols)
}


##
#' @title load_images: Place tissue images within STlist
#' @description Loads the images from tissues to the appropriate STlist slot.
#' @details
#' This function looks for `.PNG` or `.JPG` files within a folder matching the
#' sample names in an existing STlist. Then, loads the images to the STlist which
#' can be used for plotting along with other spatialGE plots.
#'
#' @param x an STlist
#' @param images image file paths or one or more folders to load images from.
#'   A named character vector (`c(sample_name = "path/to/image.png")`) maps
#'   samples to files explicitly and bypasses name matching. An unnamed vector
#'   is matched against the sample names in `x`; any element that is an existing
#'   directory is expanded to the image files it contains. Image filenames are
#'   matched to sample names by looking for `<sample_name>.` in the filename and,
#'   failing that, a normalized prefix match (case- and separator-insensitive).
#' @return an STlist with images
#'
#' @export
#'
#' @importFrom methods as is new
#
#
load_images = function(x=NULL, images=NULL){
  # Test if an STList has been input.
  x = assert_stlist(x)

  if(is.null(images)){
    stop("Please, provide a vector with images file paths.")
  }
  if(!is.character(images) || length(images) == 0){
    stop("`images` must be a character vector of image file paths or folders.")
  }

  # A named vector maps samples to files directly, which is the escape hatch when
  # filenames do not resemble the sample names at all (e.g. Visium exports that
  # name images after a GEO accession rather than the sample folder).
  explicit_map = NULL
  if(!is.null(names(images)) && any(nzchar(names(images)))){
    explicit_map = images[nzchar(names(images))]
  }
  paths = unname(images)

  # Expand any directories (one or several) into the image files they contain.
  is_dir = dir.exists(paths)
  if(any(is_dir)){
    from_dirs = unlist(lapply(paths[is_dir], list.files, full.names=TRUE), use.names=FALSE)
    paths = c(paths[!is_dir], from_dirs)
  }

  # Keep only existing files with an image extension.
  keep = file.exists(paths) & grepl('\\.(png|jpg|jpeg)$', paths, ignore.case=TRUE)
  if(!any(keep)){
    stop("None of the supplied paths are existing .png, .jpg or .jpeg image files.")
  }
  images = unique(paths[keep])

  # Tolerant key for fallback matching: lowercase, separators removed.
  norm_key = function(z) gsub('[^a-z0-9]', '', tolower(basename(z)))

  # Process each image.
  unmatched = c()
  for(i in names(x@counts)){
    fp = NULL

    # 1. Explicit sample -> file mapping.
    if(!is.null(explicit_map) && i %in% names(explicit_map)){
      cand = explicit_map[[i]]
      if(any(file.exists(cand))){
        fp = cand[file.exists(cand)][1]
      }
    }

    # 2. Documented convention: the filename contains `<sample_name>.`.
    if(is.null(fp)){
      hit = grep(paste0(i, '\\.'), images, value=TRUE)
      if(length(hit) > 1){
        warning(paste0("Multiple image files matched sample ", i,
                       ". Using the first match (", hit[1], ")."), call.=FALSE)
      }
      if(length(hit) >= 1){
        fp = hit[1]
      }
    }

    # 3. Fallback: normalized prefix match (case- and separator-insensitive).
    if(is.null(fp)){
      key = norm_key(i)
      hit = images[startsWith(norm_key(images), key)]
      if(length(hit) >= 1){
        fp = hit[1]
      }
    }

    if(is.null(fp)){
      unmatched = c(unmatched, i)
      next
    }

    # Pick the decoder from the extension, but fall back to the other one: image
    # files are frequently mislabelled (e.g. a PNG saved as .jpg), and the
    # decoder's own error does not say which file it was reading.
    if(grepl('\\.jpe?g$', fp, ignore.case=TRUE)){ # CosMx uses jpeg files
      img_obj = tryCatch(jpeg::readJPEG(fp), error=function(e) NULL)
      if(is.null(img_obj)) img_obj = tryCatch(png::readPNG(fp), error=function(e) NULL)
    } else{
      img_obj = tryCatch(png::readPNG(fp), error=function(e) NULL)
      if(is.null(img_obj)) img_obj = tryCatch(jpeg::readJPEG(fp), error=function(e) NULL)
    }

    if(is.null(img_obj)){
      warning(paste0("Could not read image file '", fp, "' for sample ", i,
                     " as PNG or JPEG; skipping this sample."), call.=FALSE)
      unmatched = c(unmatched, i)
      next
    }

    # Downsize image if too large
    if(any(dim(img_obj) > 1000)){
      img_obj = EBImage::resize(img_obj, w=1000)
    }

    x@misc[['sp_images']][[i]] = img_obj
  }

  # Report unmatched samples once, with the convention to follow.
  if(length(unmatched) > 0){
    warning(paste0("No image found for sample(s): ", paste(unmatched, collapse=', '),
                   ". Image filenames must contain '<sample_name>.' (for example '",
                   unmatched[1], ".png'), or pass a named vector mapping sample names ",
                   "to file paths, e.g. images=c(", unmatched[1], "='path/to/image.png')."),
            call.=FALSE)
  }

  return(x)
}


##
# @title create_listw_from_knn
# @param x an STlist
# @param ks
#
create_listw_from_knn = function(x=NULL, ks=NULL){
  # Define cores available
  cores = 1

  # Create distance matrix based on the coordinates of each sampled location.
  listw_list = sge_lapply(seq(names(x@spatial_meta)), function(i){
    # See the note in create_listw_from_dist() below on remove_rownames().
    coords_mtx = as.matrix(x@spatial_meta[[i]][, c('libname', 'ypos', 'xpos')] %>%
                             tibble::remove_rownames() %>%
                             tibble::column_to_rownames('libname'))
    subj_listw = spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords_mtx, k=ks, longlat=F)), style='B')
    return(subj_listw)
  }, cores=cores, preschedule=F)
  names(listw_list) = names(x@spatial_meta)
  return(listw_list)
}


##
# @title create_listw_from_dist
# @param x an STlist
# @param ks
#
create_listw_from_dist = function(x=NULL, cores=NULL){
  # Define cores available
  cores = resolve_cores(length(x@spatial_meta), cores, coerce='as_integer')

  # Create distance matrix based on the coordinates of each sampled location.
  listw_list = sge_lapply(seq(names(x@spatial_meta)), function(i){
    # remove_rownames() before column_to_rownames(): @spatial_meta is a plain
    # data.frame, so any caller that row-subset it with base `[` leaves
    # non-contiguous row names behind, which column_to_rownames() rejects with
    # "`.data` must be a data frame without row names." Row names are not used
    # here -- adj's dimnames are overwritten from libname below -- so this is
    # output-neutral. Same footgun documented for STenrich()'s domain= path.
    coords_mtx = as.matrix(x@spatial_meta[[i]][, c('libname', 'ypos', 'xpos')] %>%
                             tibble::remove_rownames() %>%
                             tibble::column_to_rownames('libname'))
    adj = as.matrix(stats::dist(coords_mtx))
    adj = adj/max(adj)
    diag(adj) = NA
    colnames(adj) = x@spatial_meta[[i]][['libname']]
    rownames(adj) = x@spatial_meta[[i]][['libname']]

    neighbours = list()
    weights = list()
    for(id in 1:nrow(adj)){
      idx = 1:nrow(adj)
      neighbours[[id]] = idx[!(rownames(adj) %in% rownames(adj)[id])]
      names(neighbours[[id]]) = idx[!(rownames(adj) %in% rownames(adj)[id])]
      weights[[id]] = as.vector(adj[id, ])
      #weights[[id]] = weights[[id]][!(rownames(adj) %in% rownames(adj)[id])]
      weights[[id]] = 1/(weights[[id]][!(rownames(adj) %in% rownames(adj)[id])])
    }
    class(neighbours) = 'nb'

    subj_listw = list(style='B',
                      neighbours=neighbours,
                      weights=weights)
    class(subj_listw) = c("listw", "nb")

    return(subj_listw)
  }, cores=cores, preschedule=F)
  names(listw_list) = names(x@spatial_meta)
  return(listw_list)
}


##
#' @title get_gene_meta: Extract gene-level metadata and statistics
#' @description Extracts gene-level metadata and spatial statistics (if already computed)
#' @details
#' This function extracts data from the `x@gene_meta` slot, optionally subsetting
#' only to those genes for which spatial statistics (Moran's I or Geary's C, see `SThet`)
#' have been calculated. The output is a data frame with data from all samples in the
#' STlist
#'
#' @param x an STlist
#' @param sthet_only logical, return only genes with spatial statistics
#' @return a data frame with gene-level data
#'
#' @examples
#' \donttest{
#' # Using included melanoma example (Thrane et al.)
#' # Download example data set from spatialGE_Data
#' thrane_tmp = tempdir()
#' unlink(thrane_tmp, recursive=TRUE)
#' dir.create(thrane_tmp)
#' lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
#' tryCatch({ # In case data is not available from network
#'   download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
#'   zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
#'   unzip(zipfile=zip_tmp, exdir=thrane_tmp)
#'   # Generate the file paths to be passed to the STlist function
#'   count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='counts')
#'   coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                             full.names=TRUE, pattern='mapping')
#'   clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
#'                           full.names=TRUE, pattern='clinical')
#'   # Create STlist
#'   library('spatialGE')
#'   melanoma <- STlist(rnacounts=count_files,
#'                      spotcoords=coord_files,
#'                      samples=clin_file)
#'   melanoma <- transform_data(melanoma)
#'   melanoma <- SThet(melanoma, genes=c('MLANA', 'TP53'), method='moran')
#'   get_gene_meta(melanoma, sthet_only=TRUE)
#' }, error = function(e) {
#'   message("Could not run example. Are you connected to the internet?")
#'   return(NULL)
#' })
#' }
#'
#' @export
#
get_gene_meta = function(x=NULL, sthet_only=FALSE){
  # Test if an STList has been input.
  x = assert_stlist(x)
  # Test is gene_meta is empty (no normalization performed on STlist)
  if(rlang::is_empty(x@gene_meta)){
    stop('No gene-level metadata available in this STlist.')
  }

  # Loop through samples
  genemeta_dfs = list()
  for(i in names(x@gene_meta)){
    genemeta_dfs[[i]] = x@gene_meta[[i]] %>%
      tibble::add_column(sample=i, .before=1)
    # Subset genes if only SThet stats are required
    if(sthet_only & any(colnames(x@gene_meta[[i]]) %in% c('moran_i', 'geary_c'))){
      genemeta_dfs[[i]] = genemeta_dfs[[i]] %>%
        dplyr::filter(!is.na(moran_i) | !is.na(geary_c))
    }
  }
  res = do.call(dplyr::bind_rows, genemeta_dfs)

  # Check if sample meta data available and add to output
  if(nrow(x@sample_meta) > 0){
    df_tmp = x@sample_meta
    colnames(df_tmp)[1] = colnames(res)[1]
    res = dplyr::left_join(df_tmp, res, by="sample", multiple='all')

    rm(df_tmp) # Clean env
  }

  if(nrow(res) == 0){
    stop('No gene-level meta available. If only spatial statistics requested, please make sure SThet was run.')
  }
  return(res)
}


## Cache for raise_err()'s err.csv lookup table, populated lazily on first
## use rather than in .onLoad -- under devtools::load_all() system.file()
## resolves differently, and the csv is edited during development, so
## reading it at package-load time would risk caching a stale/wrong path.
.err_db_cache = new.env(parent=emptyenv())

##
# @title raise_err
#
raise_err = function(err_code=NULL, samplename=NULL, ...){
  if(is.null(.err_db_cache$err_db)){
    pkg_fp = list.files(system.file(package='spatialGE'), full.names=T, pattern='err\\.csv')
    .err_db_cache$err_db = suppressWarnings(utils::read.csv(pkg_fp, header=F, quote="'"))
  }
  err_db = .err_db_cache$err_db
  str_to_print = err_db[[2]][ err_db[[1]] == err_code ]

  # An unknown or missing code used to give character(0), and stop(character(0))
  # fails with "argument is of length zero" -- which names neither the typo'd code
  # nor the function. Name it here instead.
  if(length(str_to_print) != 1 || is.na(str_to_print)){
    stop(paste0('spatialGE exception: internal error - unknown error code ',
                if(is.null(err_code)) '<NULL>' else err_code, '.'), call.=F)
  }

  # Optional sprintf() formatting. Messages that embed gene names, sample sizes or
  # thresholds otherwise have to be written inline at the call site, where they
  # drift from this catalog and lose the shared wording.
  dots = list(...)
  if(length(dots) > 0){
    str_to_print = do.call(sprintf, c(list(str_to_print), dots))
  }

  if(is.null(samplename)){
    stop(paste0('spatialGE exception: ', str_to_print), call.=F)
  } else{
    stop(paste0('spatialGE exception: ', str_to_print, ' Sample: ', samplename, '.'), call.=F)
  }
}


##
# @title calculate_vst: Calculate standardized variation (VST) using Seurat's implementation
# @param x an STlist or a matrix with spots/cells in columns
# @param samples the samples from `x` to calculate VST. Ignored if x is a matrix
# @param cores integer indicating the number of cores in parallelization
# @return an STlist with VST in @gene_meta
#
calculate_vst = function(x=NULL, samples=NULL, cores=NULL){

  # To prevent NOTES in R CMD check
  . = NULL

  if(is(x, 'STlist')){
    # Identify variable genes (Seurat's VST)
    vst_ls = sge_lapply(samples, function(i){
      # Find tops variable genes using Seurat approach
      # First check if VST hasn't been calculated. If not, calculate
      if(!any(colnames(x@gene_meta[[i]]) == 'vst.variance.standardized')){
        vst_raw = Seurat_FindVariableFeatures(sample_counts(x, i))
        # Too few non-constant genes for the loess fit. Reported as NA rather
        # than warned here: this runs inside a fork, where a warning() is not
        # reliably delivered to the caller, so the parent loop below emits it.
        if(is.null(vst_raw)) return(NA)
        vst_tmp = vst_raw %>%
          tibble::rownames_to_column(var='gene') %>%
          dplyr::select('gene', 'vst.variance.standardized') %>%
          dplyr::left_join(x@gene_meta[[i]], ., by='gene')
        return(vst_tmp)
      } else{
        return(NULL)
      }
    }, cores=cores)
    names(vst_ls) = samples

    for(i in samples){
      if(identical(vst_ls[[i]], NA)){
        # Leave @gene_meta without the column rather than failing: the caller
        # asked for a transformation, and gene_mean/gene_stdevs are still there.
        # calculate_vst() only computes when the column is absent, so whichever
        # analysis next needs it will retry rather than inherit a wrong value.
        warning(paste0('Sample ', i, ' has fewer than 4 genes with non-zero variance, ',
                       'which is too few to fit the variance-mean trend. ',
                       'Skipping vst.variance.standardized for this sample; ',
                       'functions selecting genes by variability (topgenes) will not work on it.'))
      } else if(!is.null(vst_ls[[i]])){
        x@gene_meta[[i]] = vst_ls[[i]]
      }
    }
    return(x)
  } else{
    #vst_tmp = Seurat_FindVariableFeatures(sample_counts(x, i)) %>%
    vst_raw = Seurat_FindVariableFeatures(x)
    # The matrix form exists to rank genes by variability for a caller's
    # `topgenes`. If the trend cannot be fit there is no ranking to return, so
    # this is an error rather than a skip (error when nothing matches, warn when
    # only some does). The message names STgradient because that is the only
    # caller of this branch; a second caller should add its own code.
    if(is.null(vst_raw)){
      raise_err(err_code='error0097', samplename=NULL, nrow(x))
    }
    vst_tmp = vst_raw %>%
      tibble::rownames_to_column(var='gene') %>%
      dplyr::select('gene', 'vst.variance.standardized')
    return(vst_tmp)
  }
}

