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
#' @param images a string indicating a folder to load images from
#' @return an STlist with images
#'
#' @export
#'
#' @importFrom methods as is new
#
#
load_images = function(x=NULL, images=NULL){
  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STlist.")
  }

  if(is.null(images)){
    stop("Please, provide a vector with images file paths.")
  } else{
    if(length(images) == 1 & dir.exists(images)){ # In case a directory was provided
      images = list.files(images, full.names=T)
    }
  }

  # Process each image.
  for(i in names(x@counts)){
    fp = grep(paste0(i, '\\.'), images, value=T)
    if(length(fp) == 0){
      message(paste0("Image for sample ", i, " was not found.\n"))
      next
    }

    # Return warning if more than one file matches the sample name
    if(length(fp) > 1){
      message(paste0("Multiple image files matched sample ", i, ". Using the first match (", fp[1], ")."))
    }

    if(grepl('jpg$', fp[1])){ # CosMx uses jpeg files
      img_obj = jpeg::readJPEG(fp[1])
    } else{
      img_obj = png::readPNG(fp[1])
    }

    # Downsize image if too large
    if(any(dim(img_obj) > 1000)){
      img_obj = EBImage::resize(img_obj, w=1000)
    }

    x@misc[['sp_images']][[i]] = img_obj
  }
  return(x)
}


##
# @title create_listw
# @param x an STlist
#
create_listw = function(x=NULL){
  # Define cores available
  cores = count_cores(length(x@spatial_meta))

  # Create distance matrix based on the coordinates of each sampled location.
  listw_list = parallel::mclapply(seq(names(x@spatial_meta)), function(i){
    subj_dists = as.matrix(dist(x@spatial_meta[[i]][, c('ypos', 'xpos')]))
    subj_dists[subj_dists == 0] = 0.0001
    subj_dists_inv = 1/subj_dists
    diag(subj_dists_inv) = 0
    subj_dists_inv=spdep::mat2listw(subj_dists_inv, style='B')
    return(subj_dists_inv)
  }, mc.cores=cores, mc.preschedule=F)
  names(listw_list) = names(x@spatial_meta)
  return(listw_list)
}


##
# @title create_listw_from_knn
# @param x an STlist
# @param ks
#
create_listw_from_knn = function(x=NULL, ks=NULL){
  # Define cores available
  #cores = count_cores(length(x@spatial_meta))
  cores = 1

  # Create distance matrix based on the coordinates of each sampled location.
  listw_list = parallel::mclapply(seq(names(x@spatial_meta)), function(i){
    coords_mtx = as.matrix(x@spatial_meta[[i]][, c('libname', 'ypos', 'xpos')] %>% tibble::column_to_rownames('libname'))
    subj_listw = spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(coords_mtx, k=ks, longlat=F)), style='B')
    return(subj_listw)
  }, mc.cores=cores, mc.preschedule=F)
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
  if(.Platform$OS.type == 'windows'){
    cores = 1
  }
  if(is.null(cores)){
    cores = count_cores(length(x@spatial_meta))
  } else{
    cores = as.integer(cores)
    if(is.na(cores)){
      stop('Could not recognize number of cores requested')
    }
  }

  # Create distance matrix based on the coordinates of each sampled location.
  listw_list = parallel::mclapply(seq(names(x@spatial_meta)), function(i){
    coords_mtx = as.matrix(x@spatial_meta[[i]][, c('libname', 'ypos', 'xpos')] %>% tibble::column_to_rownames('libname'))
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
  }, mc.cores=cores, mc.preschedule=F)
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
#'   #' zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
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
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STlist.")
  }
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


##
# @title raise_err
#
raise_err = function(err_code=NULL, samplename=NULL){
  pkg_fp = list.files(system.file(package='spatialGE'), full.names=T, pattern='err\\.csv')
  err_db = suppressWarnings(utils::read.csv(pkg_fp, header=F, quote="'"))
  str_to_print = err_db[[2]][ err_db[[1]] == err_code ]
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
    vst_ls = parallel::mclapply(samples, function(i){
      # Find tops variable genes using Seurat approach
      # First check if VST hasn't been calculated. If not, calculate
      if(!any(colnames(x@gene_meta[[i]]) == 'vst.variance.standardized')){
        vst_tmp = Seurat_FindVariableFeatures(x@counts[[i]]) %>%
          tibble::rownames_to_column(var='gene') %>%
          dplyr::select('gene', 'vst.variance.standardized') %>%
          dplyr::left_join(x@gene_meta[[i]], ., by='gene')
        return(vst_tmp)
      } else{
        return(NULL)
      }
    }, mc.cores=cores)
    names(vst_ls) = samples

    for(i in samples){
      if(!is.null(vst_ls[[i]])){
        x@gene_meta[[i]] = vst_ls[[i]]
      }
    }
    return(x)
  } else{
    #vst_tmp = Seurat_FindVariableFeatures(x@counts[[i]]) %>%
    vst_tmp = Seurat_FindVariableFeatures(x) %>%
      tibble::rownames_to_column(var='gene') %>%
      dplyr::select('gene', 'vst.variance.standardized')
    return(vst_tmp)
  }
}

