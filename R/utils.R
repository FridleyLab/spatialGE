##
#' @title xcell_names: Prints xCell cell type names
#' @description Prints the available cell types as provided by xCell deconvolution.
#' @details
#' Shows the names of cell types that can be used in other spatialGE functions.
#' This function access the xCell cell type names after cleaning via the janitor
#' package.
#'
#' @return The xCell cell type names after cleaning.
#'
#' @examples
#' # xcell_names()
#'
#' @export
#'
#'
xcell_names = function(){
  #xcellnames = x@cell_deconv$xCell[[1]]$cell_stdev$cell
  xcellnames = rownames(xCell::xCell.data$spill$K)
  xcellnames = janitor::make_clean_names(xcellnames)
  print(xcellnames)
}

##
# Get the number of cores to use in parallel as a function of the number of
# data to process. It will not yield a higher number of cores than half of the
# total available cores. Will default to 1 core if not Unix
# @param n, an integer representing the number of units to process
# @return cores, the number of cores to use
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
      cores = n
    }
  }
  return(cores)
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
# @paran n_cats The number of colors to produce. If NULL, assumes 5 colors.
# @return cat_cols A color palette
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
    p_palette = khroma::colour(color_pal[1])
    cat_cols = as.vector(p_palette(n_cats))
  }else if(color_pal[1] %in% rownames(RColorBrewer::brewer.pal.info)){
    cat_cols = RColorBrewer::brewer.pal(n_cats, color_pal[1])
  }else{ # Test if user provided a vector of colors.
    if(length(color_pal) >= n_cats){
      cat_cols = color_pal[1:n_cats]
    } else{
      stop('Provide enough colors to plot or an appropriate Khroma/RColorBrewer palette.')
    }
  }

  return(cat_cols)
}

##
#' @title load_images: Place tissue images within STList
#' @description Loads the images from tissues to the appropriate STList slot.
#' @details
#' This function looks for `.PNG` files within a folder matching the sample names
#' in an existing STList. Then, loads the images to the STList which can be used
#' for plotting along with quilt plots.
#'
#' @param x an STList with transformed RNA counts.
#' @param images a string indicating a folder to load images from
#' @return an STList with images
#'
#' @export
#
#
load_images = function(x=NULL, images=NULL){
  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  if(is.null(images)){
    stop("Please, provide a directory with images.")
  }

  # Get file names of images
  imageFps = list.files(images, full.names=T)

  # Process each image.
  for(i in names(x@counts)){
    fp = grep(i, imageFps, value=T)
    if(length(fp) == 0){
      cat(paste0("Image for sample ", i, " was not found."))
      next
    }

    # Return warning if more than one file matches the sample name
    if(length(fp) > 1){
      cat(paste0("Multiple image files matched sample ", i, ". Using the first match (", fp[1], ")."))
    }

    img_obj = png::readPNG(fp[1])

    # Downsize image if too large
    if(any(dim(img_obj) > 1000)){
      img_obj = EBImage::resize(img_obj, w=1000)
    }

    x@misc[['sp_images']][[i]] = img_obj
  }
  return(x)
}

