##
#' @title plot_image
#' @description
#' @details
#'
#'
#' @param x an STlist
#' @param samples a vector of numbers indicating the ST samples to plot, or their
#' sample names. If vector of numbers, it follow the order of `names(x@counts)`.
#' If NULL, the function plots all samples
#' @return a list of plots
#'
# @importFrom methods as is new
#' @import ggplot2
# @importFrom magrittr %>%
#
#' @export
#
#
plot_image = function(x=NULL, samples=NULL){

  # Define which samples to plot
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = names(x@spatial_meta)[samples]
    }
    if(length(grep(paste0(samples, collapse='|'), names(x@spatial_meta))) == 0){
      stop('The requested samples are not present in the STlist.')
    }
  }

  # Check which samples have images stored
  rm_samples = c()
  for(i in samples){
    if(is.null(x@misc[['sp_images']][[i]])){
      rm_samples = append(rm_samples, i)
    }
  }
  samples = samples[ !(samples %in% rm_samples) ]

  # Create list of plots.
  rp_list <- list()
  # Loop through each normalized count matrix.
  for (i in samples){
    # Read image if available
    if(!is.null(x@misc[['sp_images']][[i]])){
      img_obj = grid::rasterGrob(x@misc[['sp_images']][[i]])
    }

    # Get scaling factor is available
    if(is.list(x@misc[['image_scaling']][[i]]) & x@misc[['platform']] == 'visium'){
      scl_info = x@misc[['image_scaling']][[i]][['tissue_hires_scalef']]
      # Scale max and min coordinates to match pixels
      raster_extext = list(ymin=round(min(x@spatial_meta[[i]][['ypos']]) * scl_info, 0),
                           ymax=round(max(x@spatial_meta[[i]][['ypos']]) * scl_info, 0),
                           xmin=round(min(x@spatial_meta[[i]][['xpos']]) * scl_info, 0),
                           xmax=round(max(x@spatial_meta[[i]][['xpos']]) * scl_info, 0))
      img_obj[['raster']] = img_obj[['raster']][raster_extext$ymin:raster_extext$ymax, raster_extext$xmin:raster_extext$xmax]
    }

    rp_list[[paste0('image_', i)]] = ggplot() +
      annotation_custom(img_obj)
  }
  return(rp_list)
}

