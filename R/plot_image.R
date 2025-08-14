##
#' @title plot_image: Generate a ggplot object of the tissue image
#' @description Creates ggplot objects of the tissue images when available within the STlist
#' @details
#' If the STlist contains tissue images in the `@misc` slot, the `plot_image` function
#' can be used to generate ggplot objects. These ggplot objects can be plotted next to
#' quilt plots (`STplot` function) for comparative analysis.
#'
#' @param x an STlist
#' @param samples a vector of numbers indicating the ST samples to plot, or their
#' sample names. If vector of numbers, it follow the order of `names(x@counts)`.
#' If NULL, the function plots all samples
#' @return a list of plots
#'
#' @import ggplot2
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

  # Find which samples have images
  samples_tmp = samples
  for(i in samples){
    if(!(i %in% names(x@misc[['sp_images']]))){
      samples_tmp = grep(i, samples_tmp, value=T, invert=T)
    }
  }
  samples = samples_tmp
  rm(samples_tmp)

  # Create list of plots.
  rp_list <- list()

  if(length(samples) < 1){
    warning('No tissue images available for the samples in this STlist.')
  } else{
    # Loop through each normalized count matrix.
    for (i in samples){
      # Read image if available
      if(!is.null(x@misc[['sp_images']][[i]])){
        img_obj = grid::rasterGrob(x@misc[['sp_images']][[i]])
      }

      # Get scaling factor is available (VISIUM)
      if(is.list(x@misc[['image_scaling']][[i]]) & x@misc[['platform']] == 'visium'){
        scl_info_hires = x@misc[['image_scaling']][[i]][['tissue_hires_scalef']]
        scl_info_lores = x@misc[['image_scaling']][[i]][['tissue_lowres_scalef']]

        # Attempt scaling with hires scaling factor
        image_test_hires = tryCatch({
          raster_extext = list(ymin=round(min(x@spatial_meta[[i]][['ypos']]) * scl_info_hires, 0),
                               ymax=round(max(x@spatial_meta[[i]][['ypos']]) * scl_info_hires, 0),
                               xmin=round(min(x@spatial_meta[[i]][['xpos']]) * scl_info_hires, 0),
                               xmax=round(max(x@spatial_meta[[i]][['xpos']]) * scl_info_hires, 0))
          image_tmp = img_obj
          image_tmp[['raster']] = image_tmp[['raster']][raster_extext$ymin:raster_extext$ymax, raster_extext$xmin:raster_extext$xmax]
        },
        error=function(e){
          return(e)
        })

        # Attempt scaling with lowres scaling factor
        if(any(class(image_test_hires) == 'simpleError')){
          image_test_lores = tryCatch({
            raster_extext = list(ymin=round(min(x@spatial_meta[[i]][['ypos']]) * scl_info_hires, 0),
                                 ymax=round(max(x@spatial_meta[[i]][['ypos']]) * scl_info_hires, 0),
                                 xmin=round(min(x@spatial_meta[[i]][['xpos']]) * scl_info_hires, 0),
                                 xmax=round(max(x@spatial_meta[[i]][['xpos']]) * scl_info_hires, 0))
            image_tmp = img_obj
            image_tmp = image_tmp[raster_extext$ymin:raster_extext$ymax, raster_extext$xmin:raster_extext$xmax]
          },
          error=function(e){
            return(e)
          })
        }

        # Re-assign object whether modified or not
        img_obj = image_tmp
      }

      rp_list[[paste0('image_', i)]] = ggplot() +
        annotation_custom(img_obj)
    }
  }
  return(rp_list)
}

