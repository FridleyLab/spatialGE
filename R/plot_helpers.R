##
# @title quilt_p
# @description Creates a quilt plot from ST data.
# @details
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values (expression or cell scores). The data
# frame has column names 'x_pos', 'y_pos', and 'values'. It also takes a color
# palette name from the 'khroma' or RColorBrewer packages. Finally, it takes a
# name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the values to be plotted.
# @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return, a ggplot object.
#
#' @import ggplot2
#
#
quilt_p <- function(data_f=NULL, color_pal="BuRd", leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T, ptsize=0.5){

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Create plot.
  p <- ggplot2::ggplot(data=data_f, ggplot2::aes(x=xpos, y=ypos, color=values)) +
    ggplot2::geom_point(size=ptsize) +
    ggplot2::scale_color_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    ggplot2::xlab("X Position") +
    ggplot2::ylab("Y Position") +
    ggplot2::labs(color=leg_name, title=title_name) +
    ggplot2::theme_void()
    #ggplot2::theme_classic()

  if(visium){
    p <- p + ggplot2::scale_y_reverse()
    # scale_y_reverse() + coord_fixed(ratio=1.7)
  }

  p <- p + ggplot2::coord_equal() + theme(legend.position="right")

  return(p)
}


##
# @title krige_p
# @description Creates a kriging plot from ST data.
# @details
# Function to produce a "kriging plot" from a data frame with three columns:
# x coordinates, y coordinates, and predicted kriging values. The data frame has
# column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
# SpatialPolygons object to mask the predicted grid to the area of the tissue.
# It also takes a color palette name from the 'khroma' package. Finally, it
# takes a name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the kriging prediction values to be plotted.
# @param mask, an object of class SpatialPolygons containing a large polygon
# encasing all the predicted grid, and a smaller polygon drawing the concave hull
# of the tissue shape.
# @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @return, a ggplot object.
#
#' @import ggplot2
#' @import ggpolypath
#
#
krige_p <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', minvalue=minvalue, maxvalue=maxvalue,
                    visium=T){

  #requireNamespace('sf')

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Convert the SpatialPolygon mask into a data frame.
  #mask_df <- ggplot2::fortify(mask)

  # Create plot.
  p <- ggplot2::ggplot(data=data_f, ggplot2::aes(x=x_pos, y=y_pos)) +
    ggplot2::geom_raster(ggplot2::aes(fill=krige), interpolate=F) +
    ggplot2::scale_fill_gradientn(colors=p_palette, limits=c(minvalue, maxvalue), oob=scales::squish) +
    ggplot2::xlab("X Position") +
    ggplot2::ylab("Y Position") +
    ggplot2::labs(fill=leg_name, title=title_name) +
    ggplot2::geom_sf(data=mask, color='white', fill="white", linewidth=2, inherit.aes=F) +
    ggplot2::theme_void()

  # if(visium){
  #   p <- p + ggplot2::scale_y_reverse() #+ scale_x_reverse() +
  #     #coord_fixed(ratio=1.7)
  # }

  p <- p + #ggplot2::coord_equal() +
    ggplot2::theme(legend.position="right")

  return(p)
}

