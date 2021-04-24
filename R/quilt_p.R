##
#' Function to produce a "kriging plot" from a data frame with three columns:
#' x coordinates, y coordinates, and predicted kriging values. The data frame has
#' column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
#' SpatialPolygons object to mask the predicted grid to the area of the tissue.
#' It also takes a color palette name from the 'khroma' package. Finally, it
#' takes a name for the color legend title.
#'
#' @param data_f, a data with three columns: x coordinates, y coordinates, and
#' the kriging prediction values to be plotted.
#' @param mask, an object of class SpatialPolygons containing a large polygon
#' encasing all the predicted grid, and a smaller polygon drawing the concave hull
#' of the tissue shape.
#' @param color_pal, a scheme from 'khroma'.
#' @param leg_name, a short name for the legend title.
#' @param title_name, a short name for the plot title.
#' @return, a ggplot object.
#
#
require('ggplot2')
quilt_p <- function(data_f=NULL, color_pal="YlOrBr", leg_name='',
                    title_name=''){

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos, color=values)) +
    geom_point() +
    scale_color_gradientn(colours=p_palette(5)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(color=leg_name) +
    theme_classic() +
    ggtitle(title_name) +
    theme(legend.position="right")

  return(p)

}
