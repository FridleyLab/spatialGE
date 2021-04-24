##
# Function to produce a "kriging plot" from a data frame with three columns:
# x coordinates, y coordinates, and prdicted kriging values. The data frame has
# column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
# SpatialPolygons object to mask the predicted grid to the area of the tissue.
# It also takes a color palette name from the 'khroma' package. Finally, it
# takes a name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the kriging prediction values to be plotted.
# @param mask, an object of class SpatialPolygons containing a large polygon
# encasing all the predicted grid, and a smalle polygon drawing the concave hull
# of the tissue shape.
# @param color_pal, a scheme from 'khroma'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @return, a ggplot object.
#
#
require('ggplot2')
krige_p <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                    title_name=''){

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  # Convert the SpatialPolygon mask into a data frame.
  mask_df <- fortify(mask)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos)) +
    geom_raster(aes(fill=krige)) +
    scale_fill_gradientn(colors=p_palette(5)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name) +
    ggtitle(title_name) +
    ggpolypath::geom_polypath(aes(long,lat,group=group), mask_df, fill="white",
                              color='white', size=1.5) +
    theme_classic() +
    theme(legend.position="right", plot.title=element_text(size=10))

  return(p)

}
