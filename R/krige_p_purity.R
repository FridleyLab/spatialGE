##
# @title krige_p
# @description Creates a kriging plot from ST data.
# @details
# Function to produce a "kriging plot" from a data frame with three columns:
# x coordinates, y coordinates, and predicted kriging values. The data frame has
# column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
# SpatialPolygons object to mask the predicted grid to the area of the tissue.
# It also takes a color palette name from the 'khroma' package and a name for the
# color legend title. The function requires also a dataframe indicating which
# spots are tumor or stroma.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the kriging prediction values to be plotted.
# @param tumorstroma, a data frame with x,y positions and stroma/tumor assignments.
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
#
krige_p_purity <- function(data_f=NULL, tumorstroma=NULL, mask=NULL, color_pal="YlOrBr",
                           leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue,
                           visium=T){
  require('ggplot2')

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Convert the SpatialPolygon mask into a data frame.
  mask_df <- fortify(mask)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos)) +
    geom_raster(aes(fill=krige)) +
    scale_fill_gradientn(colors=p_palette, limits=c(minvalue, maxvalue), oob=scales::squish) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name, shape='', title=title_name) +
    ggpolypath::geom_polypath(aes(long,lat,group=group), mask_df, fill="white"#,
                              #color='white', size=0
                              ) +
    geom_point(data=tumorstroma[tumorstroma$cluster=='tumor', ], aes(x=xpos, y=ypos, shape=cluster), size=0.7, color='gray50') +
    scale_shape_manual(values=c(0, 3)) +
    guides(shape=guide_legend(override.aes=list(size=3))) +
    theme_classic() +
    theme(legend.position="right", plot.title=element_text(size=8), legend.text = element_text(size=10))

    if(visium){
      p <- p + scale_x_reverse() + scale_y_reverse() + coord_fixed(ratio=1.7)
    } else{
      p <- p + coord_fixed(ratio=1)
    }

  return(p)

}
