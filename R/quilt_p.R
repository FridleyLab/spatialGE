##
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values. The data frame has column names 'x_pos',
# 'y_pos', and 'values'. It also takes a color palette name from the 'khroma'
# package. Finally, it takes a name for the color legend title.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# values to be color coded.
# @color_pal, a scheme from 'khroma'.
# @leg_name, a short name for the legend title.
# @title_name, a short name for the plot title.
# @return, a ggplot object.
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
