##
# @title quilt_p
# @description Creates a quilt plot from ST data.
# @details
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values (expression or cell scores). The data
# frame has column names 'x_pos', 'y_pos', and 'values'. It also takes a color
# palette name from the 'khroma' package. Finally, it takes a name for the color
# legend title.
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
# @return, a ggplot object.
#
#
quilt_p <- function(data_f=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T){

  require('ggplot2')

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  #data_f$values <- data_f$values/max(data_f$values)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos, color=values)) +
    geom_point() +
    scale_color_gradientn(colours=p_palette(5), limits=c(minvalue, maxvalue)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(color=leg_name, title=title_name) +
    theme_classic()

    if(visium){
      p <- p + scale_x_reverse() + scale_y_reverse() + coord_fixed(ratio=1.7)
    } else {
      p <- p + coord_fixed(ratio=1)
    }

    p <- p + theme(legend.position="right")

  return(p)

}
