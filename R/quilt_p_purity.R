##
# @title quilt_p_purity
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
quilt_p_purity <- function(data_f=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T){
  require('ggplot2')

  colnames(data_f) = c('Y', 'X', 'Value', 'Tumor_Stroma')
  data_f$Tumor_Stroma = as.factor(data_f$Tumor_Stroma)

  # Creates color palette function.
  #p_palette <- khroma::colour(color_pal)
  p_palette = color_parse(color_pal)

  # Create plot.
  p1 <- ggplot() +
    geom_point(data=data_f, aes(x=X, y=Y, fill=Value, shape=Tumor_Stroma, color=Value), stroke=1, size=0.5) +
    scale_fill_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    scale_color_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    scale_shape_manual(values=c(3, 15)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name, shape='tumor/stroma',
         title=title_name, color=leg_name) +
    guides(shape=guide_legend(override.aes=list(color='black', stroke=0.5, size=2))) +
    theme_classic()

  if(visium){
    p1 <- p1 + scale_x_reverse() + scale_y_reverse() + coord_fixed(ratio=1.7)
  } else{
    p1 <- p1 + coord_fixed(ratio=1)
  }

  p1 <- p1 + theme(legend.position="right")

  return(p1)
}
