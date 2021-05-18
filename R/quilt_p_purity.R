##
#' @title quilt_p_purity
#' @description Creates a quilt plot from ST data.
#' @details
#' Function to produce a "quilt plot" from a data frame with three columns:
#' x coordinates, y coordinates, and values (expression or cell scores). The data
#' frame has column names 'x_pos', 'y_pos', and 'values'. It also takes a color
#' palette name from the 'khroma' package. Finally, it takes a name for the color
#' legend title.
#'
#' @param data_f, a data with three columns: x coordinates, y coordinates, and
#' the values to be plotted.
#' @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
#' @param leg_name, a short name for the legend title.
#' @param title_name, a short name for the plot title.
#' @return, a ggplot object.
#
#
quilt_p_purity <- function(data_f=NULL, color_pal="YlOrBr", leg_name='',
                    title_name=''){

  require('ggplot2')

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  data_f$values <- data_f$values/max(data_f$values)

  # Create plot.
  p <- ggplot() +
    #geom_point(data=data_f, aes(x=x_pos, y=y_pos, fill=values, shape=as.factor(cluster), color=as.factor(cluster)), stroke=0.2) +
    geom_point(data=data_f, aes(x=x_pos, y=y_pos, fill=values, shape=as.factor(cluster), color=values), stroke=1, size=0.5) +
    scale_fill_gradientn(colours=p_palette(5)) +
    #scale_color_manual(values=c('white', 'gray50')) +
    scale_color_gradientn(colours=p_palette(5)) +
    #scale_shape_manual(values=c(21, 22)) +
    scale_shape_manual(values=c(3, 15)) +
    #guides(color=guide_legend(override.aes=list(color='black', stroke=0.5, size=2))) +
    xlab("X Position") +
    ylab("Y Position") +
    #labs(fill=leg_name, shape='tumor/stroma', color='tumor/stroma', title=title_name) +
    labs(fill=leg_name, shape='tumor/stroma', title=title_name, color=leg_name) +
    theme_classic() +
    scale_x_reverse() +
    scale_y_reverse() +
    coord_fixed() +
    theme(legend.position="right")

  return(p)

}
