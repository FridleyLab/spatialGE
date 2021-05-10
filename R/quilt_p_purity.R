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

  tumor_df <- data_f[data_f$cluster == 1, ]
  stroma_df <- data_f[data_f$cluster == 2, ]

  # Create plot.
  p <- ggplot() +
    geom_point(data=tumor_df, aes(x=x_pos, y=y_pos, fill=values, stroke=0), shape=24) +
    geom_point(data=stroma_df, aes(x=x_pos, y=y_pos, fill=values, stroke=0), shape=22) +
    scale_fill_gradientn(colours=p_palette(5)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(color=leg_name) +
    theme_classic() +
    ggtitle(title_name) +
    coord_fixed() +
    theme(legend.position="right")

  return(p)

}
