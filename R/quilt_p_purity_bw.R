##
# @title quilt_p_purity_bw
# @description Creates a tumor/stroma quilt plot from ST data without showing expression values.
# @details
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values (expression or cell scores). The data
# frame has column names 'x_pos', 'y_pos', and 'values'.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the values to be plotted.
# @return, a ggplot object.
#
#
quilt_p_purity_bw <- function(data_f=NULL){

  require('ggplot2')

  p2 <- ggplot() +
    geom_point(data=data_f, aes(x=x_pos, y=y_pos, color=as.factor(cluster), shape=as.factor(cluster)), size=0.7) +
    scale_shape_manual(values=c(3, 15)) +
    scale_color_manual(values=c('gray60', 'black')) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(shape='tumor/stroma', color='tumor/stroma', title='ESTIMATE\ntumor/stroma') +
    guides(color=guide_legend(override.aes=list(size=2))) +
    theme_classic() +
    scale_x_reverse() +
    scale_y_reverse() +
    coord_fixed()

  #purity_p_list <- list()
  #purity_p_list[['p1']] <- p1
  #purity_p_list[['p2']] <- p2

  #purity_gb <- gridExtra::arrangeGrob(grobs=purity_p_list, heights = c(1, 1), layout_matrix = rbind(c(1), c(2)))

  #return(purity_gb)
  return(p2)
}
