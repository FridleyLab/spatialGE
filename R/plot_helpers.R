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
# @ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return, a ggplot object.
#
#
quilt_p <- function(data_f=NULL, color_pal="YlOrBr", leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T, ptsize=NULL){
  require('ggplot2')

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos, color=values)) +
    geom_point(size=ptsize) +
    scale_color_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(color=leg_name, title=title_name) +
    theme_classic()

    if(visium){
      p <- p + #scale_x_reverse() +
        scale_y_reverse() + coord_fixed(ratio=1.7)
    } else {
      p <- p + coord_fixed(ratio=1)
    }

    p <- p + theme(legend.position="right")

  return(p)
}

##
# @title quilt_p_purity
# @description Creates a quilt plot from ST data.
# @details
# Function to produce a "quilt plot" from a data frame with four columns:
# x coordinates, y coordinates, values (expression or cell scores), and tumor/stroma
# classes. It also takes a color palette name from the 'khroma' or RColorBrewer
# packages. Finally, it takes a name for the color legend title.
#
# @param data_f, a data with four columns: x coordinates, y coordinates, values to be
# plotted, and tumor/stroma classes.
# @param color_pal, a scheme from 'khroma'. Default is 'YlOrBr'.
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return, a ggplot object.
#
#
quilt_p_purity <- function(data_f=NULL, color_pal="YlOrBr", leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T, ptsize=0.5){
  require('ggplot2')

  colnames(data_f) = c('Y', 'X', 'Value', 'Tumor_Stroma')
  data_f$Tumor_Stroma = as.factor(data_f$Tumor_Stroma)

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Create plot.
  p1 <- ggplot() +
    geom_point(data=data_f, aes(x=X, y=Y, fill=Value, shape=Tumor_Stroma, color=Value), stroke=1, size=ptsize) +
    scale_fill_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    scale_color_gradientn(colours=p_palette, limits=c(minvalue, maxvalue)) +
    scale_shape_manual(values=c(3, 15)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name, shape='tumor/stroma',
         title=title_name, color=leg_name) +
    guides(shape=guide_legend(override.aes=list(color='black', stroke=0.5, size=2), nrow=2, byrow=T, title=NULL)) +
    theme_classic()

  if(visium){
    p1 <- p1 + #scale_x_reverse() +
      scale_y_reverse() + coord_fixed(ratio=1.7)
  } else{
    p1 <- p1 + coord_fixed(ratio=1)
  }

  p1 <- p1 + theme(legend.position="right")

  return(p1)
}

##
# @title quilt_p_purity_bw
# @description Creates a tumor/stroma quilt plot from ST data without showing expression values.
# @details
# Function to produce a "quilt plot" from a data frame with three columns:
# x coordinates, y coordinates, and values (expression or cell scores). The data
# frame has column names 'x_pos', 'y_pos', and 'values'.
#
# @param data_f, a data with four columns: x coordinates, y coordinates,values to be
# plotted, and tumor/stroma classes. The third column is not used.
# @param title_name, a short name for the plot title.
# @param visium, whether or not to reverse axes for Visium slides.
# @ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return, a ggplot object.
#
#
quilt_p_purity_bw <- function(data_f=NULL, title_name='', visium=T, ptsize=NULL){
  require('ggplot2')

  colnames(data_f) = c('Y', 'X', 'Value', 'Tumor_Stroma')
  data_f$Tumor_Stroma = as.factor(data_f$Tumor_Stroma)

  p2 <- ggplot() +
    geom_point(data=data_f, aes(x=X, y=Y, color=Tumor_Stroma, shape=Tumor_Stroma), size=ptsize) +
    scale_shape_manual(values=c(3, 15)) +
    scale_color_manual(values=c('gray60', 'black')) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(shape='tumor/stroma', color='tumor/stroma', title=) +
    guides(color=guide_legend(override.aes=list(size=2))) +
    theme_classic()


    if(visium){
      p2 <- p2 + #scale_x_reverse() +
        scale_y_reverse() + coord_fixed(ratio=1.7)
    } else{
      p2 <- p2 + coord_fixed(ratio=1)
    }

  return(p2)
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
#
krige_p <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', minvalue=minvalue, maxvalue=maxvalue,
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
    labs(fill=leg_name, title=title_name) +
    ggpolypath::geom_polypath(aes(long,lat,group=group), mask_df, fill="white"
    ) +
    theme_classic()

  if(visium){
    p <- p + scale_y_reverse() + #scale_x_reverse() +
      coord_fixed(ratio=1.7)
  } else{
    p <- p + coord_fixed(ratio=1)
  }

  p <- p + theme(legend.position="right")

  return(p)

}
