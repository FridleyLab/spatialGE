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
quilt_p <- function(data_f=NULL, color_pal="YlOrBr", leg_name='', title_name='', minvalue=minvalue, maxvalue=maxvalue, visium=T, ptsize=0.5){
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
quilt_p_purity_bw <- function(data_f=NULL, title_name='', visium=T, ptsize=0.5){
  require('ggplot2')

  colnames(data_f) = c('Y', 'X', 'Value', 'Tumor_Stroma')
  data_f$Tumor_Stroma = as.factor(data_f$Tumor_Stroma)

  p2 <- ggplot() +
    geom_point(data=data_f, aes(x=X, y=Y, color=Tumor_Stroma, shape=Tumor_Stroma), size=ptsize) +
    scale_shape_manual(values=c(3, 15)) +
    scale_color_manual(values=c('gray60', 'black')) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(shape='tumor/stroma', color='tumor/stroma', title=title_name) +
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

##
# @title krige_p_pvals
# @description Creates a kriging plot from ST data with xCell pvalues.
# @details
# Function to produce a "kriging plot" from a data frame with three columns:
# x coordinates, y coordinates, and predicted kriging values. The data frame has
# column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
# SpatialPolygons object to mask the predicted grid to the area of the tissue.
# It also takes a color palette name from the 'khroma' package. Finally, it
# takes a name for the color legend title. Finally, the function overlays dots
# for spots in which significant p-values (<0.05) were observed for a cell type.
#
# @param data_f, a data with three columns: x coordinates, y coordinates, and
# the kriging prediction values to be plotted.
# @param mask, an object of class SpatialPolygons containing a large polygon
# encasing all the predicted grid, and a smaller polygon drawing the concave hull
# of the tissue shape.
# @param color_pal, a scheme from 'khroma'. Defaul is 'YlOrBr'
# @param leg_name, a short name for the legend title.
# @param title_name, a short name for the plot title.
# @param x, an STList from where coordinates will be taken.
# @param plot_who, an integer indicating the spatial array to be plotted.
# @param cell, a cell name from a deconvoluted matrix.
# @param minvalue, the minimum value of gene expression or cell score. Used for
# standardization.
# @param maxvalue, the maximum value of gene expression or cell score. Used for
# standardization.
# @param visium, whether or not to reverse axes for Visium slides.
# @return, a ggplot object.
#
#
krige_p_pvals <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                          title_name='', x=NULL, plot_who=NULL, cell=NULL,
                          minvalue=minvalue, maxvalue=maxvalue, visium=T){

  require('ggplot2')

  # Creates color palette function.
  p_palette = color_parse(color_pal)

  # Convert the SpatialPolygon mask into a data frame.
  mask_df <- fortify(mask)

  # Create data freame of points to be plotted where xCells scores were significant (p<0.05).
  cell_pvals <- unlist(
    x@cell_deconv$xCell[[plot_who]]$pvals[x@cell_deconv$xCell[[plot_who]]$pvals$cell_names == cell, ])
  cell_coords <- x@coords[[plot_who]]
  pvals_df <- dplyr::bind_cols(cell_coords, pval=cell_pvals[-1])
  pvals_sign <- pvals_df[pvals_df$pval < 0.05, ]

  # Create plot.
  p <- ggplot(data=data_f, aes(x=x_pos, y=y_pos)) +
    geom_raster(aes(fill=krige)) +
    scale_fill_gradientn(colors=p_palette, limits=c(minvalue, maxvalue), oob=scales::squish) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name) +
    ggtitle(title_name) +
    ggpolypath::geom_polypath(aes(long,lat,group=group), mask_df, fill="white"#,
                              #color='white', size=1.5
    ) +
    geom_point(data=pvals_sign, aes(x=xpos, y=ypos), shape=1, size=0.7, color='gray50') +
    theme_classic() +
    theme(legend.position="right", plot.title=element_text(size=8))

  if(visium){
    p <- p + scale_x_reverse() + scale_y_reverse() + coord_fixed(ratio=1.7)
  } else{
    p <- p + coord_fixed(ratio=1)
  }

  return(p)

}

##
# @title krige_p_purity
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
