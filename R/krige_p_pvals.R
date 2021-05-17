##
#' @title krige_p_pvals
#' @description Creates a kriging plot from ST data with xCell pvalues.
#' @details
#' Function to produce a "kriging plot" from a data frame with three columns:
#' x coordinates, y coordinates, and predicted kriging values. The data frame has
#' column names 'x_pos', 'y_pos', and 'krige'. The function also takes a
#' SpatialPolygons object to mask the predicted grid to the area of the tissue.
#' It also takes a color palette name from the 'khroma' package. Finally, it
#' takes a name for the color legend title. Finally, the function overlays dots
#' for spots in which significant p-values (<0.05) were observed for a cell type.
#'
#' @param data_f, a data with three columns: x coordinates, y coordinates, and
#' the kriging prediction values to be plotted.
#' @param mask, an object of class SpatialPolygons containing a large polygon
#' encasing all the predicted grid, and a smaller polygon drawing the concave hull
#' of the tissue shape.
#' @param color_pal, a scheme from 'khroma'. Defaul is 'YlOrBr'
#' @param leg_name, a short name for the legend title.
#' @param title_name, a short name for the plot title.
#' @param x, an STList from where coordinates will be taken.
#' @param plot_who, an integer indicating the spatial array to be plotted.
#' @param cell, a cell name from a deconvoluted matrix.
#' @return, a ggplot object.
#
#
krige_p_pvals <- function(data_f=NULL, mask=NULL, color_pal="YlOrBr", leg_name='',
                    title_name='', x=NULL, plot_who=NULL, cell=NULL){

  require('ggplot2')

  # Creates color palette function.
  p_palette <- khroma::colour(color_pal)

  data_f$krige <- data_f$krige/max(data_f$krige)

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
    scale_fill_gradientn(colors=p_palette(5)) +
    xlab("X Position") +
    ylab("Y Position") +
    labs(fill=leg_name) +
    ggtitle(title_name) +
    ggpolypath::geom_polypath(aes(long,lat,group=group), mask_df, fill="white"#,
                              #color='white', size=1.5
                              ) +
    geom_point(data=pvals_sign, aes(x=X2, y=X3), shape=1, alpha=0.25, size=0.5) +
    coord_fixed() +
    theme_classic() +
    theme(legend.position="right", plot.title=element_text(size=8))

  return(p)

}
