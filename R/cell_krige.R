# This function performs spatial interpolation of deconvoluted cell scores. This
# function takes a STList and a list of cell names or the token 'top' for the 5
# cells with the highest variation. It also calculates spatial heterogeneity
# measures for the genes. The function can perform ordinary or universal kriging.
# The result can be plotted using the plot_cell_krige() function.
#
# @param x, a STList with normalized counts
# @plot_who, a vector of cell names or 'top'. If 'top', kriging for the 5 cells
# with highest standard deviation is estimated.
# @univ, a logical stating whether or not to perform universal or ordinary kriging.
# @res, a number to adjust the resolution of the plot. Fractions of 1 lead to
# more resolution.
# @return x, a STList including an spatial interpolation object.
#
#
require('concaveman')
require('geoR')
require('RColorBrewer')

cell_krige <- function(x=NULL, plot_who='top', univ=F, res=0.1){

  # If plot_who='top', get names of 5 cells with the highest standard deviation.
  if(length(plot_who) == 1){
    if(plot_who == 'top'){
      cells <- x@cell_stdev$cell[order(x@cell_stdev$cell_stdevs, decreasing=T)][1:5]
    }else{
      cells <- plot_who
    }
  }else{
    cells <- plot_who
  }

  # Loop through cells.
  for(cell in cells){

    # Test that the cell name is present in the deconvolution matrix.
    if(!any(x@cell_deconv$deconv_matrix[[1]] == cell)){
      stop(paste(cell, "is not a cell in the", x@cell_deconv$deconv_method,
                 "deconvolution matrix."))
    }

    # Extract abundance/score data for a given cell.
    cell_abund <- x@cell_deconv$deconv_matrix[
      x@cell_deconv$deconv_matrix[[1]] == cell, -1]

    # Transpose abundance/score data to turn it into a column. Then turn library
    # names into a column and assign column names (first row).
    cell_abund <- as.data.frame(t(cell_abund))
    cell_abund <- cell_abund %>% rownames_to_column(., var='position')
    #colnames(gene_expr) <- gene_expr[1,]
    colnames(cell_abund)[2] <- 'cell_data'
    #gene_expr <- as.data.frame(gene_expr[-1,])

    # Sort cell data using the order in the mapping file. Then add
    # coordinates data to expression data frame.
    cell_abund <- cell_abund[match(x@coords[[1]], cell_abund[[1]]), ]
    cell_geo_df <- cbind(x@coords[2:3], as.numeric(cell_abund[[2]]))

    # Create concave hull to use as delimiter of sampled area. Needs to be
    # done before converting data frame to spatial object.
    conc_hull <- concaveman(as.matrix(cell_geo_df[1:2]))

    # Create geodata object from cell and coordinate data
    cell_geo <- as.geodata(cell_geo_df, coords.col=c(1,2), data.col=3)

    # Create a grid finer than the sampled locations to predict locations.
    cell_geo_grid <-expand.grid(
      seq((min(x@coords[[2]])-1), (max(x@coords[[2]])+1), by=res),
      seq((min(x@coords[[3]])-1), (max(x@coords[[3]])+1), by=res)
    )

    # Add concave hull to geodata.
    cell_geo$borders <- conc_hull

    # OC <- output.control(simulations=TRUE, n.pred=10,
    #                      quantile=c(0.1, 0.25, 0.5, 0.75, 0.9),
    #                      threshold = 350)

    # Test if slot for gene kriging is already present. Else, create it.
    if(is.null(x@cell_krige[[cell]])){
      x@cell_krige[[cell]] <- list(ord=NULL,
                                   univ=NULL)
    }

    # Create controls for either ordinary or universal kriging, and perform
    # estimation.
    if(univ == F){
      # NOTE: Need to check how to decide on init.cov.pars
      cell_geo_lhood <- likfit(cell_geo, trend='cte', ini.cov.pars=c(1, 0.15))

      # Specify control (and output) parameters for ordinary kriging.
      KC <- krige.control(obj.model=cell_geo_lhood)

      # Perform ordinary kriging.
      cell_krig <- krige.conv(cell_geo, locations=cell_geo_grid, krige=KC)

      x@cell_krige[[cell]]$ord <- cell_krig

    }
    else{
      # NOTE: Need to use regression analysis of variogram to get values for
      # nugget.
      cell_geo_lhood <- likfit(cell_geo, trend='cte',
                               ini.cov.pars=c(1000, 500), nug=100)

      # Specify control (and output) parameters for universal kriging.
      KC <- krige.control(type.krige="OK", obj.m=cell_geo_lhood,
                          trend.d="cte",
                          trend.l="cte")

      # Perform universal kriging.
      cell_krig <- krige.conv(cell_geo, locations=cell_geo_grid, krige=KC
                                   #output = OC
                                   )

      x@cell_krige[[cell]]$univ <- cell_krig
    }

    # Calculate spatial heterogeneity statistics.
    x <- cell_moran_I(x, cells=cell)
    x <- cell_geary_C(x, cells=cell)
    x <- cell_getis_Gi(x, cells=cell)

  }

  return(x)

}
