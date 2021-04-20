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
require('rlang')
cell_krige <- function(x=NULL, cells='top', univ=F, res=0.2, who=NULL){

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who <- c(1:length(x@cell_deconv))
  }

  # Test if deconvoluted data are available.
  if (is_empty(x@cell_deconv)) {
    stop(paste("There are no deconvolution results in this STList."))
  }

  # Loop through each deconvolution results table.
  for (i in who) {

    # If cells='top', get names of 10 cell types with the highest standard deviation.
    if(length(cells) == 1){
      if(cells == 'top'){
        cells <- x@cell_stdev[[i]]$cell[order(x@cell_stdev[[i]]$cell_stdevs, decreasing=T)][1:10]
      }
    }

    # Append "stromal score"
    cells <- append(cells, 'stroma_score')

    # Loop through cells.
    for(cell in cells){

      # Test that the cell name is present in the deconvolution matrix.
      if(!any(x@cell_deconv[[i]]$deconv_matrix[[1]] == cell)){
        if(cell != 'stroma_score'){
          cat(paste(cell, "is not a cell in the", x@cell_deconv$deconv_method,
                    "deconvolution matrix."))
          next
        }
      }

      # Test if slot for cell kriging is already present. Else, create it.
      if(!is.null(x@cell_krige[[cell]])){
        if(length(x@cell_krige[[cell]]) < i){
          x@cell_krige[[cell]][[i]] <- list(ord=NULL,
                                            univ=NULL)
        }else{
          next
        }
      }else{
        x@cell_krige[[cell]] <- list()
        x@cell_krige[[cell]][[i]] <- list(ord=NULL,
                                          univ=NULL)
      }

      # If 'stroma score' is being kriged, then get data from corresponding slot.
      if(cell != 'stroma_score'){
        # Extract abundance/score data for a given cell.
        cell_abund <- x@cell_deconv[[i]]$transf_deconv_matrix[
          x@cell_deconv[[i]]$transf_deconv_matrix[[1]] == cell, -1]
      }else{
        cell_abund <- x@cell_deconv[[i]]$transf_tumorstroma[
          x@cell_deconv[[i]]$transf_tumorstroma[[1]] == cell, -1]
        cell_abund <- sqrt(cell_abund)
      }

      # Transpose abundance/score data to turn it into a column. Then turn library
      # names into a column and assign column names (first row).
      cell_abund <- as.data.frame(t(cell_abund))
      cell_abund <- cell_abund %>% tibble::rownames_to_column(., var='position')
      #colnames(gene_expr) <- gene_expr[1,]
      colnames(cell_abund)[2] <- 'cell_data'
      #gene_expr <- as.data.frame(gene_expr[-1,])

      # Sort cell data using the order in the mapping file. Then add
      # coordinates data to expression data frame.
      cell_abund <- cell_abund[match(x@coords[[i]][[1]], cell_abund[[1]]), ]
      cell_geo_df <- cbind(x@coords[[i]][2:3], as.numeric(cell_abund[[2]]))

      # Create concave hull to use as delimiter of sampled area. Needs to be
      # done before converting data frame to spatial object.
      x@prediction_border[[i]] <- concaveman::concaveman(as.matrix(cell_geo_df[1:2]))

      # Create geodata object from cell and coordinate data
      cell_geo <- geoR::as.geodata(cell_geo_df, coords.col=c(1,2), data.col=3)

      # Create a grid finer than the sampled locations to predict locations.
      cell_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res),
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res)
      )

      # Store prediction grid in STList.
      if(univ){
        x@cell_krige[[cell]][[i]]$univ_grid <- cell_geo_grid
      }else{
        x@cell_krige[[cell]][[i]]$ord_grid <- cell_geo_grid
      }

      # OC <- output.control(simulations=TRUE, n.pred=10,
      #                      quantile=c(0.1, 0.25, 0.5, 0.75, 0.9),
      #                      threshold = 350)

      # Create controls for either ordinary or universal kriging, and perform
      # estimation.
      if(univ == F){
        # NOTE: Need to check how to decide on init.cov.pars
        cell_geo_lhood <- geoR::likfit(cell_geo, trend='cte', ini.cov.pars=c(1, 0.15))

        # Specify control (and output) parameters for ordinary kriging.
        KC <- geoR::krige.control(obj.model=cell_geo_lhood)

        # Perform ordinary kriging.
        cell_krig <- geoR::krige.conv(cell_geo, locations=cell_geo_grid, krige=KC)

        x@cell_krige[[cell]][[i]][['ord']] <- cell_krig

      }else if(univ == T){
        # NOTE: Need to use regression analysis of variogram to get values for
        # nugget.
        cell_geo_lhood <- geoR::likfit(cell_geo, trend='cte',
                                 ini.cov.pars=c(1000, 500), nug=100)

        # Specify control (and output) parameters for universal kriging.
        KC <- geoR::krige.control(type.krige="OK", obj.m=cell_geo_lhood,
                            trend.d="cte",
                            trend.l="cte")

        # Perform universal kriging.
        cell_krig <- geoR::krige.conv(cell_geo, locations=cell_geo_grid, krige=KC
                                #output = OC
        )

        x@cell_krige[[cell]][[i]][['ord']] <- cell_krig
      }

      if(cell != 'stroma_score'){
        # Calculate spatial heterogeneity statistics.
        x <- cell_moran_I(x, cells=cell, subj=i)
        x <- cell_geary_C(x, cells=cell, subj=i)
        x <- cell_getis_Gi(x, cells=cell, subj=i)
      }

    }

  }

  return(x)

}
