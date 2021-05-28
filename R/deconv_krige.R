##
#' @title deconv_krige: Spatial interpolation of gene expression deconvolution scores
#' @description Performs spatial interpolation ('kriging') of cell scores from
#' gene expression deconvolution in spatially-resolved transcriptomics data.
#' @details
#' This function takes a STList and a vector of xCell cell names, or the token 'top' for
#' the 10 cells with the highest standard deviation. It also calculates spatial heterogeneity
#' statistics for the cell scores. The function can perform ordinary or universal kriging.
#' The result can be plotted using the plot_deconv_krige() function.
#'
#' @param x, an STList with transformed xCell scores.
#' @param cells, a vector of cell names or 'top'. If 'top' (default), kriging for the 10
#' cells with highest standard deviation is estimated.
#' @param univ, a logical stating whether or not to perform universal or ordinary kriging.
#' Default is FALSE (ordinary kriging).
#' @param res, a double to adjust the resolution of the plot. Fractions of 1 lead to
#' more resolution. Default is res=0.2.
#' @param who, the spatial arrays for which kriging will be performed. If NULL (Default),
#' all arrays are kriged.
#' @return x, an STList including spatial interpolations.
#' @export
#
#
deconv_krige <- function(x=NULL, cells='top', univ=F, res=0.2, who=NULL){

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who <- c(1:length(x@cell_deconv$xCell))
  }

  # Test if deconvoluted data are available.
  if (rlang::is_empty(x@cell_deconv)) {
    stop(paste("There are no deconvolution results in this STList."))
  }

  # Loop through each deconvolution matrix.
  for (i in who) {

    # If cells='top', get names of 10 cell types with the highest standard deviation.
    if(length(cells) == 1){
      if(cells == 'top'){
        cells <- x@cell_deconv$xCell[[i]]$cell_stdev$cell[order(x@cell_deconv$xCell[[i]]$cell_stdev, decreasing=T)][1:10]
      }
    }

    # Loop through cells.
    for(cell in cells){

      # Test that the cell name is present in the deconvolution matrix.
      if(!any(x@cell_deconv$xCell[[i]]$sqrt_scores[[1]] == cell)){
        cat(paste(cell, "is not a cell in the", names(x@cell_deconv)[2], "deconvolution matrix."))
        next
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

      # Extract abundance/score data for a given cell.
      cell_abund <- x@cell_deconv$xCell[[i]]$sqrt_scores[x@cell_deconv$xCell[[i]]$sqrt_scores[[1]] == cell, -1]

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
      cell_geo_df <- cbind(x@coords[[i]][c(3,2)], as.numeric(cell_abund[[2]]))

      # Create concave hull to use as delimiter of sampled area. Needs to be
      # done before converting data frame to spatial object.
      x@prediction_border[[i]] <- concaveman::concaveman(as.matrix(cell_geo_df[1:2]))

      # Create geodata object from cell and coordinate data
      cell_geo <- geoR::as.geodata(cell_geo_df, coords.col=c(1,2), data.col=3)

      # Create a grid finer than the sampled locations to predict locations.
      cell_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res),
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res)
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

      # Provide info progress.
      cat(paste0("Performing spatial interpolation ('kriging') of ", cell, " for subject ", i, "...\n"))

      # Create controls for either ordinary or universal kriging, and perform
      # estimation.
      if(univ == F){
        # NOTE: Need to check how to decide on init.cov.pars
        cell_geo_lhood <- geoR::likfit(cell_geo, trend='cte',
                                       ini.cov.pars=c(1, 0.15),
                                       messages=F)

        # Specify control (and output) parameters for ordinary kriging.
        KC <- geoR::krige.control(obj.model=cell_geo_lhood)
        OC <- geoR::output.control(messages=F)

        # Perform ordinary kriging.
        cell_krig <- geoR::krige.conv(cell_geo, locations=cell_geo_grid,
                                      krige=KC, output=OC)

        x@cell_krige[[cell]][[i]][['ord']] <- cell_krig

      }else if(univ == T){
        # NOTE: Need to use regression analysis of variogram to get values for
        # nugget.
        cell_geo_lhood <- geoR::likfit(cell_geo, trend='cte',
                                 ini.cov.pars=c(1000, 500), nug=100,
                                 messages=F)

        # Specify control (and output) parameters for universal kriging.
        KC <- geoR::krige.control(type.krige="OK", obj.m=cell_geo_lhood,
                            trend.d="cte",
                            trend.l="cte")
        OC <- geoR::output.control(messages=F)

        # Perform universal kriging.
        cell_krig <- geoR::krige.conv(cell_geo, locations=cell_geo_grid,
                                      krige=KC, output=OC)

        x@cell_krige[[cell]][[i]][['ord']] <- cell_krig
      }

        # Calculate spatial heterogeneity statistics.
      x <- cell_moran_I(x, cells=cell, subj=i)
      x <- cell_geary_C(x, cells=cell, subj=i)
      x <- cell_getis_Gi(x, cells=cell, subj=i)

    }

  }

  return(x)

}
