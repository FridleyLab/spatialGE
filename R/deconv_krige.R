##
#' @title deconv_krige: Spatial interpolation of gene expression deconvolution scores
#' @description Performs spatial interpolation ('kriging') of cell scores from
#' gene expression deconvolution scores.
#' @details
#' This function takes a STList and a vector of xCell cell names and performs
#' spatial interpolation ("kriging") of deconvoluted cell scores. It also calculates
#' spatial heterogeneity statistics. The function can perform ordinary or universal kriging.
#' If cells='top', then the 10 cell types with the highest standard deviation for each
#' tissue array are interpolated. The function supports kriging via the Python package
#' PyKrige or geoR. The first option is faster but requires a Python environment to be
#' set up.
#'
#' @param x an STList with transformed xCell scores.
#' @param cells a vector of cell names or 'top'. If 'top' (default), kriging is
#' estimated for the 10 cells with highest standard deviation for each spatial arrays.
#' @param univ a logical stating whether or not to perform universal kriging.
#' Default is FALSE (ordinary kriging).
#' @param res a double to adjust the resolution of the plot. Fractions of 1 lead to
#' more resolution, but longer run times. Default is 0.5 for Visium ST arrays.
#' @param who the spatial arrays for which kriging will be performed. If NULL (Default),
#' all arrays are kriged.
#' @param method the method from which deconvolution scores will be taken from. The
#' default is 'xCell' and the only method supported (for now).
#' @param python a logical, whether or not to use the Python implementation. If FALSE,
#' geoR is used.
#' @return x, an STList including spatial interpolations.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' melanoma <- deconv_krige(melanoma, cells=c('b_cells', 'i_dc'), who=2, python=F)
#'
#' # Running python implementation on the most variable cell types.
#' # melanoma <- deconv_krige(melanoma, cells='top', who=2, python=T)
#'
#' @export
#
#
deconv_krige = function(x=NULL, cells='top', univ=F, res=NULL, who=NULL, method='xcell', python=T){

  require("magrittr")

  # geoR implementation of universal kriging to be implmented. Probably allow users to
  # specify parameters from variogram
  if(python == F && univ == T){
    stop('Currently, universal kriging is only available using Python kriging (PyKrige)')
  }

  # TEMPORARY: This check due to STList getting too heavy on memory after one cell type.
  # if(nrow(x@coords[[1]]) > 1007 && python == F){
  #   if(length(cells) > 1){
  #     stop('For large arrays (e.g. Visium), one cell type at a time can be interpolated.')
  #   }
  # }

  # Detect the deconvolution method requested.
  if(tolower(method) == 'xcell'){
    method = 'xCell'
  } else if(tolower(method) == 'ssgsea'){
    method = 'ssGSEA'
  } else{
    stop('Please, specify a deconvolution method to process.')
  }

  # Test if deconvoluted data are available.
  if (rlang::is_empty(x@cell_deconv[[method]])) {
    stop(paste("There are no deconvolution results in this STList."))
  }

  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who <- c(1:length(x@cell_deconv[[method]]))
  }

  # Test that a cell name was entered.
  if (is.null(cells)) {
    stop("Please, enter one or more cell names to plot.")
  }

  # Get requested list of deconvoluted matrices.
  deconv_list <- x@cell_deconv[[method]]

  # If cells='top', get names of 10 cells with the highest standard deviation.
  if(length(cells) == 1 && cells == 'top' && method == 'xCell'){
    cells = c()
    for(i in who){
      cells = append(cells, deconv_list[[i]]$cell_stdev$cell[order(deconv_list[[i]]$cell_stdev$cell_stdevs, decreasing=T)][1:10])
    }
    # Get unique cell names from most variable.
    cells = unique(cells)
  } else if(length(cells) == 1 && cells == 'top' && method != 'xCell'){
    stop('The token \"top\" is only enabled for xCell.')
  }

  # Test if list with kriging exists for each cell If not create it.
  for(cell in cells){
    if(is.null(x@cell_krige[[cell]]) && rlang::is_empty(x@cell_krige[[cell]])){
      x@cell_krige[[cell]] = list()
      for(i in length(x@counts)){
        x@cell_krige[[cell]][[i]] = list()
      }
    }
  }

  # Loop through each deconvolution matrix.
  for (i in who) {
    # Specify resolution if not input by user. Original ST slide had 1007 spots.
    if(is.null(res)){
      if(nrow(x@coords[[i]]) > 1007){
        res = 0.5
      } else{
        res = 0.2
      }
      x@deconv_krige_data[['res']] = res
    } else if(!is.null(x@deconv_krige_data[['res']])){
      if(x@deconv_krige_data[['res']] != res){
        cat('Kriging results will be overwritten. Previosuly \"kriged\" cell scores will not be available.')
        x@deconv_krige_data[['res']] = res
      }
    }

    # Store kriging type.
    if(univ){
      x@deconv_krige_data[['krige_type']] = 'universal'
    }else{
      x@deconv_krige_data[['krige_type']] = 'ordinary'
    }

    # Create lists to store prediction grids and borders.
    if(is.null(x@deconv_krige_data[['krige_border']])){
      x@deconv_krige_data[['krige_border']] = list()
      x@deconv_krige_data[['krige_grid']] = list()
      for(k in length(x@counts)){
        x@deconv_krige_data[['krige_border']][[k]] = list()
        x@deconv_krige_data[['krige_grid']][[k]] = list()
      }
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@deconv_krige_data[['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@coords[[i]][c(3, 2)]))

    # Get cells present in specific sample.
    subsetcells_mask = cells %in% deconv_list[[i]]$sqrt_scores$cell_names
    subsetcells = cells[subsetcells_mask]

    # Get cells not present.
    notcells = cells[!subsetcells_mask]

    if(!rlang::is_empty(notcells)){
      cat(paste(paste(notcells, collapse=', '), ": Not present in the deconvoluted data for subject", i, "."))
    }

    # Create grid for PyKrige or geoR.
    if(python == T){
      gridx = seq(min(x@coords[[i]][[3]]), max(x@coords[[i]][[3]]), res)
      gridy = seq(min(x@coords[[i]][[2]]), max(x@coords[[i]][[2]]), res)
      gridx = gridx[-length(gridx)]
      gridy = gridy[-length(gridy)]
      cell_geo_grid = expand.grid(
        seq(min(gridx), max(gridx), by=res),
        seq(min(gridy), max(gridy), by=res)
      )
      # Store prediction grid and kriging algorithm in STList.
      x@deconv_krige_data[['algorithm']] = 'pykrige'
      #x@gene_krige_data[['krige_grid']] = list()
      x@deconv_krige_data[['krige_grid']][[i]] = cell_geo_grid
    } else if(python == F && univ == F){ #  Create grid for geoR estimation.
      cell_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res),
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res)
      )
      # Store prediction grid and kriging algorithm in STList.
      x@deconv_krige_data[['algorithm']] = 'geor'
      #x@gene_krige_data[['krige_grid']] = list()
      x@deconv_krige_data[['krige_grid']][[i]] = cell_geo_grid
    }

    # Define number of available cores to use.
    cores = 1
    if(.Platform$OS.type == 'unix'){
      avail_cores = parallel::detectCores()
      if(avail_cores > (length(subsetcells) + 1)){
        cores = (length(subsetcells) + 1)
      } else if( (avail_cores <= (length(subsetcells) + 1)) && avail_cores > 1){
        cores = avail_cores - 1
      }
    }

    # Prepare inputs for kriging.
    # Loop through cells
    kriging_list = parallel::mclapply(seq_along(subsetcells), function(j){
      # Get transformed scores.
      cell_abund = deconv_list[[i]]$sqrt_scores[deconv_list[[i]]$sqrt_scores[[1]] == subsetcells[j], -1]
      cell_abund = as.data.frame(t(cell_abund))
      cell_abund = cell_abund %>%
        tibble::rownames_to_column(., var="position")
      # Match order of library names in counts data and coordinate data.
      cell_abund = cell_abund[match(x@coords[[i]][[1]], cell_abund[[1]]), ]
      cell_geo_df = cbind(x@coords[[i]][c(3, 2)], as.numeric(cell_abund[[2]]))
      colnames(cell_geo_df)[3] = "cell_data"
      # Call the requested Kriging algorithm
      if(python == T){
        #gridx_red = gridx[-length(gridx)]
        #gridy_red = gridy[-length(gridy)]
        # Call PyKrige implementation.
        kriging_res = krige_py(gridx=gridx, gridy=gridy, geo_df=cell_geo_df, univ=univ)
      } else{
        # Create geodata object from expression and coordinate data
        cell_geo <- geoR::as.geodata(cell_geo_df, coords.col=c(1,2), data.col=3)
        kriging_res = krige_geor(geodata=cell_geo, locations=x@deconv_krige_data[['krige_grid']][[i]], univ=univ)
      }
      return(kriging_res)
    }, mc.cores=cores, mc.preschedule=T)

    # Store kriging results in STList.
    for(j in 1:length(subsetcells)){
      x@cell_krige[[subsetcells[j]]][[i]] = kriging_list[[j]]
    }

    # Test if list with spatial statistics exists for each gene. If not create it.
    for(cell in cells){
      if(is.null(x@cell_het[[cell]]) && rlang::is_empty(x@cell_het[[cell]])){
        x@cell_het[[cell]] = list()
        for(k in length(x@counts)){
          x@cell_het[[cell]][[k]] = list()
        }
      }
    }

    # Estimate spatial heterogeneity statistics.
    # moran_list = parallel::mclapply(seq_along(subsetcells), function(j){
    #   temp_x = cell_moran_I(x, cells=subsetcells[j], subj=i)
    #   temp_x = temp_x@cell_het[[subsetcells[j]]][[i]]$morans_I
    #   return(temp_x)
    # }, mc.cores=cores, mc.preschedule=T)
    #
    # geary_list = parallel::mclapply(seq_along(subsetcells), function(j){
    #   temp_x = cell_geary_C(x, cells=subsetcells[j], subj=i)
    #   temp_x = temp_x@cell_het[[subsetcells[j]]][[i]]$gearys_C
    #   return(temp_x)
    # }, mc.cores=cores, mc.preschedule=T)
    #
    # getis_list = parallel::mclapply(seq_along(subsetcells), function(j){
    #   temp_x = cell_getis_Gi(x, cells=subsetcells[j], subj=i)
    #   temp_x = temp_x@cell_het[[subsetcells[j]]][[i]]$getis_ord_Gi
    #   return(temp_x)
    # }, mc.cores=cores, mc.preschedule=T)
    #
    # # Store spatial statistics in STList.
    # for(j in 1:length(subsetcells)){
    #   x@cell_het[[subsetcells[j]]][[i]]$morans_I = moran_list[[j]]
    #   x@cell_het[[subsetcells[j]]][[i]]$gearys_C = geary_list[[j]]
    #   x@cell_het[[subsetcells[j]]][[i]]$getis_ord_Gi = getis_list[[j]]
    # }
  }
  return(x)
}
