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
#' # melanoma <- deconv_krige(melanoma, cells=c('b_cells', 'i_dc'), who=2, python=F)
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
      cells = append(cells, deconv_list[[i]]$cell_var$cell[order(deconv_list[[i]]$cell_var$cell_stdevs, decreasing=T)][1:10])
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
      for(i in 1:length(x@tr_counts)){
        x@cell_krige[[cell]][[names(x@tr_counts[i])]] = list()
      }
    }
  }

  # Store kriging type.
  if(univ){
    x@misc[['cell_krige_type']] = 'universal'
  }else{
    x@misc[['cell_krige_type']] = 'ordinary'
  }

  # Specify resolution if not input by user
  if(is.null(res)){
    res = 0.5
  }
  x@misc[['cell_krige_res']] = res

  # Give warning about kriging res < 0.5 for large matrices (e.g. Visium)
  if(res < 0.5 & (nrow(x@coords[[1]]) > 1000)){
    cat('Kriging at the requested resolution might take some time.\n')
  }

  # Create lists to store prediction grids and borders.
  if(is.null(x@misc[['krige_border']])){
    x@misc[['krige_border']] = list()
    x@misc[['cell_krige_grid']] = list()
    for(k in length(x@tr_counts)){
      x@misc[['krige_border']][[k]] = list()
      x@misc[['cell_krige_grid']][[k]] = list()
    }
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in who){
    subsetcells_mask = cells %in% deconv_list[[i]][['sqrt_scores']][['cell_names']]
    subsetcells = cells[subsetcells_mask]
    combo = dplyr::bind_rows(combo, expand.grid(names(x@tr_counts[i]), subsetcells))

    # Get genes not present.
    notcells = cells[!subsetcells_mask]

    if(!rlang::is_empty(notcells)){
      cat(paste(paste(notcells, collapse=', '), ": Not present in the transformed counts for sample ", names(x@tr_counts[i]), ".\n"))
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@misc[['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@coords[[i]][c(3, 2)]))

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
      x@misc[['cell_krige_grid']][[i]] = cell_geo_grid
    } else if(python == F && univ == F){ #  Create grid for geoR estimation.
      cell_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res),
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res)
      )
      x@misc[['cell_krige_grid']][[i]] = cell_geo_grid
    }
  }

  # Store prediction grid and kriging algorithm in STList.
  if(python){
    x@misc[['cell_krige_algorithm']] = 'pykrige'
  } else{
    x@misc[['cell_krige_algorithm']] = 'geor'
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  kriging_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get transformed counts.
    cell_score = deconv_list[[i]][['sqrt_scores']][deconv_list[[i]][['sqrt_scores']][['cell_names']] == j, -1]
    cell_score = as.data.frame(t(cell_score))
    cell_score = cell_score %>%
      tibble::rownames_to_column(., var="position")
    # Match order of library names in counts data and coordinate data.
    cell_score = cell_score[match(x@coords[[i]][[1]], cell_score[[1]]), ]
    cell_geo_df = cbind(x@coords[[i]][c(3, 2)], as.numeric(cell_score[[2]]))
    colnames(cell_geo_df)[3] = "cell_score"
    # Call the requested Kriging algorithm
    if(python == T){
      # Call PyKrige implementation.
      gridx = seq(min(x@coords[[i]][[3]]), max(x@coords[[i]][[3]]), res)
      gridy = seq(min(x@coords[[i]][[2]]), max(x@coords[[i]][[2]]), res)
      gridx = gridx[-length(gridx)]
      gridy = gridy[-length(gridy)]
      kriging_res = krige_py(gridx=gridx, gridy=gridy, geo_df=cell_geo_df, univ=univ)
    } else{
      # Create geodata object from expression and coordinate data
      cell_geo <- geoR::as.geodata(cell_geo_df, coords.col=c(1,2), data.col=3)
      kriging_res = krige_geor(geodata=cell_geo, locations=x@misc[['cell_krige_grid']][[grep(i, names(deconv_list))]], univ=univ)
    }
    return(kriging_res)
  }, mc.cores=cores, mc.preschedule=F)
  names(kriging_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(kriging_list)[i], split = '&&'))
    x@cell_krige[[combo_name[2]]][[combo_name[1]]] = kriging_list[[i]]
  }

  return(x)
}

