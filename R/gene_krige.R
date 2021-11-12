##
#' @title gene_krige: Spatial interpolation of gene expression
#' @description Performs spatial interpolation ('kriging') of transformed gene
#' counts in spatially-resolved transcriptomics data.
#' @details
#' This function takes a STList and a vector of gene names and generates spatial
#' interpolation of gene expression values. It also calculates spatial heterogeneity
#' statistics. The function can perform ordinary or universal kriging. If genes='top',
#' then the 10 genes  with the highest standard deviation for each tissue array are
#' interpolated. The function supports kriging via the Python package PyKrige or
#' geoR. The first option is faster but requires a Python environment to be set up.
#'
#' @param x an STList with transformed RNA counts.
#' @param genes a vector of HUGO names or 'top'. If 'top' (default), kriging of
#' the 10 genes with highest standard deviation in each spatial array is estimated.
#' @param univ a logical stating whether or not to perform universal kriging.
#' Default is FALSE (ordinary kriging).
#' @param res a double to adjust the resolution of the plot. Fractions of 1 lead to
#' more resolution, but longer run times. Default is 0.5 for Visium ST arrays.
#' @param who, the spatial arrays for which kriging will be performed. If NULL (Default),
#' all arrays are kriged.
#' @param python a logical, whether or not to use the Python implementation. If FALSE,
#' geoR is used.
#' @return x, a STList including spatial interpolations.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # Kriging of spatial arrays 2 and 3 for 2 genes.
#' melanoma <- gene_krige(melanoma, genes=c('CD74', 'SOX10'), who=c(2, 3), python=F)
#'
#' #' # Kriging of spatial array 1 for the 10 most variable genes and using the Python implementation.
#' # melanoma <- gene_krige(melanoma, genes='top', who=1, python=T)
#'
#' @export
#
#
gene_krige = function(x=NULL, genes='top', univ=F, res=NULL, who=NULL, python=T){

  require("magrittr")

  # geoR implementation of universal kriging to be implmented. Probably allow users to
  # specify parameters from variogram
  if(python == F && univ == T){
    stop('Currently, universal kriging is only available using Python kriging (PyKrige)')
  }

  # TEMPORARY: This check due to STList getting too heavy on memory after one gene.
  # if(nrow(x@coords[[1]]) > 1007 && python == F){
  #   if(length(genes) > 1){
  #     stop('For large arrays (e.g. Visium), one gene at a time can be interpolated.')
  #   }
  # }

  # Test that transformed counts are available
  if(rlang::is_empty(x@tr_counts)) {
    stop("There are not transformed counts in this STList.")
  }

  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who = c(1:length(x@tr_counts))
  }

  # Test that a gene name was entered.
  if(is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1 && genes == 'top'){
    genes = c()
    for(i in who){
      genes = append(genes, x@gene_var[[i]]$gene[order(x@gene_var[[i]]$gene_stdevs, decreasing=T)][1:10])
    }
    # Get unique genes from most variable.
    genes = unique(genes)
  }

  # Test if list with kriging exists for each gene. If not create it.
  for(gene in genes){
    if(is.null(x@gene_krige[[gene]]) && rlang::is_empty(x@gene_krige[[gene]])){
      x@gene_krige[[gene]] = list()
      for(i in 1:length(x@tr_counts)){
        x@gene_krige[[gene]][[names(x@tr_counts[i])]] = list()
      }
    }
  }

  # Store kriging type.
  if(univ){
    x@misc[['gene_krige_type']] = 'universal'
  }else{
    x@misc[['gene_krige_type']] = 'ordinary'
  }

  # Specify resolution if not input by user
  if(is.null(res)){
    res = 0.5
  }
  x@misc[['gene_krige_res']] = res

  # Give warning about kriging res < 0.5 for large matrices (e.g. Visium)
  if(res < 0.5 & (nrow(x@coords[[1]]) > 1000)){
    cat('Kriging at the requested resolution might take some time.\n')
  }

  # Create lists to store prediction grids and borders.
  if(is.null(x@misc[['krige_border']])){
    x@misc[['krige_border']] = list()
    x@misc[['gene_krige_grid']] = list()
    for(k in length(x@tr_counts)){
      x@misc[['krige_border']][[k]] = list()
      x@misc[['gene_krige_grid']][[k]] = list()
    }
  }

  # Generate combination of sample x gene to for.
  combo = tibble::tibble()
  for(i in who){
    subsetgenes_mask = genes %in% x@tr_counts[[i]]$gene
    subsetgenes = genes[subsetgenes_mask]
    combo = dplyr::bind_rows(combo, expand.grid(names(x@tr_counts[i]), subsetgenes))

    # Get genes not present.
    notgenes = genes[!subsetgenes_mask]

    if(!rlang::is_empty(notgenes)){
      cat(paste(paste(notgenes, collapse=', '), ": Not present in the transformed counts for sample ", names(x@tr_counts[i]), ".\n"))
    }

    # Create concave hull to "cookie cut" border kriging surface.
    x@misc[['krige_border']][[i]] = concaveman::concaveman(as.matrix(x@coords[[i]][c(3, 2)]))

    # Create grid for PyKrige or geoR.
    if(python == T){
      gridx = seq(min(x@coords[[i]][[3]]), max(x@coords[[i]][[3]]), res)
      gridy = seq(min(x@coords[[i]][[2]]), max(x@coords[[i]][[2]]), res)
      gridx = gridx[-length(gridx)]
      gridy = gridy[-length(gridy)]
      gene_geo_grid = expand.grid(
        seq(min(gridx), max(gridx), by=res),
        seq(min(gridy), max(gridy), by=res)
      )
      x@misc[['gene_krige_grid']][[i]] = gene_geo_grid
    } else if(python == F && univ == F){ #  Create grid for geoR estimation.
      gene_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res),
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res)
      )
      x@misc[['gene_krige_grid']][[i]] = gene_geo_grid
    }
  }

  # Store prediction grid and kriging algorithm in STList.
  if(python){
    x@misc[['gene_krige_algorithm']] = 'pykrige'
  } else{
    x@misc[['gene_krige_algorithm']] = 'geor'
  }

  # Store decompressed transformed counts list in separate object
  tr_counts = x@tr_counts
  for(mtx in who){
    tr_counts[[mtx]] = expandSparse(x@tr_counts[[mtx]])
  }

  # Define cores available
  cores = count_cores(nrow(combo))
  # Loop through combinations of samples x genes
  kriging_list = parallel::mclapply(seq_along(1:nrow(combo)), function(i_combo){
    i = as.vector(unlist(combo[i_combo, 1]))
    j = as.vector(unlist(combo[i_combo, 2]))

    # Get transformed counts.
    gene_expr = tr_counts[[i]][tr_counts[[i]]$gene == j, -1]
    gene_expr = as.data.frame(t(gene_expr))
    gene_expr = gene_expr %>%
      tibble::rownames_to_column(., var="position")
    # Match order of library names in counts data and coordinate data.
    gene_expr = gene_expr[match(x@coords[[i]][[1]], gene_expr[[1]]), ]
    gene_geo_df = cbind(x@coords[[i]][c(3, 2)], as.numeric(gene_expr[[2]]))
    colnames(gene_geo_df)[3] = "gene_expr"
    # Call the requested Kriging algorithm
    if(python == T){
      # Call PyKrige implementation.
      gridx = seq(min(x@coords[[i]][[3]]), max(x@coords[[i]][[3]]), res)
      gridy = seq(min(x@coords[[i]][[2]]), max(x@coords[[i]][[2]]), res)
      gridx = gridx[-length(gridx)]
      gridy = gridy[-length(gridy)]
      kriging_res = krige_py(gridx=gridx, gridy=gridy, geo_df=gene_geo_df, univ=univ)
    } else{
      # Create geodata object from expression and coordinate data
      gene_geo <- geoR::as.geodata(gene_geo_df, coords.col=c(1,2), data.col=3)
      kriging_res = krige_geor(geodata=gene_geo, locations=x@misc[['gene_krige_grid']][[i]], univ=univ)
    }
    return(kriging_res)
  }, mc.cores=cores, mc.preschedule=F)
  names(kriging_list) = paste(combo[[1]], combo[[2]], sep='&&')

  # Store kriging results in STList.
  for(i in 1:nrow(combo)){
    combo_name = unlist(strsplit(names(kriging_list)[i], split = '&&'))
    x@gene_krige[[combo_name[2]]][[combo_name[1]]] = kriging_list[[i]]
  }

  return(x)
}

