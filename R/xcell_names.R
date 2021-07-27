##
#' @title xcell_names: Prints xCell cell type names
#' @description Prints the available cell types as provided by xCell deconvolution.
#' @details
#' Shows the names of cell types that can be used in other spatialGE functions.
#' This function access the xCell cell type names after cleaning via the janitor
#' package.
#'
#' @return The xCell cell type names after cleaning.
#'
#' @examples
#' # xcell_names()
#'
#' @export
#'
#'
xcell_names = function(){
  require('magrittr')

  xcellnames = x@cell_deconv$xCell[[1]]$cell_stdev$cell

  print(xcellnames)
}
