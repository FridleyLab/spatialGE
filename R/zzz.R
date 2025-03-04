##
#
#

# .onLoad = function (libname, pkgname) {
#   ns = topenv()
#   ns$datafile = system.file("inst/extdata/melanoma_thrane", package="spatialGE")
# }


thrane_files = function() {
  #ns = topenv()
  # Get the full path to  data file within the package
  #ns$thrane_files = system.file("inst/extdata/melanoma_thrane", package='spatialGE')
  thrane_files = system.file("extdata/melanoma_thrane", package='spatialGE')
  return(thrane_files)
  # Load the data using the path
  #data(file=thrane_files, envir=environment())
}
