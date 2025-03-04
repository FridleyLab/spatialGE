##
#
#

# .onLoad = function (libname, pkgname) {
#   ns = topenv()
#   ns$datafile = system.file("inst/extdata/melanoma_thrane", package="spatialGE")
# }


.onLoad = function(libname, pkgname) {
  # Get the full path to  data file within the package
  thrane_files = system.file("inst/extdata/melanoma_thrane", package=pkgname)

  # Load the data using the path
  #data(file=thrane_files, envir=environment())
}
