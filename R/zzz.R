##
#
#

.onLoad = function (libname, pkgname) {
  ns = topenv()
  ns$datafile = system.file("extdata/melanoma_thrane", package="spatialGE")
}
