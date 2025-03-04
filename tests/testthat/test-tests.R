##
# Unit tests via testthat
#

# Create STlist from Thrane et al. data
# FUTURE DEV: Maybe, would be better to have a package (spatialGEData) to perform
# tests with other ST technologies

# library('spatialGE')
data_files = list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
count_files = grep("counts", data_files, value=T)
coord_files = grep("mapping", data_files, value=T)
clin_file = grep("thrane_clinical", data_files, value=T)
melanoma = STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples

# Test that resulting STlist is an S4 object
testthat::test_that("Data input checks output single character string", {
  testthat::expect_s4_class(melanoma, 'STlist')
})
