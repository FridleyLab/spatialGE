##
# Unit tests via testthat
#

# Create STlist from Thrane et al. data
# FUTURE DEV: Maybe, would be better to have a package (spatialGEData) to perform
# tests with other ST technologies

# library('spatialGE')
data_files <- system.file("extdata", 'melanoma_thrane', package="spatialGE")
count_files <- list.files(data_files, full.names=T, pattern='counts')
coord_files <- list.files(data_files, full.names=T, pattern='mapping')
clin_file <- list.files(data_files, full.names=T, pattern='clinical')
melanoma = STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples

# Test that resulting STlist is an S4 object
testthat::test_that("Data input checks output single character string", {
  testthat::expect_s4_class(melanoma, 'STlist')
})
