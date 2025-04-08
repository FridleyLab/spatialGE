##
# Unit tests via testthat
#

# Create STlist from Thrane et al. data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)
download.file('https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=',
			  destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
unzip(zipfile=zip_tmp, exdir=thrane_tmp)
count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='counts')
coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='mapping')
clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='clinical')
# library('spatialGE')
melanoma = STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples

# Test that resulting STlist is an S4 object
testthat::test_that("Data input checks output single character string", {
  testthat::expect_s4_class(melanoma, 'STlist')
})
