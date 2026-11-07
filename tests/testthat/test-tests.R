##
# Unit tests via testthat
#

# Create STlist from Thrane et al. data
testthat::test_that("STlist() builds correctly from the Thrane et al. network data", {
  testthat::skip_on_cran()
  testthat::skip_if_offline()

  thrane_tmp = tempdir()
  unlink(thrane_tmp, recursive=TRUE)
  dir.create(thrane_tmp)

  lk = 'https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
  zip_tmp = tryCatch({ # In case data is not available from network
    download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb', quiet=TRUE)
    list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
  }, error = function(e) character(0))

  testthat::skip_if(length(zip_tmp) == 0, "Could not download Thrane et al. test data.")

  unzip(zipfile=zip_tmp, exdir=thrane_tmp)
  count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='counts')
  coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='mapping')
  clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'), full.names=TRUE, pattern='clinical')
  melanoma = STlist(rnacounts=count_files[c(1,2)], spotcoords=coord_files[c(1,2)], samples=clin_file) # Only first two samples

  # Test that resulting STlist is an S4 object
  testthat::expect_s4_class(melanoma, 'STlist')
})

