# summarize_STlist: Generates a data frame with summary statistics

Produces a data frame with counts per gene and counts per ROI/spot/cell

## Usage

``` r
summarize_STlist(x = NULL)
```

## Arguments

- x:

  an STlist

## Value

a data frame

## Details

The function creates a table with counts per gene and counts per region
of interest (ROI), spot, or cell in the samples stored in the STlist

## Examples

``` r
# \donttest{
# Using included melanoma example (Thrane et al.)
# Download example data set from spatialGE_Data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)
lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
tryCatch({ # In case data is not available from network
  download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
  zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
  unzip(zipfile=zip_tmp, exdir=thrane_tmp)
  # Generate the file paths to be passed to the STlist function
  count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='counts')
  coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='mapping')
  clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                          full.names=TRUE, pattern='clinical')
  # Create STlist
  library('spatialGE')
  melanoma <- STlist(rnacounts=count_files,
                     spotcoords=coord_files,
                     samples=clin_file) # Only first two samples
  summarize_STlist(melanoma)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?.")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?.
#> NULL
# }
```
