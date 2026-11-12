# get_gene_meta: Extract gene-level metadata and statistics

Extracts gene-level metadata and spatial statistics (if already
computed)

## Usage

``` r
get_gene_meta(x = NULL, sthet_only = FALSE)
```

## Arguments

- x:

  an STlist

- sthet_only:

  logical, return only genes with spatial statistics

## Value

a data frame with gene-level data

## Details

This function extracts data from the `x@gene_meta` slot, optionally
subsetting only to those genes for which spatial statistics (Moran's I
or Geary's C, see `SThet`) have been calculated. The output is a data
frame with data from all samples in the STlist

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
                     samples=clin_file)
  melanoma <- transform_data(melanoma)
  melanoma <- SThet(melanoma, genes=c('MLANA', 'TP53'), method='moran')
  get_gene_meta(melanoma, sthet_only=TRUE)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```
