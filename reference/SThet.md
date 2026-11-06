# SThet: Computes global spatial autocorrelation statistics on gene expression

Computes the global spatial autocorrelation statistics Moran's I and/or
Geary's C for a set of genes

## Usage

``` r
SThet(
  x = NULL,
  genes = NULL,
  samples = NULL,
  method = "moran",
  k = NULL,
  overwrite = TRUE,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  an STlist

- genes:

  a vector of gene names to compute statistics

- samples:

  the samples to compute statistics

- method:

  The spatial statistic(s) to estimate. It can be set to 'moran',
  'geary' or both. Default is 'moran'

- k:

  the number of neighbors to estimate weights. By default NULL, meaning
  that spatial weights will be estimated from Euclidean distances. If an
  positive integer is entered, then the faster k nearest-neighbors
  approach is used. Please keep in mind that estimates are not as
  accurate as when using the default distance-based method.

- overwrite:

  logical indicating if previous statistics should be overwritten.
  Default to FALSE (do not overwrite)

- cores:

  integer indicating the number of cores to use during parallelization.
  If NULL, the function uses half of the available cores at a maximum.
  The parallelization uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) and
  works only in Unix systems

- verbose:

  logical, whether to print text to console

## Value

an STlist containing spatial statistics

## Details

The function computes global spatial autocorrelation statistics (Moran's
I and/or Geary's C) for the requested genes and samples. Then
computation uses the package `spdep`. The calculated statistics are
stored in the STlist, which can be accessed with the `get_gene_meta`
function. For visual comparative analysis, the function `compare_SThet`
can be used afterwards.

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
