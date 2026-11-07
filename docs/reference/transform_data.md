# transform_data: Transformation of spatial transcriptomics data

Applies data transformation methods to spatial transcriptomics samples
within an STlist

## Usage

``` r
transform_data(
  x = NULL,
  method = "log",
  scale_f = 10000,
  sct_n_regr_genes = 3000,
  sct_min_cells = 5,
  cores = NULL
)
```

## Arguments

- x:

  an STlist with raw count matrices.

- method:

  one of `log` or `sct`. If `log`, log-normalization is performed. If
  `sct`, then the SCTransform method is applied by calling
  [`sctransform::vst`](https://rdrr.io/pkg/sctransform/man/vst.html)

- scale_f:

  the scale factor used in logarithmic transformation

- sct_n_regr_genes:

  the number of genes to be used in the regression model during
  SCTransform. The function
  [`sctransform::vst`](https://rdrr.io/pkg/sctransform/man/vst.html)
  makes a random gene selection based on this number

- sct_min_cells:

  The minimum number of spots/cells to be used in the regression model
  fit by
  [`sctransform::vst`](https://rdrr.io/pkg/sctransform/man/vst.html)

- cores:

  integer indicating the number of cores to use during parallelization.
  If NULL, the function uses half of the available cores at a maximum.
  The parallelization uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) and
  works only in Unix systems.

## Value

x an updated STlist with transformed counts.

## Details

This function takes an STlist with raw counts and performs data
transformation. The user has the option to select between log
transformation after library size normalization (`method='log'`), or
SCTransform (`method='sct'`). In the case of logarithmic transformation,
a scaling factor (10^4 by default) is applied. The function uses
parallelization using "forking" (not available in Windows OS). Note that
the method `sct` returns a matrix with less genes as filtering is done
for low expression genes.

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
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```
