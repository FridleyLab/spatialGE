# pseudobulk_dim_plot: Plot PCA of pseudobulk samples

Generates a PCA plot after computation of "pseudobulk" counts

## Usage

``` r
pseudobulk_dim_plot(
  x = NULL,
  color_pal = "muted",
  plot_meta = NULL,
  dim = "pca",
  pcx = 1,
  pcy = 2,
  ptsize = 5
)
```

## Arguments

- x:

  an STlist with pseudobulk PCA results in the `@misc` slot (generated
  by `pseudobulk_samples`)

- color_pal:

  a string of a color palette from khroma or RColorBrewer, or a vector
  of color names or HEX values. Each color represents a category in the
  variable specified in `plot_meta`

- plot_meta:

  a string indicating the name of the variable in the sample metadata to
  color points in the PCA plot

- dim:

  one of `umap` or `pca`. The dimension reduction to plot

- pcx:

  integer indicating the principal component to plot in the x axis

- pcy:

  integer indicating the principal component to plot in the y axis

- ptsize:

  the size of the points in the PCA plot. Passed to the `size` aesthetic
  from `ggplot2`

## Value

a ggplot object

## Details

Generates a Principal Components Analysis plot to help in initial data
exploration of differences among samples. The points in the plot
represent "pseudobulk" samples. This function follows after usage of
`pseudobulk_samples`.

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
  melanoma <- pseudobulk_samples(melanoma)
  pseudobulk_dim_plot(melanoma, plot_meta='patient')
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```
