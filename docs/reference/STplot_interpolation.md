# STplot_interpolation: Visualize gene expression surfaces

Produces a gene expression surface from kriging interpolation of ST
data.

## Usage

``` r
STplot_interpolation(
  x = NULL,
  genes = NULL,
  top_n = 10,
  samples = NULL,
  color_pal = "BuRd"
)
```

## Arguments

- x:

  an STlist containing results from `gene_krige` for the genes selected.

- genes:

  a vector of gene names (one or several) to plot. If 'top', the 10
  genes with highest standard deviation from each spatial sample are
  plotted.

- top_n:

  an integer indicating how many top genes to perform kriging. Default
  is 10.

- samples:

  a vector indicating the spatial samples to plot. If vector of numbers,
  it follows the order of `names(x@counts)`. If NULL, the function plots
  all samples

- color_pal:

  a color scheme from `khroma` or `RColorBrewer`.

## Value

a list of plots

## Details

This function produces a gene expression surface plot via kriging for
one or several genes and spatial samples

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
  melanoma <- gene_interpolation(melanoma, genes=c('MLANA', 'COL1A1'), samples='ST_mel1_rep2')
  kp = STplot_interpolation(melanoma, genes=c('MLANA', 'COL1A1'), samples='ST_mel1_rep2')
  ggpubr::ggarrange(plotlist=kp)
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```
