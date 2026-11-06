# STdiff_volcano: Generates volcano plots from STdiff results

Generates volcano plots of differential expression results from STdiff

## Usage

``` r
STdiff_volcano(
  x = NULL,
  samples = NULL,
  clusters = NULL,
  pval_thr = 0.05,
  color_pal = NULL
)
```

## Arguments

- x:

  the output of `STdiff`

- samples:

  samples to create plots

- clusters:

  names of the clusters to generate comparisons

- pval_thr:

  the p-value threshold to color genes with differential expression

- color_pal:

  the palette to color genes by significance

## Value

a list of ggplot objects

## Details

The function generated volcano plots (p-value vs. log-fold change) for
genes tested with `STdiff`. Colors can be customized to show
significance from spatial and non-spatial models
