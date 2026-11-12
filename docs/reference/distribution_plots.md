# per_unit_counts: Generates distribution plots of spot/cell meta data or gene expression

Generates violin plots, boxplots, or density plots of variables in the
spatial meta data or of gene expression

## Usage

``` r
distribution_plots(
  x = NULL,
  plot_meta = NULL,
  genes = NULL,
  samples = NULL,
  data_type = "tr",
  color_pal = "okabeito",
  plot_type = "violin",
  ptsize = 0.5,
  ptalpha = 0.5
)
```

## Arguments

- x:

  an STlist

- plot_meta:

  vector of variables in `x@spatial_meta` to plot distributions. If
  'total_counts', the function plots the counts per spot/cell. If
  'total_genes', the function plots the number of genes per spot/cell
  are plotted

- genes:

  vector of genes to plot expression distribution. If used in
  conjunction with `plot_meta`, the expression values are grouped using
  that variable

- samples:

  samples to include in the plot. Default (NULL) includes all samples

- data_type:

  one of 'tr' or 'raw', to plot transformed or raw counts

- color_pal:

  a string of a color palette from `khroma` or `RColorBrewer`, or a
  vector with colors

- plot_type:

  one of "violin", "box", or "density" (violin plots, box plots, or
  density plots respectively). If `plot_meta` and `gene` are used
  together, then density plots are disabled

- ptsize:

  the size of points in the plots

- ptalpha:

  the transparency of points (violin/box plot) or curves (density plots)

## Value

a list containing `ggplot2` objects

## Details

The function allows to visualize the distribution of spot/cell total
counts, total genes, or expression of specific genes across all samples
for comparative purposes. It also allows grouping of gene expression
values by categorical variables (e.g., clusters).
