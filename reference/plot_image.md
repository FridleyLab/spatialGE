# plot_image: Generate a ggplot object of the tissue image

Creates ggplot objects of the tissue images when available within the
STlist

## Usage

``` r
plot_image(x = NULL, samples = NULL)
```

## Arguments

- x:

  an STlist

- samples:

  a vector of numbers indicating the ST samples to plot, or their sample
  names. If vector of numbers, it follow the order of `names(x@counts)`.
  If NULL, the function plots all samples

## Value

a list of plots

## Details

If the STlist contains tissue images in the `@misc` slot, the
`plot_image` function can be used to generate ggplot objects. These
ggplot objects can be plotted next to quilt plots (`STplot` function)
for comparative analysis.
