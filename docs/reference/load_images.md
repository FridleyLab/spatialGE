# load_images: Place tissue images within STlist

Loads the images from tissues to the appropriate STlist slot.

## Usage

``` r
load_images(x = NULL, images = NULL)
```

## Arguments

- x:

  an STlist

- images:

  a string indicating a folder to load images from

## Value

an STlist with images

## Details

This function looks for `.PNG` or `.JPG` files within a folder matching
the sample names in an existing STlist. Then, loads the images to the
STlist which can be used for plotting along with other spatialGE plots.
