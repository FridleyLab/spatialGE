# as_Seurat: Convert an STlist to a Seurat object

Converts an STlist into a Seurat (v5) object, one image/slice per
sample, so spatialGE results can be taken into Seurat-based workflows.

## Usage

``` r
as_Seurat(x = NULL, assay = "Spatial", samples = NULL)
```

## Arguments

- x:

  an STlist

- assay:

  the assay name to create. Default `'Spatial'`

- samples:

  the samples to include. If NULL (default), all samples

## Value

a Seurat object

## Details

Each sample in the STlist becomes one FOV in the returned object, named
after the sample. Raw counts go to the `counts` layer and, when
[`transform_data()`](https://fridleylab.github.io/spatialGE/reference/transform_data.md)
has been run, transformed counts go to the `data` layer. Spot-level
columns from `@spatial_meta` – including `STclust` cluster assignments –
are carried into `@meta.data`, prefixed by sample where a merge would
otherwise collide.

Cell names are prefixed with the sample name (`<sample>_<libname>`)
because spot barcodes are only unique within a sample in the general
case, and a merged Seurat object requires globally unique cell names.
