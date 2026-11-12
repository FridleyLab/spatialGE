# spatialGE

![spatialGE logo](reference/figures/logo.png)

An R package for visualizing and analyzing spatially-resolved
transcriptomics data from spot-based platforms (e.g., 10X Visium) and
imaging-based platforms (e.g., Xenium, CosMx). Data and results from
multiple tissue sections are stored in a single **STlist** (Spatial
Transcriptomics List) object.

## Features

- **Visualization:** `STplot`, `gene_interpolation`,
  `STplot_interpolation` for exploring gene expression in spatial
  context.
- **Spatial autocorrelation:** `SThet`, `compare_SThet` to quantify
  spatial uniformity (Moran’s I, Geary’s C) and relate it to sample
  metadata.
- **Tissue domain detection:** `STclust` for spatially-informed
  hierarchical clustering.
- **Gene set enrichment:** `STenrich` to detect gene sets with spatial
  patterns, optionally within a tissue domain.
- **Spatial gradients:** `STgradient` to find genes whose expression
  varies with distance to a tissue domain.
- **Differential expression:** `STdiff` for mixed models with spatial
  covariance, plus non-spatial Wilcoxon and t-tests.
- **Spatial niche detection:** `STniche` to identify tissue niches from
  the local spatial organization of gene sets, with `plot_stniche` for
  visualization (see below).

![spatialGE workflow](reference/figures/workflow.png)

## Spatial niche detection with STniche

`STniche` detects niches, local microenvironments where particular cell
types and molecular programs co-occur, from spot-based (Visium) and
single-cell (CosMx) data. Niche labels are shared across samples
analyzed together, so they are directly comparable. `plot_stniche`
visualizes the results.

``` r

res <- STniche(x = stlist, gene_sets = gene_sets,
               samples = names(stlist@spatial_meta),
               ks = 2:10, knns = 22, out_i = TRUE)

plot_stniche(res[[1]], type = "clusters")
```

The method, all `plot_stniche` plot types, and worked Visium (GSE175540)
and CosMx (ccRCC) examples are in the STniche vignette:
[`vignette("vignette_stniche", package = "spatialGE")`](https://fridleylab.github.io/spatialGE/articles/vignette_stniche.md),
also available under Tutorials below.

## Supported input

| Platform | How to load |
|----|----|
| 10X Visium | [`STlist()`](https://fridleylab.github.io/spatialGE/reference/STlist.md) on one Space Ranger folder per sample |
| Xenium | Load with [`Seurat::LoadXenium()`](https://satijalab.org/seurat/reference/ReadXenium.html), then [`STlist()`](https://fridleylab.github.io/spatialGE/reference/STlist.md) on the Seurat object |
| CosMx | [`STlist()`](https://fridleylab.github.io/spatialGE/reference/STlist.md) on a Seurat object |
| Any other platform | [`STlist()`](https://fridleylab.github.io/spatialGE/reference/STlist.md) on a counts matrix plus a coordinates table (generic format) |

``` r

# Generic format: counts (genes x spots/cells) and coordinates (libname, ypos, xpos)
st <- STlist(rnacounts  = list(sample1 = counts),
             spotcoords = list(sample1 = coords),
             samples    = "sample1")
```

Step-by-step import instructions for each platform are in the *Basic
functionality* tutorial (section *Loading data from other platforms*):
[`vignette("basic_functions_vignette", package = "spatialGE")`](https://fridleylab.github.io/spatialGE/articles/basic_functions_vignette.md).

## What’s new in 2.0.0

- `STniche` and `plot_stniche` for spatial niche detection on Visium and
  CosMx data.
- Subset an `STlist` by sample with `x[i]`; `length(x)` returns the
  number of samples.
- Seurat v5 support, and
  [`as_Seurat()`](https://fridleylab.github.io/spatialGE/reference/as_Seurat.md)
  for converting an `STlist` back to Seurat.
- Optional on-disk (parquet) storage for data larger than memory
  ([`STlist_to_disk()`](https://fridleylab.github.io/spatialGE/reference/STlist_to_disk.md),
  [`STlist_to_memory()`](https://fridleylab.github.io/spatialGE/reference/STlist_to_memory.md)).
- Parallel processing on Windows.
- Large-data performance: `filter_data`, `STdiff`, `STgradient` and
  `SThet` now scale to Xenium- and CosMx-sized samples.
- [`filter_data()`](https://fridleylab.github.io/spatialGE/reference/filter_data.md)
  reports what it removed.

Version 2.0.0 also includes behavior changes that may affect existing
code (for example [`dim()`](https://rdrr.io/r/base/dim.html) on an
`STlist`, errors for unmatched genes or domains, and the order of spots)
and many bug fixes. These, together with earlier releases, are listed in
[NEWS.md](https://fridleylab.github.io/spatialGE/NEWS.md).

## Installation

``` r

options(timeout = 9999999) # Avoid connection timeouts when downloading from GitHub
devtools::install_github("fridleylab/spatialGE")
```

## Documentation

- Tutorials (basic functionality and data input, STdiff,
  STenrich/STgradient, STniche):
  <https://fridleylab.github.io/spatialGE/>
- Source code: <https://github.com/FridleyLab/spatialGE>
- **[spatialGE-Web](https://spatialge.moffitt.org):** a point-and-click
  web application for using spatialGE without coding. It supports Visium
  outputs and csv/tsv expression files paired with coordinate files.

## Citation

If you use spatialGE, please cite:

Ospina, O. E., Wilson, C. M., Soupir, A. C., Berglund, A., Smalley, I.,
Tsai, K. Y., Fridley, B. L. (2022). spatialGE: quantification and
visualization of the tumor microenvironment heterogeneity using spatial
transcriptomics. *Bioinformatics*, 38, 2645-2647.
<https://doi.org/10.1093/bioinformatics/btac145>

Ospina, O. E., Manjarres-Betancur, R., Gonzalez-Calderon, G., Soupir, A.
C., Smalley, I., Tsai, K. Y., Markowitz, J., Khaled, M. L., Vallebuona,
E., Berglund, A. E., Eschrich, S. A., Yu, X., Fridley, B. L. (2025).
spatialGE is a user-friendly web application that facilitates spatial
transcriptomics data analysis. *Cancer Research*, 85(5), 848-858.
<https://doi.org/10.1158/0008-5472.CAN-24-2346>

When using `STdiff`, please also cite:

Ospina, O. E., Soupir, A. C., Manjarres-Betancur, R., Gonzalez-Calderon,
G., Yu, X., Fridley, B. L. (2024). Differential gene expression analysis
of spatial transcriptomic experiments using spatial mixed models.
*Scientific Reports*, 14, 10967.
<https://doi.org/10.1038/s41598-024-61758-0>
