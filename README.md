# spatialGE

<img align="left" src="man/figures/logo.png" height="100" width="90" alt="spatialGE logo" />

An R package for visualizing and analyzing spatially-resolved transcriptomics data from spot-based
platforms (e.g., 10X Visium) and imaging-based platforms (e.g., Xenium, CosMx). Data and results
from multiple tissue sections are stored in a single **STlist** (Spatial Transcriptomics List) object.

## Features

- **Visualization:** `STplot`, `gene_interpolation`, `STplot_interpolation` for exploring gene expression in spatial context.
- **Spatial autocorrelation:** `SThet`, `compare_SThet` to quantify spatial uniformity (Moran's I, Geary's C) and relate it to sample metadata.
- **Tissue domain detection:** `STclust` for spatially-informed hierarchical clustering.
- **Gene set enrichment:** `STenrich` to detect gene sets with spatial patterns, optionally within a tissue domain.
- **Spatial gradients:** `STgradient` to find genes whose expression varies with distance to a tissue domain.
- **Differential expression:** `STdiff` for mixed models with spatial covariance, plus non-spatial Wilcoxon and t-tests.
- **Spatial niche detection:** `STniche` to identify tissue niches from the local spatial organization of gene sets, with `plot_stniche` for visualization (see below).

<p align="center">
<img src="man/figures/workflow.png" width="700" alt="spatialGE workflow">
</p>

## Spatial niche detection with STniche

`STniche` detects niches, local microenvironments where particular cell types and molecular programs
co-occur, in three steps:

1. **Gene set scoring:** each spot or cell receives an activity score for every gene set of interest.
2. **Spatial autocorrelation:** a local Moran's I statistic (reported as a Z-score, Zi) measures whether
   each gene set is spatially concentrated around that location.
3. **Clustering:** Zi profiles are clustered with model-based clustering, and the number of niches is
   chosen by BIC.

It works on spot-based (Visium) and single-cell (CosMx) data. When run on several samples at once, niche
labels (C1, C2, ...) are shared and directly comparable across samples.

```r
library(spatialGE)

# gene_sets: a named list of character vectors, one per gene set
res <- STniche(x = stlist, gene_sets = gene_sets,
               samples = names(stlist@spatial_meta),
               ks = 10, knns = 22, out_i = TRUE, cores = 4)

st <- res[[1]]  # STlist with niche assignments
zi <- res[[2]]  # Local Moran's Z (Zi) matrix
```

`plot_stniche` visualizes the results through its `type` argument:

| `type` | Shows |
|---|---|
| `"clusters"` | Niche map of each sample |
| `"Zi"` | Local Moran's Z of one gene set (requires `zi` and `gene_set`) |
| `"celltype_heatmap"` | Fraction of each niche's cells belonging to each cell type |
| `"celltype_composition"` | Stacked bars of cell-type composition per niche |
| `"multi_sample_alluvial"` | Niche correspondence between samples through shared gene-set programs (requires `zi`) |

```r
plot_stniche(st, type = "clusters", samples = "sample1")
plot_stniche(st, type = "Zi", zi = zi, gene_set = "TLS_core", samples = "sample1")
```

The STniche vignette applies these functions to two case studies: Visium sections from renal cell
carcinoma (GSE175540), where niches are compared against pathologist-annotated tertiary lymphoid
structures, and CosMx clear cell RCC fields of view, where niches are characterized by cell-type
composition and immune signaling in IO-naive and IO-exposed tumors.

## What's new in 2.0.0

**New**

- `STniche` and `plot_stniche` for spatial niche detection on Visium and CosMx data.
- Subset an `STlist` by sample with `x[i]` (name, position, negative position, or logical); `length(x)` now returns the number of samples. `x[i, j]` and `x[[i]]` raise an error by design.
- Seurat v5 support (uses `DefaultAssay()`, sparse import) and `as_Seurat()` for converting back.
- Optional on-disk (parquet) storage for data larger than memory: `STlist_to_disk()`, `STlist_to_memory()`, `is_backed()`, `validate_store()`. It adds capacity, not speed.
- Parallel processing on Windows (PSOCK backend; force it elsewhere with `options(spatialGE.psock = TRUE)`).
- `filter_data()` now reports what it removed (`verbose = TRUE`), or that nothing was filtered.
- Scales to Xenium- and CosMx-sized samples (200 genes x 50,000 cells benchmark):
  - `filter_data` 27x faster, per-gene SD 42x, log transformation 33x (identical output).
  - `STdiff` t-test/Wilcoxon 17-27x faster (default `mm` path unchanged).
  - `STgradient` distances computed in O(n) memory (previously ~20 GB at 50,000 cells); output is bit-identical.
  - `SThet` (default `k = NULL`) no longer builds a weights list (previously ~1.2 GB for 4,096 spots).

**Changes that may affect existing code**

- `dim()` returns `c(n_genes, n_spots)` for a single sample and `NULL` for multiple samples; use `summarize_STlist()` for per-sample dimensions.
- `STlist` objects are validated at construction.
- Unmatched genes, annotations, or domains now raise an error instead of being silently ignored (`SThet(genes=)`, `STgradient(annot=, ref=)`, `gene_interpolation(genes=)`). Partial matches give a message naming what was dropped.
- `STenrich()` requires `annot` and `domain` together. In 1.2.x, passing only one silently ran the *unrestricted* analysis; regenerate those results.
- `transform_data()` now populates `vst.variance.standardized` in `@gene_meta` (for `log` and `sct`), which now precedes `moran_i` and `geary_c`; select columns by name. Samples with fewer than 4 non-constant genes get a warning and no VST column, instead of failing.
- A Seurat object's platform is now inferred from its image class (it was always assumed Visium), which affects `STgradient()`'s default `nb_dist_thr` for Xenium/CosMx imports.
- `cores=` is now honored by `STlist()` and `transform_data()`. Use `cores = 1` when exact `sct` reproducibility matters on large samples.
- Spots are now ordered deterministically by coordinates (the third coordinate column, normally `xpos` then `ypos`). Values are identical when matched by name, but column order of `@counts`/`@tr_counts` and row order of `@spatial_meta` can change; address spots by `libname`, not position.
- Some results depend on spot order by construction (as in 1.2.x, but now fixed by coordinates rather than file layout): `STenrich()` p-values, `SThet(k=)` Moran's I (up to ~0.02), `STgradient(distsumm = 'avg')` correlations (up to ~0.004), and density/jitter in `plot_counts()` and violin plots. For order-independent values, use `SThet()` without `k`, `distsumm = 'min'`, and a higher `reps` in `STenrich()`.
- Negligible numeric shifts: `SThet()` default (~2e-15) and `gene_stdevs` (~1e-16). `@misc$sthet$listws` is no longer populated on the default path.
- Validation errors now carry the `spatialGE exception:` prefix; match on that rather than on exact wording.

**Fixes**

- `pseudobulk_heatmap()` previously displayed the *least* variable genes (CV was computed as mean/sd). Heatmaps from 1.2.x should be regenerated.
- `STenrich(annot=, domain=)` now works; it failed in 1.2.1 and 1.2.2. A single-gene set no longer crashes it.
- Objects returned by `gene_interpolation()` no longer fail `validObject()`.
- `load_images()` matches images by normalized sample name across all directories and warns on unreadable files, instead of dropping samples.
- `STplot_interpolation()` works with its default `genes = NULL`.
- `STdiff()` gives a clear error for numeric or single-group `annot=`, and no longer emits dplyr's many-to-many warning.
- `transform_data(method = 'sct')` errors now name the sample and its dimensions; invalid `cores=` is caught on every path.
- `pseudobulk_samples()` drops constant genes with a warning instead of failing.
- `SThet()` works on a row-subset `@spatial_meta`; `STlist()` accepts a Seurat object with `samples = NULL`; `filter_data()` can keep exactly one spot; `Seurat_FindVariableFeatures()` accepts non-`dgCMatrix` input; `STgradient()` no longer leaks `pval_comment` between genes.
- All 14 documented examples now run.

**Earlier releases:** 1.2.1 added C++ speedups and GSVA scoring for `STenrich`, domain-restricted enrichment, log distances in `STgradient`, and `spatial_metadata()`/`tissue_names()`. 1.2.0 introduced `STenrich`, `STgradient`, and `STdiff`. 1.1.0 changed `STlist` in ways incompatible with 1.0 objects.

See [NEWS.md](NEWS.md) for the full list of changes.

## Installation

```r
options(timeout = 9999999) # Avoid connection timeouts when downloading from GitHub
devtools::install_github("fridleylab/spatialGE")
```

## Documentation

- Tutorials: https://fridleylab.github.io/spatialGE/
- Source code: https://github.com/FridleyLab/spatialGE
- **[spatialGE-Web](https://spatialge.moffitt.org):** a point-and-click web application for using spatialGE without coding. It supports Visium outputs and csv/tsv expression files paired with coordinate files.

## Citation

If you use spatialGE, please cite:

Ospina, O. E., Wilson, C. M., Soupir, A. C., Berglund, A., Smalley, I., Tsai, K. Y., Fridley, B. L. (2022).
spatialGE: quantification and visualization of the tumor microenvironment heterogeneity using spatial
transcriptomics. *Bioinformatics*, 38, 2645-2647. https://doi.org/10.1093/bioinformatics/btac145

Ospina, O. E., Manjarres-Betancur, R., Gonzalez-Calderon, G., Soupir, A. C., Smalley, I., Tsai, K. Y.,
Markowitz, J., Khaled, M. L., Vallebuona, E., Berglund, A. E., Eschrich, S. A., Yu, X., Fridley, B. L. (2025).
spatialGE is a user-friendly web application that facilitates spatial transcriptomics data analysis.
*Cancer Research*, 85(5), 848-858. https://doi.org/10.1158/0008-5472.CAN-24-2346

When using `STdiff`, please also cite:

Ospina, O. E., Soupir, A. C., Manjarres-Betancur, R., Gonzalez-Calderon, G., Yu, X., Fridley, B. L. (2024).
Differential gene expression analysis of spatial transcriptomic experiments using spatial mixed models.
*Scientific Reports*, 14, 10967. https://doi.org/10.1038/s41598-024-61758-0

<!-- Visitor-counter badge, commented out because R CMD check --as-cran flags its URL
     as possibly invalid. Uncomment to restore.
<p align="center">
<a href="https://clustrmaps.com/site/1bwt8">
<img src="https://clustrmaps.com/map_v2.png?cl=a0d9e6&w=a&t=tt&d=MaZTayTGh9YYcDsI3LkJtKTJFTPygJFYO9LnIVaLr2Q&co=062f6e" />
</a>
</p>
-->