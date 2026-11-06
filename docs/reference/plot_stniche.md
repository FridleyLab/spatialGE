# Plot STniche results

Fifteen plot types: clusters, Zi, ranking, composition, noise, violin,
zi_heatmap, celltype_heatmap, celltype_composition, celltype_stackedbar,
cluster_annotation, cluster_matching, celltype_spatial, niche_signature,
multi_sample_alluvial.

## Usage

``` r
plot_stniche(
  x = NULL,
  samples = NULL,
  type = c("clusters", "Zi", "ranking", "composition", "noise", "violin", "zi_heatmap",
    "celltype_heatmap", "celltype_composition", "cluster_annotation", "cluster_matching",
    "celltype_spatial", "celltype_stackedbar", "niche_signature",
    "multi_sample_alluvial"),
  zi = NULL,
  gene_set = NULL,
  palette = NULL,
  ptsize = 0.9,
  stroke = 0,
  reverse_y = FALSE,
  overlay_samples = FALSE,
  title = NULL,
  legend_name = NULL,
  pearson_thresh = 0.6,
  jaccard_thresh = 0.15,
  top_n_spots = 200L,
  pie_top_k = 3L,
  out_dir = NULL,
  highlight_types = NULL,
  highlight_cols = NULL,
  pool_samples = FALSE,
  bar_scale = c("percent", "count"),
  top_n_celltypes = NULL,
  min_cells = 1L,
  label_min_pct = 5
)
```

## Arguments

- x:

  An STlist object.

- samples:

  Character vector of sample names (or indices) to plot. Default = all
  samples. For `type = "multi_sample_alluvial"` a named list of sample
  groups is also accepted.

- type:

  Character. One of the supported plot types (see details).

- zi:

  Matrix of Local Moran Z values (rows = gene sets, columns =
  "sample&&libname"). Required for several plot types.

- gene_set:

  Character. Gene set / pathway name required for `type` in
  `c("Zi", "ranking", "violin")`.

- palette:

  Colour palette. Can be a named character vector, a function that
  returns colours, or a length-3 vector/named vector for Zi gradients.

- ptsize:

  Numeric. Point size for spatial plots.

- stroke:

  Numeric. Point stroke width (0 = no stroke).

- reverse_y:

  Logical. Reverse the y-axis.

- overlay_samples:

  Logical. Overlay multiple samples in one plot when possible.

- title:

  Character. Plot title (overrides the default title).

- legend_name:

  Character. Legend title.

- pearson_thresh:

  Numeric. Pearson correlation threshold (currently unused; reserved for
  future matching logic).

- jaccard_thresh:

  Numeric. Jaccard index threshold (currently unused; reserved for
  future matching logic).

- top_n_spots:

  Integer. Number of top spots (currently unused; reserved for future
  ranking logic).

- pie_top_k:

  Integer. Number of top gene sets shown in pie charts.

- out_dir:

  Character. Output directory (currently unused; reserved for optional
  file export).

- highlight_types:

  Cell types to highlight (only for `type = "celltype_spatial"`). Either
  an unnamed character vector, or a *named* vector whose names are
  display labels and whose values are regular expressions, e.g.
  `c(CD8 = "CD8", M2 = "M2|CD163")`. Patterns are matched with
  `grepl(..., ignore.case = TRUE)` against `cell_type`.

- highlight_cols:

  Named character of colours per highlighted type (names must equal the
  display labels).

- pool_samples:

  Logical; for `type = "celltype_stackedbar"`, pool all selected samples
  into one set of bars (default `FALSE`: one panel per sample).

- bar_scale:

  `"percent"` (default; each bar sums to 100) or `"count"` (absolute
  number of cells per niche) for `type = "celltype_stackedbar"`.

- top_n_celltypes:

  Optional integer; keep the `top_n_celltypes` most abundant cell types
  and collapse the rest into "Other" (`type = "celltype_stackedbar"`).

- min_cells:

  Minimum number of cells a niche needs to be drawn
  (`type = "celltype_stackedbar"`).

- label_min_pct:

  Segments at least this large (in percent) are labelled with their
  value (`type = "celltype_stackedbar"`).

## Value

A `ggplot` object, a named list of them (one per sample), a
`ComplexHeatmap` object/list, or a data frame, depending on `type`.
`"celltype_stackedbar"` and `"niche_signature"` return a single
`ggplot`.

## Details

**Why a stacked bar of cell types per niche?** A niche is only
biologically interpretable if its cellular make-up is known. The
`celltype_heatmap` shows the same information as a grid, which is
efficient for many cell types but makes it hard to judge how much of a
niche is occupied by each type. `type = "celltype_stackedbar"` draws one
bar per niche whose segments sum to 100\\ (i) the dominant cell type(s)
of every niche, and therefore its purity, can be read directly; (ii)
niches can be compared side by side within a sample; and (iii) the same
niche can be compared across samples (e.g. IO-naive vs IO-exposed) when
several samples are supplied, which exposes sample-specific composition
shifts that pooled summaries hide. Because the bars are compositional, a
high share of one type is not evidence of enrichment relative to the
tissue background; interpret them together with `celltype_heatmap` and
with the niche's absolute size (printed as `n` under each bar).
