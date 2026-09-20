# spatialGE 2.0.0

## Behavior changes

* **Requesting a gene, annotation or reference domain that does not exist is now
  an error rather than a silent no-op.** `SThet(genes=)`, `STgradient(annot=)`,
  `STgradient(ref=)` and `gene_interpolation(genes=)` previously accepted a name
  matching nothing and returned without complaint, so a typo looked like a run
  that legitimately found nothing. The rule is now: error when *nothing* matches,
  and message naming what was dropped when only some does. Code that relied on the
  no-op will now stop; that reliance was almost certainly accidental.

* **`STenrich()` now requires `annot` and `domain` together.** Domain-restricted
  mode was guarded by `!is.null(annot) & !is.null(domain)`, so supplying one
  without the other skipped the entire block: the **unrestricted** analysis ran and
  returned results a caller would reasonably read as domain-restricted, with
  nothing downstream able to detect the discrepancy. Supplying one alone is now an
  error. If you have results from 1.2.x produced by passing `annot` without
  `domain`, they are not domain-restricted and are worth regenerating.

* **`SThet()`'s default (`k=NULL`) statistics moved by ~2e-15.** The
  distance-weighted complete graph is no longer materialized as a `listw`; Moran's
  I and Geary's C are computed directly from the coordinates in C++. `Szero` is
  bit-identical, and the residual difference is summation order inside the
  weighted dot product -- four orders of magnitude below the ~2e-2 spot-order
  sensitivity already documented for `SThet(k=)`. The reason for the change is
  capacity: the old `listw` needed 1.2 GB for a single 4,096-spot sample and about
  181 GB extrapolated to 50,000 cells, i.e. it could not run on Xenium or CosMx
  data at all. It now completes in ~39 s using ~100 KB. `@misc$sthet$listws` is no
  longer populated on this path; `SThet(k=)` still builds a real `spdep` object.

* **`gene_stdevs` in `@gene_meta` moved by ~1e-16.** Per-gene standard deviations
  are computed without densifying the matrix, which changes the summation order.
  Nothing downstream consumes `gene_stdevs` -- `topgenes` selection uses
  `vst.variance.standardized` -- and cluster assignments, all `vst.*` columns and
  `@tr_counts` were verified bit-identical.

* **A Seurat object's platform is no longer assumed to be Visium.** `STlist()`
  hardcoded `@misc$platform = 'visium'` for any Seurat input; it is now inferred
  from the image class. This matters because `platform` selects `STgradient()`'s
  default `nb_dist_thr` and branches in `plot_image()` and
  `STplot_interpolation()`, so importing an imaging-based (Xenium/CosMx) Seurat
  object previously produced results computed with Visium-appropriate distance
  thresholds.

* `pseudobulk_heatmap()` displays a different set of genes, because it was
  selecting the wrong ones. The heatmap keeps the `hm_display_genes` most
  variable genes, ranked by coefficient of variation, but the CV was computed
  upside down as `mean/sd` instead of `sd/mean`. Sorting that decreasing put the
  *highest* `mean/sd` first -- that is, the **lowest** true CV -- so the plot
  showed the most uniformly-expressed genes rather than the most variable ones,
  the exact opposite of the intent. Any heatmap produced by 1.2.x is showing the
  wrong gene set and is worth regenerating.

* `cores=` is now honored by `STlist()` and `transform_data()`. Both accepted the
  argument and then ignored it, so they always ran at the auto-detected core
  count no matter what was passed. Passing `cores=` now actually changes the
  number of forked workers. Two consequences: wall-clock times for these two
  functions will differ from 1.2.x for any explicit `cores=`, and because forked
  children re-seed their own RNG, `transform_data(method='sct')` is only
  reproducible at `cores=1` once a sample exceeds `sctransform::vst`'s internal
  sampling threshold. Pass `cores=1` where exact SCT reproducibility matters.

* Spot/cell ordering within an `STlist` is now deterministic. `process_lists()`
  contained a coordinate sort that silently had no effect, so through 1.2.2
  spots simply retained whatever order they appeared in in the input files --
  meaning two files describing the same tissue with rows in a different order
  produced differently-ordered `STlist`s. The sort now applies, so ordering
  depends on the coordinates themselves rather than on input file layout.

  Note what it sorts by: the **third column of the coordinate data**, whichever
  that is, not `xpos` specifically. With the usual `libname`, `ypos`, `xpos`
  column order that means `xpos`, and because the sort is stable the result reads
  as `xpos` then `ypos`. Supply the coordinate columns in a different order and
  the sort key changes with them.

  No function signatures, arguments, or return types changed. Counts, transformed
  counts, and the `vst.*` gene-level statistics are bit-identical once matched by
  name. What changes for them is order: the column order of `@counts`/`@tr_counts`,
  the row order of `@spatial_meta`, and therefore the row order of results built
  from them. Code that addresses spots by `libname` is unaffected; code that indexes
  them positionally, or that depends on how ties break in an `adj_p_value`-sorted
  result, may need updating.

  Four results are **not** merely reordered, because they depend on spot order by
  construction rather than by layout. All four behaved this way in 1.2.x too --
  there, the order came from your input file, so the same tissue supplied with rows
  in a different order gave different numbers. Ordering is now determined by the
  coordinates, which removes that file-layout dependence, but the underlying
  sensitivity remains and is worth knowing about:

  - `STenrich()` p-values. The permutation null is drawn over spot *positions*, so
    a fixed `seed` selects a different random subset when the spots are in a
    different order. The observed statistic is unaffected; only the null, and
    therefore `p_value`/`adj_p_value`, move.
  - `SThet(k = )` Moran's I. With `k` given, neighbours come from a k-nearest
    neighbour graph, and on a regular grid the k-th distance is tied between
    several spots; the tie is broken by position. On our test data Moran's I moves
    by up to 0.02 -- large enough to matter. `SThet()` without `k` uses a
    distance-based weight list and is unaffected (agreement to ~1e-15).
  - `STgradient(distsumm = 'avg')` Spearman correlations. Averaging each spot's
    distances to the reference domain produces exact ties on a regular grid, and a
    change in summation order of ~1e-16 is enough to flip one, which flips a rank.
    Correlations move by up to 0.004. `distsumm = 'min'` does not create the same
    ties.
  - `plot_counts()` density curves, which plot a seeded random `distrib_subset` of
    spots selected by position. Cosmetic, but the curve is not identical. The
    point jitter in `distribution_plots(plot_type = 'violin')` behaves the same way.

  If you need a value that does not depend on spot order at all, prefer
  `SThet()` without `k` and `STgradient(distsumm = 'min')`, and raise `reps` in
  `STenrich()` so the permutation estimate is less sensitive to which subset is
  drawn.

## New features

* **Parallelism on Windows.** `resolve_cores()` returned `1L` unconditionally on
  Windows and `sge_lapply()` forwarded only to `parallel::mclapply`, which has no
  `fork()` there -- so the package was entirely serial. A non-forking (PSOCK)
  backend is now used on Windows, and can be forced anywhere with
  `options(spatialGE.psock=TRUE)`. `mclapply` remains the default where `fork()`
  exists. The `spatialGE.cores`/`SPATIALGE_CORES` override is also now consulted
  before any platform decision; it was previously ignored on Windows.

* **Seurat v5 support, and conversion in both directions.** Seurat input was v4
  only: v5 stores counts in an `Assay5`'s layers, so the old access returned
  `NULL`, and v5's `VisiumV2`/`FOV` images carry no `@coordinates` frame. Both are
  handled now, the assay defaults to `DefaultAssay()` rather than a hardcoded
  `'Spatial'`, and the import is sparse-to-sparse instead of densifying every gene
  by every cell across all slices up front. `as_Seurat()` is new and provides the
  reverse conversion.

* **Optional on-disk storage.** `STlist_to_disk(x, root)` moves count matrices
  into per-sample parquet files so data larger than memory can be held and read in
  gene-subset slices; `STlist_to_memory()` reverses it, exactly. Opt-in and off by
  default -- an STlist built the normal way is unchanged. `is_backed()` and
  `validate_store()` support it. Note this buys capacity, not speed: aggregates
  still compute in memory, so a backed STlist is currently slower on those paths.

* **Large-data performance.** Several operations no longer densify sparse
  matrices. Measured on a 200-gene x 50,000-cell synthetic sample: `filter_data()`
  27x faster, per-gene standard deviations 42x, and log transformation 33x
  (bit-identical output).

* **`STgradient()`'s distance step no longer materializes the n x n distance
  matrix.** `stgradient_dists_one_sample()` built `as.matrix(stats::dist(coords))`
  -- 934 MB allocated at 4,096 spots and about 20 GB at 50,000 cells, so the
  function could not run on Xenium or CosMx data at all. The distances are now
  computed in C++ in O(n) memory. Output is bit-identical for both
  `distsumm='min'` and `distsumm='avg'` (the latter reproduces R's two-pass
  corrected `mean()`), so no results change; measured 1.6x faster and about
  3,100x less allocation at 4,096 spots, and 50,000 cells now completes instead
  of failing to allocate.

* **`STdiff()`'s `t_test` and `wilcoxon` paths are 17-27x faster.**
  `stdiff_mean_test()` re-filtered and re-selected the whole expression frame
  once per gene x cluster combination; the two groups depend only on the
  annotation, so the split is now computed with base indexing and the gene column
  is pulled directly. Output is unchanged -- the same `t.test`/`wilcox.test`
  calls on the same vectors. The default `mm` path is not affected.


* Set the environment variable `SPATIALGE_CORES`, or the option
  `spatialGE.cores`, to force every parallelized call site to run serially and
  in-process without editing any code. Intended for profiling (`profvis` cannot
  see inside `parallel::mclapply`'s forks) and for debugging a worker that fails
  only under forking. The option takes precedence over the environment variable.
  This is a development aid and is deliberately not part of the documented API.

## Bug fixes

* `pseudobulk_samples()` no longer fails on a gene whose pseudobulk total is
  constant across samples. Such a gene produced an all-`NaN` column after scaling
  and `prcomp()` rejected it outright, so any dataset containing an all-zero gene
  -- routine after a permissive `filter_data()` -- could not be analysed. Those
  genes are now dropped with a warning.

* `SThet()` no longer fails when `@spatial_meta` has been row-subset with base
  `[`. Both spatial-weight builders fed the frame to
  `tibble::column_to_rownames()`, which rejects non-contiguous row names, so the
  call died with an error naming neither `SThet` nor the cause. This is the same
  defect that made `STenrich()`'s `domain=` path unreachable through 1.2.2, in two
  call sites that fix did not cover.

* `STenrich()` no longer crashes on a gene set with exactly one gene present.
  Subsetting dropped the matrix to a dimensionless vector, which
  `DelayedMatrixStats::colMeans2()` rejects -- even though `min_genes=1` explicitly
  permits single-gene sets.

* `Seurat_FindVariableFeatures()` no longer fails on any input that is not already
  a `dgCMatrix`. It called `as.sparse()`, which is `SeuratObject::as.sparse()` and
  was never imported, so the call had been unresolvable since the code was
  vendored -- masked because every caller happened to pass a `dgCMatrix`, and
  because the name was listed in `utils::globalVariables()`, suppressing the check
  NOTE that would have caught it.

* `STlist()` accepts a Seurat object with `samples=NULL`. The sample-name parsing
  tested an argument derived from `samples` before reaching the Seurat branch, so
  the normal way to call it raised "argument is of length zero".

* `filter_data()` can now produce a sample with exactly one surviving spot. The
  dense intermediate collapsed to a vector at one column, which then broke
  row-name handling.

* `transform_data()`, `SThet()`, `STenrich()` and `gene_interpolation()` report
  invalid input as spatialGE errors naming the function and argument, instead of
  surfacing internal failures such as `object 'tr_results' not found` or
  `attempt to select less than one element in get1index`.


* `STenrich(annot=, domain=)` now works. Domain-restricted gene set testing was
  announced in 1.2.1 but never functioned: restricting to a domain subsets the
  spatial metadata rows, which leaves stray row names that
  `tibble::column_to_rownames()` rejects, so the call always failed with
  ``​`.data` must be a data frame without row names.`` Every domain-restricted
  call in 1.2.1 and 1.2.2 errored out. The unrestricted and `score_type='gsva'`
  paths were unaffected.
* All 14 documented function examples now run. A doubled roxygen comment marker
  (`#'   #' zip_tmp = ...`) commented out the assignment that the following
  `unzip()` call depended on, so every example failed with
  `object 'zip_tmp' not found`. Each example's `tryCatch` reported this as
  "Could not run example. Are you connected to the internet?", misattributing a
  code defect to the network, while `R CMD check` still reported
  `checking examples ... OK`.
* `STgradient()` no longer leaks state between genes. `pval_warn` was only assigned
  inside a branch that requires more than one spot in the fitted range, so a gene
  that fell below that threshold left the previous gene's value in place. Affects
  the `pval_comment` column only.
* `.github/` is no longer shipped inside the package tarball (added to
  `.Rbuildignore`).

## Minor

* The error raised when an argument that should be an `STlist` is `NULL` or of the
  wrong type now reads `cannot get a slot (...) from an object of type ...` rather
  than ``no applicable method for `@` applied to an object of class ...``, because
  the shared argument-resolution helper reads slots via `methods::slot()`. Both
  messages are raw R errors and neither is informative; nothing is caught or
  handled differently.

* `STdiff()` no longer emits dplyr's "Detected an unexpected many-to-many
  relationship" warning under dplyr 1.1 and later. The gene x annotation table it
  builds is deliberately many-to-many -- every gene of a sample combines with every
  annotation of that sample -- but it was declared with the deprecated
  `multiple='all'`, which newer dplyr checks separately from `relationship=` and
  therefore ignores. The join now declares `relationship='many-to-many'`. The
  joined rows, and every result built from them, are unchanged.

# spatialGE 1.2.2

# spatialGE 1.2.1

* Added C++ code to speed up `STenrich` computations (Thank you, Dr. Soupir!).
* The `STenrich` function can now calculate gene set enrichment scores (via `GSVA`) in addition to gene set average expression.
* The `STenrich` function can now test for gene set "hot-spots" within an specific tissue domain (`annot` and `domain` arguments).
* Spatial distances can now be log-transformed for `STgradient` (`log` argument).
* Example data sets are now in the spatialGE_Data GitHub repository.
* New functions `spatial_metadata` and `tissue_names` to quickly access the names of spot/cell annotations and sample names.
* Re-assessed package dependencies.
* Use of `DelayedArray` for some calculations.
* Several bug fixes.

# spatialGE 1.2.0

* Functions `STenrich`, `STgradient`, and `STdiff` available.
* Added C++ code for Seurat's implementation of `FindVariableFeatures`

# spatialGE 1.1.0

* Multiple, significant changes to STlist (incompatible with STlist objects from version 1.0).
