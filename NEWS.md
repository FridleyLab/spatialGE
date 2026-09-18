# spatialGE 2.0.0

## Behavior changes

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

* Set the environment variable `SPATIALGE_CORES`, or the option
  `spatialGE.cores`, to force every parallelized call site to run serially and
  in-process without editing any code. Intended for profiling (`profvis` cannot
  see inside `parallel::mclapply`'s forks) and for debugging a worker that fails
  only under forking. The option takes precedence over the environment variable.
  This is a development aid and is deliberately not part of the documented API.

## Bug fixes

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
