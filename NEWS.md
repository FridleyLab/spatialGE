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
  produced differently-ordered `STlist`s. The sort now applies, ordering by
  `xpos` then `ypos`, so ordering depends on the coordinates themselves rather
  than on input file layout.

  No function signatures, arguments, or return types changed, and no expression
  values changed -- every count, transformed count, and gene-level statistic is
  identical once reordered. What changes is order: the column order of
  `@counts`/`@tr_counts`, the row order of `@spatial_meta`, and therefore the row
  order of results built from them. Code that addresses spots by `libname` is
  unaffected; code that indexes them positionally, or that depends on how ties
  break in an `adj_p_value`-sorted result, may need updating.

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
* `.github/` is no longer shipped inside the package tarball (added to
  `.Rbuildignore`).

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
