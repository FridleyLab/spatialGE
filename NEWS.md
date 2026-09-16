# spatialGE 2.0.0

## Behavior changes

* Spot/cell ordering within an `STlist` has changed. `process_lists()` sorted
  coordinates in a way that silently had no effect, so spots retained their
  input file order; the sort now actually applies, ordering by `xpos` then
  `ypos`. No function signatures, arguments, or return types changed, and no
  expression values changed -- every count, transformed count, and gene-level
  statistic is identical once reordered. What does change is order: the column
  order of `@counts`/`@tr_counts`, the row order of `@spatial_meta`, and
  therefore the row order of results built from them. Code that addresses
  spots by `libname` is unaffected; code that indexes them positionally, or
  that depends on how ties break in a `adj_p_value`-sorted result, may need
  updating.

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
