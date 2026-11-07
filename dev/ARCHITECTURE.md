# spatialGE architecture

Working notes for developing on `alex-refactor`. This documents the shape the
profiling/fidelity refactor left behind, not aspirational design -- if code and
doc disagree, trust the code and fix this file.

## STlist slot contract

Defined in `R/classDefinitions.R`. All list-valued slots are **sample-major**:
named by sample, in the same key order as `@counts` unless noted otherwise in
fidelity rule 7 below.

| Slot | Type | Contents |
|---|---|---|
| `counts` | `list` of matrix-like | raw per-spot/cell x gene counts, one entry per sample |
| `spatial_meta` | `list` of `data.frame` | per-spot/cell `libname`/`xpos`/`ypos` + any spot-level annotations. **A plain `data.frame`, not a tibble** -- see the row-names warning below |
| `gene_meta` | `list` of `tibble` | per-gene stats (mean/variance, VST, Moran's I / Geary's C once `SThet` has run) |
| `sample_meta` | `tibble` | one row per sample, user-supplied sample-level metadata |
| `tr_counts` | `list` of matrix-like | transformed counts (log or SCT), populated by `transform_data()` |
| `gene_krige` | `list` | kriging interpolation results (`gene_interpolation()`); **not** sample-major the way the others are -- see the Phase 7 accessor note below |
| `misc` | `list` | grab-bag of parameters/derived objects that don't fit the slots above -- see the key census next |

**`spatial_meta` row names are a live footgun.** Because it is a plain
`data.frame` and not a tibble, base row-subsetting (`d[keep, ]`) leaves
non-contiguous row names behind (`1, 5, 9, ...`), which
`tibble::column_to_rownames()` then rejects with ``​`.data` must be a data frame
without row names.`` That was a real, released bug in `STenrich()`'s
`domain=`-restricted path -- it made the entire domain-restricted branch
unreachable while the unrestricted branch worked fine, because only the
restricted branch subsets rows. `dplyr` verbs (`filter`, joins) reset row names
and are safe; `dplyr::select()` and base `[` do **not**. If you row-subset
`spatial_meta` and then hand it to `column_to_rownames()`, insert
`tibble::remove_rownames()` first.

Four call sites beyond `STenrich()` reach `column_to_rownames()` from
`spatial_meta`. Two are now defended outright: `create_listw_from_dist()` and
`create_listw_from_knn()` (`R/utils.R`) call `remove_rownames()` first, which is
output-neutral there because both overwrite the resulting dimnames from
`libname` anyway. Before that, **any** STlist whose `spatial_meta` had been
row-subset with base `[` made `SThet()` fail outright -- both weightings, with
an error naming neither `SThet` nor the real cause. Pinned by
`tests/testthat/test-numeric-sthet.R`.

`stgradient_dists_one_sample()` no longer reaches `column_to_rownames()`
at all: the item-1b rewrite takes row positions with `which()` and hands the
coordinate matrix to C++, so it is safe on a row-subset frame now.

The remaining one (`R/STclust.R`'s `calculate_dist_matrices`) is still only
safe because its caller happens to pass the whole, unsubset frame via
`sample_spatial_meta()`. That is an invariant nothing enforces, so treat it as
load-bearing when touching that caller.

Internal, non-exported accessors for the sample-major slots live in
`R/utils_accessors.R` (`sample_counts()`, `sample_spatial_meta()`,
`sample_tr_counts()`, `sample_gene_meta()`). They are a seam, not a validation
layer -- they don't check that `sample` exists, matching how `x@slot[[sample]]`
already behaved at every call site they replaced. Route new sample-major slot
reads through them; that's what would let `gene_krige` become sample-major
without touching every caller again.

## `misc` key census

`@misc` has no schema; these are the keys currently written and where:

| Key | Writer | Readers |
|---|---|---|
| `sp_images` | `load_images()` (`R/utils.R`) | `plot_image.R`, `plot_spatial_geneset.R` |
| `platform` | `STlist()` | `STgradient.R`, `plot_image.R`, `STplot_interpolation.R`, `STplot.R` |
| `transform` | `transform_data()` | `plot_spatial_expression.R` |
| `sthet` | `SThet()` (stores the per-sample `listw`/`nb` objects, among other things) | consumed within `SThet.R`/`compare_SThet.R` |
| `gene_krige` | `gene_interpolation()` | `STplot_interpolation.R` |
| `pbulk_pca`, `pbulk_pca_var`, `pbulk_umap`, `scaled_pbulk_mtx` | `pseudobulk_samples()` | `pseudobulk_samples.R` (self-consumed across its own sub-functions) |
| `image_scaling` | `plot_image.R` | `plot_image.R` |

Adding a new `misc` key is fine, but add a row here in the same commit.

`store` is written by `STlist_to_disk()` and read by `validate_store()` /
`STlist_to_memory()` (`R/store_parquet.R`). It holds `root`, `sample_order`,
`schema_ver` and `portable`.

## Optional on-disk backing (parquet COO)

**Opt-in, reversible, and off by default.** `STlist_to_disk(x, root)` replaces
each `@counts`/`@tr_counts` element with an `sge_manifest` describing a parquet
file; `STlist_to_memory(x)` reverses it. Every other slot is untouched. Nothing
changes for an STlist built the normal way.

It works because all matrix reads go through `R/utils_accessors.R`:
`sample_counts()`/`sample_tr_counts()` materialize a manifest on demand (and
forward `genes=` so a subset read touches only the row groups it needs), while
`sample_dim()`/`sample_genes()`/`sample_cells()` answer from the manifest with
zero I/O. That split is the whole reason the accessor layer was widened before
the store existed — `base::nrow()` returns `NULL` on a manifest rather than
erroring, so an unrouted dim read would have silently produced `NULL`.

Layout: one file per sample per slot, schema
`(gene_idx int32, cell_idx int32, value double)`, **sorted by
`(gene_idx, cell_idx)`**. Sorting by gene gives row-group min/max statistics so
`filter(gene_idx %in% ...)` prunes whole row groups — the dominant access pattern
here. Including `cell_idx` makes a per-gene aggregate visit values in the order
`Matrix::rowMeans` walks a row, which is what would let an arrow-side aggregate
agree bit-for-bit. Ordering comes free from `Matrix::t()` rather than an explicit
sort: a transposed `dgCMatrix` has the original genes as its columns.

**The manifest is the only authority for dimensions, dimnames and order.** COO
omits structural zeros, so an all-zero gene has *no rows at all* — and
`filter_data()`'s defaults are all `>= 0`, so such genes are retained and real.
Reconstructing dimnames from the data would silently shrink and reorder the
matrix, which fidelity rule 13(c) makes a genuine regression. Pinned by a test.

Constraints honored, each for a specific reason:

- **Sample names are never used as paths.** `store_id` is `<index>-<hash>`.
  Fidelity rule 11 documents that real names contain `"`, `$`, backticks and
  backslashes; paths add `/`, `:`, leading dots, macOS NFC/NFD, the 255-byte
  component limit, and case-insensitive collision on APFS.
- **Writes go to a temp file then rename**, so a reader never sees a partial
  parquet file.
- **`arrow`'s own CPU thread pool is clamped** by `sge_clamp_arrow_threads()`
  under `_R_CHECK_LIMIT_CORES_` and by the `spatialGE.cores` override. It sizes
  itself from `detectCores()` independently of anything here, so it is otherwise a
  second unguarded source of threads (fidelity rule 12).
- **No arrow handle is ever stored in the object or crossed over a fork.** Only
  paths are held; a handle is opened per read.
- **`portable` is FALSE for a store under `tempdir()`.** A backed STlist is not
  self-contained — `saveRDS()` records paths — so the flag marks an object that
  cannot survive the session.

**What this buys and does not buy.** It buys capacity: an STlist larger than RAM
can be held and read in gene-subset slices, and the in-R footprint of the fixture
drops from 888 KB to 72 KB. It does **not** yet speed up the aggregate paths —
`colSums`/`rowMeans`-style work still materializes and computes in memory rather
than pushing down into arrow queries. That is follow-on work, and until it lands
a backed STlist is slower than an in-memory one on those paths.

Parity is the property under test: `tests/testthat/test-store_parquet.R` runs
`SThet`, `STclust` and `summarize_STlist` both ways and asserts identical
results, and round-trips every matrix asserting `identical()`.

## Call graph / file layout

Post Phase-8 split, each analysis family is `entry -> prepare -> compute ->
assemble` across separate files rather than one monolith:

- `STdiff.R` (entry/validation) -> `STdiff_prepare.R` (VST + grid construction)
  -> `STdiff_fit.R` (the two parallel fit phases) -> `STdiff_assemble.R`
  (bind_rows + column contract).
- `STList.R` (entry: `STlist()`, `filter_data()` callers, etc.) ->
  `stlist_readers.R` (per-format import) -> `stlist_process.R` (post-import
  normalization) -> `sparse_utils.R` (shared sparse-matrix helpers).
- `detect_input.R` -> `detect_input_cases.R` (one named predicate per input
  type). Deliberately keeps **two independent, non-overlapping input-type
  vocabularies** side by side (one `STlist()` branches on, one it doesn't) --
  see the source comments before trying to unify them.
- `filter_data.R` -> `filter_data_stages.R` (5 named filter stages).
- `STenrich.R` -> `stenrich_helpers.R`; `gene_interpolation.R` ->
  `gene_interpolation_helpers.R`. Same split, smaller scope.

Every former `mclapply` anonymous closure is now a named, top-level, directly
callable unit taking explicit arguments instead of closing over the whole
`STlist` (Phase 7):

| Family | Named units |
|---|---|
| `transform_data` | `log_transf_one_sample`, `sct_transf_one_sample`, `gene_meta_one_sample` |
| `STclust` | `stclust_topgenes_one_sample`, `stclust_one_sample` |
| `SThet` | `sthet_stat_one_sample` |
| `STgradient` | `stgradient_one_sample` -> `stgradient_dists_one_sample`, `stgradient_expr_one_sample`, `stgradient_fit_one_gene` |
| `STdiff` | `stdiff_nonspatial_one_sample`, `stdiff_spatial_one_cluster` (fits themselves: `non_spatial_de()`/`spatial_de()` in `STdiff_fit.R`) |
| `STenrich` | `stenrich_perm_one_combo`, `stenrich_score_one_combo` |
| `gene_interpolation` | `interp_one_sample` -> `interp_krige_one_gene` |
| `STList` importers | `read_one_visium_sample`, `read_one_xenium_sample`, `read_one_cosmx_slide`, `read_one_counts_file`, `read_one_coords_file` |

This is what makes `bench_units.R` (below) and targeted `profvis` calls
possible -- point either at `stclust_one_sample()` directly instead of the
whole `STclust()` pipeline.

## Parallelism model

Single choke point: `R/utils_parallel.R`.

- `resolve_cores(n, cores=NULL, coerce=)` -- the boilerplate every call site
  used to duplicate (Windows forces 1 core; `cores=NULL` auto-detects via
  `count_cores(n)`; otherwise coerces the user-supplied value per `coerce`).
  Always call this instead of re-deriving cores by hand.
- `sge_lapply(X, FUN, ..., cores, preschedule)` -- forwards straight to
  `parallel::mclapply` rather than reimplementing its dispatch, so `mclapply`'s
  own semantics (serial `lapply` fallback under `cores < 2` or
  `length(X) < 2`, `.check_ncores()`, `mc.reset.stream()`) are inherited, not
  duplicated.
- `sge_bpparam(cores)` -- the `BiocParallel` equivalent; returns
  `SerialParam()` whenever cores resolve below 2, since a stray
  `MulticoreParam(workers=1)` still forks.

**Global override:** `options(spatialGE.cores=1)` or env var `SPATIALGE_CORES`,
checked by `sge_cores_override()` and honored by all three functions above.
Setting either forces every parallelized call site serial and in-process, with
no source edits -- this is what makes `profvis()` see past `mclapply` into the
named units above.

All five `sge_lapply` call sites in `R/stlist_readers.R` (lines 129, 136, 299,
470, 558) pass `cores=cores` and `preschedule=F`, with `STlist()` resolving
`cores` once at `R/STList.R:145` and threading it down. An earlier version of
this document claimed two of them passed no `mc.cores` and therefore fell
through to `getOption("mc.cores", 2L)`; that stopped being true with the "Fix
cores being silently ignored in STList.R and transform_data.R" change, and
fidelity rule 3 below is corrected accordingly. `read_matrices_fps()` and
`read_cosmx_input()` still never call `resolve_cores()` themselves, so a direct
call with `cores=NULL` would fall through to `mclapply`'s own default -- but
`STlist()` never does that.

## Naming conventions

- `sge_*` -- the parallelism-abstraction layer (`R/utils_parallel.R`).
- `resolve_*` -- boilerplate-unification helpers that resolve one argument
  from possibly-several sources (`resolve_cores`, `resolve_samples`).
- `assert_*` -- validation helpers that `stop()`/`raise_err()` on failure
  (`R/utils_validate.R`).
- `*_one_sample`, `*_one_cluster`, `*_one_combo`, `*_one_gene` -- Phase 7 named
  units of work, one per unit of parallelization.
- `sample_*` -- Phase 7 slot accessors (`R/utils_accessors.R`).
- `bl_*` -- baseline/regression-harness helpers (`dev/baseline/`), never
  shipped in the package (Rbuildignored).
- `stlist_digest()` -- the one function that walks every `STlist` slot in
  fixed order for golden comparison; extend it, don't duplicate it, if a new
  slot needs baseline coverage.

## Profiling switches

`dev/profiling/` (git-tracked, Rbuildignored):

- `profile_all.R` -- sets `options(spatialGE.cores=1)`, then runs
  `profvis::profvis()` per exported function against the fixture; HTML output
  under `dev/profiling/out/` (gitignored). Plot-returning functions are run
  through `force_plots()` (calls `ggplot_build()` on every returned plot) so
  rendering cost is actually captured.
- `rprof_summary.R` -- `Rprof()`/`summaryRprof()` text tables per pipeline
  step, diffable across commits, under `dev/profiling/out/rprof/`. Uses
  `Rprof(memory.profiling=TRUE)` since some hot spots (full-matrix
  `expandSparse`) are allocation problems that don't show up in CPU time.
- `bench_units.R` -- `bench::mark()` timings of the Phase 7 named units
  directly, in isolation from the orchestration/parallelization layer.
- `coverage.R` -- `covr::package_coverage()` over the fixture harness plus the
  testthat suite, so undercovered code is known rather than assumed.
- `FINDINGS.md` -- measured hot spots and caveats, written from the output of
  the four scripts above. Update it if a Phase 10-style optimization changes
  the picture.

All four scripts require `options(spatialGE.cores=1)` (or
`SPATIALGE_CORES=1`) to be useful -- without it, every profile bottoms out in
`mclapply` again.

## Regression / baseline workflow

`dev/baseline/` (git-tracked, Rbuildignored) plus `tests/golden/` (gitignored
-- never committed, regenerated locally):

- `dev/baseline/pipeline.R` (`bl_pipeline(cores, verbose)`) -- runs every
  exported function over the fixture matrix and returns one named list of
  results.
- `dev/baseline/capture_baseline.R` -- runs `bl_pipeline()` at **both**
  `cores=1` and `cores=3` (3, not 2 -- `mclapply` short-circuits below 2
  elements, so 3 samples at 3 cores is the smallest chunking that's actually
  non-degenerate), and saves each result (or captured error) to
  `tests/golden/<fn>__cores<N>.rds`.
- `dev/baseline/compare_baseline.R` (`bl_compare()`) -- structural comparison
  against a captured golden object; returns `NULL` on a match or a located
  diff otherwise. Never `identical()` (compares attribute order/environment
  identity, both irrelevant here) -- uses `waldo::compare(tolerance=NULL)`
  after `normalize_for_baseline()` strips call/environment noise per type
  (`dgCMatrix`, `listw`/`nb`, `sf`, `sp`, ggplot internals, etc. -- see the
  fidelity rules and Phase 3 notes in `plans/this-is-an-r-velvety-music.md`
  for the exact per-type recipe if extending this).
- `tests/testthat/test-golden.R` -- gated on `Sys.getenv("SPATIALGE_GOLDEN")
  == "true"` (never runs on CRAN or by default); three tests -- `cores=1` vs
  golden, `cores=3` vs golden, and a `cores=1`-vs-`cores=3` invariance check.

**Regenerating the baseline** (only when a change is deliberately
output-changing, per the Phase 10/11 commit discipline below):

```sh
Rscript dev/baseline/capture_baseline.R
```

**Running the regression gate:**

```sh
Rscript -e 'devtools::load_all(); testthat::test_local()'                       # plain, every commit
SPATIALGE_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_local()' # golden, before merging
```

**Commit discipline this repo follows for output-affecting changes**
(established during the Phase 10/11 work, still the expectation going
forward): a change that moves any golden-covered output gets its own commit,
and that commit both makes the change *and* regenerates `tests/golden/` via
`capture_baseline.R` -- so `git log` always has a single commit to bisect to
for "when did this number change and why", and the working tree is never left
with a stale baseline.

## Fidelity rules

Load-bearing constraints on this parallelism model, each backed by reading
`parallel`'s source or a specific call site. Breaking any of these silently
changes output rather than erroring:

1. **Forked children discard the RNG.** `mc.set.stream()` deletes
   `.Random.seed` per child, so children re-seed from time+PID.
   `transform_data(method='sct')` is not reproducible at `cores>=2` with >=2
   samples once gene count exceeds `sctransform::vst`'s sampling threshold
   (3000); `pseudobulk_samples(calc_umap=TRUE)` is unseeded and multithreaded
   via `uwot::umap`. `STenrich` *is* safe -- its `set.seed(seed)` sits inside
   the worker before the only RNG consumer, same seed every combo row.
2. **Forked and serial error semantics genuinely differ.** Preschedule mode
   runs `try(lapply(...))` per chunk, so one failing element poisons every
   element sharing that fork with only a `warning()`; plain `lapply` (the
   `cores=1` path) aborts on the first failure instead. `non_spatial_de()`
   (`R/STdiff_fit.R`) now guards against this (Phase 11.12) by checking
   `class(res_mod) == 'HLfit'` before calling `summary.HLfit()`, but the
   underlying asymmetry itself is inherent to `mclapply` and is not something
   to "fix" — design fixture data so every fit converges instead.
3. **CORRECTED (was: "`stlist_readers.R`'s two `mclapply` calls take no
   `mc.cores`, on purpose").** All five reader call sites pass `cores=cores`
   today; `STlist()` resolves it once at `R/STList.R:145`. The rule as written
   would have argued against a change that had already landed -- see the
   parallelism-model section above for what is actually true.
4. **`mc.preschedule` is per-call-site and never standardized.**
   `preschedule=FALSE` has a third set of semantics (per-element
   `mcparallel`, errors confined to one element, plus a `child.res[[1]]`
   unwrap when `is.list(child.res)`).
5. **`resolve_cores` preserves value *and type***, not just fork count.
   `coerce='as_integer'` can produce `NA` (caller must check and raise);
   `coerce='ceiling'` returns a double, letting a bad input error from
   `ceiling()` itself rather than with a spatialGE message -- this is
   intentional call-site-preserving behavior, not a bug. `count_cores()`
   returns `detectCores()/2` **unrounded**, and non-`mclapply` consumers
   (`calculate_vst(cores=)`, `create_listw_from_dist(cores=)`,
   `MulticoreParam(workers=)`) rely on that non-integral double surviving.
   `STenrich()` re-derives `cores` multiple times from different unit counts
   within one call -- it is not one value in that function, don't collapse it.
6. **Never relocate a `<<-`.** Any `<<-` in this codebase resolves to its
   enclosing function's own frame; moving that body into a top-level helper
   silently redirects the assignment into the package namespace or
   `globalenv()` instead of erroring.
7. **Sample order is user-visible.** It flows into `mclapply`'s `X`, then into
   `rbind`/`bind_rows` row order, then into how ties are broken by
   `adj_p_value`. Each call site's order source (`names(x@spatial_meta)` vs
   `names(x@counts)` vs `unique(combo[[1]])` vs the user's `samples=`) must be
   preserved exactly, not just made internally consistent.
8. **Column position and type are API.** `add_column(.before=)`,
   `relocate(.after=)`, and early-row type promotion in incremental
   `bind_rows` accumulation (starting from an empty `tibble()`, with a column
   that's logical `NA` in some rows and double in others) are all real
   contract, not incidental.
9. **`expandSparse` subsetting must use a logical mask, not character
   indexing** (`m[genes, ]` reorders to `genes`' order instead of the matrix's
   own row order). `data.frame()` also errors on duplicate row names, so
   full-matrix vs subset densification can differ in whether an error occurs;
   `check.names=FALSE` is load-bearing throughout.
10. **Timing output (the `cat`ted `difftime` timers) is inherently
    non-reproducible.** Compare return values only in any regression check.
11. **The worker progress messages are 18 inline `system(sprintf('echo "%s"',
    msg))` calls, not a helper.** There is no `sge_worker_msg` function -- an
    earlier version of this document named one that never existed. They sit
    inline in `transform_data.R`, `STclust.R`, `STdiff_fit.R`,
    `stenrich_helpers.R`, `gene_interpolation_helpers.R` and
    `stlist_readers.R`. Each shell-interpolates `"`, `$`, backticks, `\` --
    real sample names can contain these. Also writes to *process* stdout, which
    `capture.output` cannot see from inside a fork -- don't build any harness
    around stdout
    for this.
12. **`_R_CHECK_LIMIT_CORES_` makes `.check_ncores()` `stop()`/`warn()` when
    more than 2 cores are requested.** R CMD check sets this on CRAN's check
    machines. `count_cores()` (Phase 11.12) clamps its auto-detected value to
    at most 2 whenever this env var is active and non-`"false"`, so a
    many-core dev/CI machine can't trip a check that only exists to catch
    exactly that. Only the auto-detect path is clamped; an explicit
    user-supplied `cores=` is left alone. Baseline capture and comparison
    must run with this env var set identically on both sides, or diffs
    against golden are spurious rather than real.
13. **Four code paths compute a value that depends on spot *order*, not just on
    spot content.** Established by `dev/xver/30_order_probe.R`, which permutes
    spot order inside a single build and re-runs -- so these are properties of
    the algorithms, not of any refactor, and they were equally true in 1.2.2.
    Magnitudes measured on the standard fixture:

    | Path | Mechanism | Moves by |
    |---|---|---|
    | `STenrich` `p_value`/`adj_p_value` | `computeSubsampleSums` draws `Rcpp::sample()` over row *positions* of `coords_mat`, so a fixed `seed` picks different spots | full p-value resolution |
    | `SThet(k=)` `moran_i` | `knearneigh` breaks tied k-th distances by position; a regular grid is saturated with such ties | ~2e-2 |
    | `STgradient(distsumm='avg')` Spearman | averaged distances tie exactly; a ~1e-16 summation-order change flips a tie, which flips a rank | ~4e-3 |
    | `plot_counts()` density | `set.seed(subset_seed)` then `df[sample(1:nrow(df), ...), ]` -- again positions | visibly |

    Consequences for this repo: (a) do not treat a diff in these fields as a
    regression without running the order probe first; (b) `SThet()` without `k`
    and `STgradient(distsumm='min')` are the order-insensitive variants
    (agreement ~1e-15), so prefer them in any new fixture-based test that needs
    a stable number; (c) `@counts`, `@tr_counts` and the `vst.*` statistics are
    genuinely bit-identical under reordering when matched by name, so a diff
    there *is* real.

## Cross-version comparison harness

`dev/xver/` (untracked, not part of the package) compares an installed `main`
against an installed `alex-refactor` over 17 scenarios and classifies every
differing leaf path. It exists because `tests/golden/` can only ever compare this
tree against *itself*: the baseline was captured from the same working tree before
the `R/` edits, so it is blind to anything already baked into that capture, and it
reaches only one input format, one parameter set per function, three samples, and
no warnings.

Read `dev/xver/REPORT.md` for results and `dev/xver/README.md` for the design.
Two findings there were load-bearing for `dev/baseline/`, and both are now
addressed -- recorded because the reasoning still constrains any extension:

- **`bl_compare_stlist()` used to compare only 2 of the 10 `@misc` keys**
  (`sthet$listws` and `gene_krige`), so a change to `scaled_pbulk_mtx` -- the
  direct input to the `pseudobulk_heatmap` gene ranking -- or to `sp_images`, the
  entire output of `load_images()`, was invisible. All ten are compared now, and a
  `misc_uncompared_keys` guard fails the gate when a key appears with no
  comparison branch -- so adding a `@misc` key without also adding it here and to
  the census above is an error rather than a silent hole.
- **`pseudobulk_heatmap()` returns a `HeatmapList`** (it ends in
  `ComplexHeatmap::draw()`), so the displayed matrix is at
  `@ht_list[[1]]@matrix`, not `@matrix`. `bl_compare()` now has a
  `bl_compare_heatmap_signature()` branch comparing that matrix plus the
  dendrogram-determined row/column order, rather than falling through to `waldo`
  on grobs and environments -- noise that is part of why the inverted-CV bug went
  unnoticed for two releases. That branch is currently **unreachable on the
  standard fixture**: `pseudobulk_samples()` requires >= 4 samples and the fixture
  has 3, so all three `pseudobulk_*` steps capture an error string instead of
  output.

**The golden gate is same-machine, same-toolchain, and bit-exact**
(`tolerance = NULL` throughout). The pipeline reaches BLAS/LAPACK via `prcomp`,
`stats::loess`, `gstat` kriging, `spaMM::fitme` and `uwot::umap`, so results
differ in the low bits across BLAS implementations and across version bumps in
those dependencies. That is why `tests/golden/` is gitignored rather than
committed, and why the gate must not run on a CI runner whose platform differs
from the one that ran `capture_baseline.R` -- it would fail on float noise, and
the only way to green it would be loosening the tolerance, which removes the
drift detection that is the entire point. The platform-portable net is the
committed `tests/testthat/test-numeric-*.R` tier, asserting literal values at
explicit tolerances.

Captured results are stored **reduced**: `bl_reduce()` replaces each `ggplot` and
`HeatmapList` with exactly the signature `bl_compare()` compares, taking the
baseline from 13.9 MB to 4.0 MB. The reduction is idempotent and applied to both
sides, so a golden file captured before it existed still compares correctly. The
tradeoff -- widening a signature function requires re-running
`capture_baseline.R`, because a reduced file cannot answer a question its
signature does not already capture.
