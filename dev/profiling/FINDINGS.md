# Phase 9 profiling findings

> **STATUS: historical record of the pre-Phase-10 state. The line numbers are
> stale.**
>
> All six "Phase 10 targets" below have shipped (`f2e414a`, `74e5935`,
> `96bf081`, `8090078`, `1ca7962`, `330a69d`), so every `file:line` here points
> at code that has since moved -- the duplicate `cor.test` is cited at
> `R/STgradient.R:433` and following, but `spearman_cor_safe()` now exists at
> `R/STgradient.R:399`. Read this for the reasoning and the measured numbers,
> not as a current map.
>
> **It also has a blind spot worth naming, because it is why this file missed
> the package's worst scaling defect.** Every measurement here was taken on the
> standard fixture: 3 samples, 120 spots, ~200 genes. At that size a quadratic
> term is invisible. `create_listw_from_dist()` -- `SThet()`'s *default*
> weighting -- builds a complete graph, and had no benchmark and no test
> coverage at all. See "Scaling defects invisible at fixture size" at the end of
> this file.

Generated from three complementary tools, all run against the synthetic
fixture (`inst/testdata/make_fixtures.R`):

- `dev/profiling/profile_all.R` + `dev/profiling/rprof_summary.R` -- sampling
  profiles (`Rprof`) of the full multi-sample pipeline, one file per pipeline
  step, in `dev/profiling/out/rprof/*.txt`.
- `dev/profiling/bench_units.R` -- `bench::mark()` timings of the individual
  Phase 7 worker functions in isolation (no orchestration/parallelization
  layer), in `dev/profiling/out/bench_units.csv`.
- `dev/profiling/coverage.R` -- `covr` line coverage of `bl_pipeline()` (both
  `cores=1` and `cores=3`) plus the testthat suite, in
  `dev/profiling/out/coverage.txt`.

Caveat that applies throughout: the synthetic fixture is small (a handful of
samples/spots/genes), so the rprof sampling profiles are noisy -- several
steps have `sampling.time` under 0.05s, i.e. under 50 samples total, which is
too few for the by-self breakdown to reliably attribute cost to a specific
call. Where that's the case below, the `bench_units.csv` numbers (fixed
iteration counts, not sampling) are the more trustworthy source, and the
rprof number is noted only as corroboration or is explicitly flagged as
inconclusive.

## Coverage summary

Overall: **59.41%** line coverage (`dev/profiling/out/coverage.txt`).

The lowest-covered files are almost entirely visualization/reporting code
that the profiling/test harnesses don't call because they produce plots
rather than data (`plot_spatial_geneset.R`, `plot_spatial_meta.R`,
`pseudobulk_samples.R`, `distribution_plots.R`, `plot_image.R`,
`plot_counts.R`) plus alternate-input-format detection paths exercised by
real files rather than the in-memory fixture (`detect_input_cases.R`,
`detect_input.R`, `import_visium.R`, `import_xenium.R`). None of the six
Phase 10 targets are in zero-coverage code -- every target function is
exercised by the fixture, so the timings below reflect real executions, not
untested code paths.

One incidental finding: `STenrich()` fails on the synthetic fixture
(`dev/profiling/out/rprof/STenrich.txt` shows `ok: FALSE`, `` error: `.data`
must be a data frame without row names. ``). This is a pre-existing bug
unrelated to Phase 9/10 (already known from `test-golden.R`); it means
`bench_units.R`'s two `stenrich_*` benchmarks had to reconstruct `STenrich()`'s
internal combo/coords/delayed-matrix setup by hand rather than calling
`STenrich()` itself. Worth carrying into Phase 11's bug-fix pass.

## Phase 10 targets

### 1. `STclust.R` hclust hoist

- **Current location:** `R/STclust.R:361`, inside `get_hier_clusters_ks()`'s
  `for(k in ks) {...}` loop. `weighted_dists[[w]]` does not change across
  `k`, so `hclust()` (or the surrounding distance-matrix setup) is
  recomputed once per `k` for no reason. Contrast with
  `get_hier_clusters_dtc()` (`R/STclust.R:317-320`), which calls `hclust()`
  once per weight with no inner loop, since DynamicTreeCut determines the
  cluster count itself.
- **Measured:** `bench_units.csv` timed both with `ws=0.025`:
  `get_hier_clusters_ks(ks=c(2,3))` -> median 2.07ms/call;
  `get_hier_clusters_dtc()` -> median 5.03ms/call (dtc is slower overall
  because `cutreeDynamic`'s algorithm has more overhead than `cutree`, not
  because of the loop -- so the two numbers aren't directly comparable as a
  before/after).
- **Assessment:** confirmed structurally redundant, but the fixture doesn't
  isolate a standalone `hclust()`-only cost, so the *absolute* savings on
  this fixture can't be quoted precisely. What's clear is that the
  redundancy scales linearly with `length(ks)`: the benchmark above only
  used `ks=c(2,3)` (1 redundant call); real usage often passes a wider `ks`
  range (e.g. `2:10`), which would mean 9 redundant `hclust()` calls per
  weight instead of 1. Worth doing regardless of this fixture's small
  absolute number, since the fix (hoist `hclust()` before the `for(k in ks)`
  loop, index into the resulting dendrogram per `k`) is low-risk and
  output-neutral.

### 2. `expandSparse` subsetting before decompression

Four confirmed live call sites decompress a full sparse sample matrix with
`expandSparse()` and only subset to the genes actually needed *after*
decompression -- the opposite of `R/plot_spatial_expression.R:75`, which
already does this correctly (`expandSparse(counts[[i]][rownames(counts[[i]])
%in% genes, , drop=FALSE])`) and can serve as the template for the fix.

- `R/STdiff_fit.R:156` (`stdiff_run_spatial_tests`) and `R/STdiff_fit.R:202`
  (`stdiff_nonspatial_one_sample`).
- `R/STgradient.R:367` (`stgradient_expr_one_sample`).
- `R/plot_spatial_geneset.R:48` (inside a `for(i in samples)` loop; gene-set
  membership is only checked afterward, per gene set, inside the inner
  loop -- subsetting to the union of all requested gene-set genes before
  `expandSparse()` would shrink every per-set operation, not just avoid one
  redundant full decompression).

Ruled out: `R/filter_data_stages.R:135` (`filter_stage_counts`) also calls
`expandSparse()` without subsetting, but genuinely needs the full matrix --
it computes per-spot total counts/genes across *all* genes, so there's
nothing to subset to. Not a target.

- **Measured:** `expandSparse` itself is cheap in isolation (`bench_units.csv`:
  median 0.35ms/call, 927 itr/sec) -- the cost is the dense allocation size,
  which scales with genes-in-sample x spots-in-sample regardless of how few
  genes the caller actually needs. `stgradient_expr_one_sample` (which wraps
  one such call) measures at 9.27ms median, 2.8MB/call; `stdiff_nonspatial_one_sample`
  measures 753ms median / 26.2MB per call, and its variance across genes is
  driven substantially by this same full-matrix decompression happening
  once per sample rather than once per (sample, gene-subset).
- **Assessment:** confirmed, output-neutral fix (subset the sparse matrix by
  a `rownames(...) %in% wanted_genes` mask before `expandSparse()`, exactly
  as `plot_spatial_expression.R:75` already does). Savings scale with how
  small the wanted-gene subset is relative to the full gene panel.

### 3. `gene_interpolation.R` prediction grid

- **Current location:** `R/gene_interpolation_helpers.R:48` (a `match()`
  call) and `:60` (`sf::st_make_grid()`), both inside
  `interp_krige_one_gene()` -- i.e. both currently run once *per gene*, even
  though the prediction grid geometry only depends on the sample's spot
  coordinates and the interpolation resolution, not on which gene is being
  interpolated.
- **Measured:** `bench_units.csv`: `interp_krige_one_gene` -> median 5.58ms/call,
  2.27MB/call (n=5 iterations, since kriging is comparatively expensive);
  `interp_one_sample` (the per-sample wrapper that calls
  `interp_krige_one_gene` once per gene) -> median 12.0ms/call for this
  fixture's small gene set. `dev/profiling/out/rprof/gene_interpolation.txt`'s
  sampling.time is only 0.093s total (too few samples to break down further
  reliably), but its top self-time entries (`dyn.load` 45%, `context_eval`
  16%, `lazyLoadDBfetch` 15%) are consistent with per-call package/namespace
  dispatch overhead rather than the kriging math itself -- exactly the kind
  of fixed per-call cost that hoisting the grid construction out of the
  per-gene loop would amortize across genes instead of paying it every time.
- **Assessment:** confirmed target. With real gene panels (hundreds to
  thousands of genes vs. this fixture's handful), hoisting `match()` and
  `st_make_grid()` to run once per sample instead of once per gene should
  produce a much larger relative win than this fixture can demonstrate.

### 4. `err.csv` memoization

- **Current location:** `R/utils.R:327-338` (`raise_err()`). Every call
  re-reads `err.csv` from disk via `list.files(system.file(...), pattern=
  'err\\.csv')` + `read.csv()` -- there is no caching across calls in the
  same session.
- **Measured:** `raise_err()` is not exercised by the fixture or test suite
  (`R/utils.R:337` appears in the coverage.txt zero-coverage listing), so
  there's no timing data from this profiling pass. It's only reached on
  package-internal validation failures (error paths), so it is not a
  per-iteration hot path in normal successful runs -- the cost is one
  `read.csv()` of a small CSV, paid once per *raised error*, not once per
  sample/gene/etc.
- **Assessment:** legitimate but low-priority relative to the other five
  targets: it's simple and safe to do (memoize the CSV in an environment on
  first read), but it only matters for code paths that repeatedly raise and
  catch spatialGE errors in a loop (e.g. validation retried across many
  samples), which isn't exercised anywhere in the current fixture or
  baseline pipeline. No measured evidence of it being a real hot spot here;
  worth doing for correctness-of-intent (Phase 10 named it) but not
  expecting a measurable win on typical single-error usage.

### 5. `STList.R` -> `stlist_process.R` colSums densification

- **Current location:** `R/stlist_process.R:60-61`, inside `process_lists()`
  (called from `STlist()`):
  ```r
  coords_df_list[[name_i]][['total_counts']] = as.vector(colSums(as.matrix(counts_df_list[[name_i]])))
  coords_df_list[[name_i]][['total_genes']] = as.vector(colSums(as.matrix(counts_df_list[[name_i]]) != 0))
  ```
  Each line densifies the sparse counts matrix via `as.matrix()` (two
  separate dense copies, one per line) purely to call base `colSums()`;
  `Matrix::colSums()` operates on the sparse representation directly.
- **Measured:** `dev/profiling/out/rprof/STlist.txt`'s sampling.time is only
  0.005s (5 samples total -- too small to be conclusive on its own), but
  `process_lists` is nonetheless the single largest self-time entry in that
  profile (0.001s / 20% self, and it's the only frame with `total.pct=100`
  in the by-total table, i.e. everything sampled during `STlist()` traces
  through it). Consistent with, but not strong independent confirmation of,
  the target.
- **Assessment:** confirmed by code inspection (two avoidable dense
  conversions of the full counts matrix per sample, at `STlist()`
  construction time) even though the fixture is too small for rprof to
  isolate the cost precisely. Fix is mechanical and output-neutral: replace
  both `as.matrix(...)` + `colSums(...)` with `Matrix::colSums(...)`
  directly on the sparse matrix (and its `!= 0` variant), dropping both
  dense copies. Impact scales with matrix size, so this matters far more on
  real gene panels (thousands of genes) than on the fixture.

### 6. `STgradient.R` duplicate `cor.test()`

- **Current location:** `R/STgradient.R:433/438`, `:459/464`, `:486/491`,
  all inside `stgradient_fit_one_gene()`'s three mutually-exclusive
  branches (`out_rm & !robust`, `robust`, plain). Read the full function
  body this pass (`R/STgradient.R:405-498`) to confirm: only *one* of the
  three branches executes per call (they're `if`/`else if`/`else`), so this
  is not three redundant executions per call -- it's the same
  probe-then-retry pattern duplicated three times as source code:
  ```r
  cor_res = tryCatch({cor.test(x, y, method='spearman')}, warning=function(w){return(w)})
  if(any(class(cor_res) == 'simpleWarning')){
    if(grepl('standard deviation is zero', cor_res$message)){ pval_warn = 'zero_st_deviation' }
    cor_res = cor.test(x, y, method='spearman', exact=F)
  }
  ```
  The actual "duplicate work" is within a single call: `cor.test()` runs
  once to *probe* for a tie/zero-variance warning, then a second time with
  `exact=F` if the probe warned -- i.e. up to 2x `cor.test()` calls per
  gene, but only when ties are present (common with sparse/zero-inflated
  spot-level expression, where many spots share distance-to-reference or
  expression values).
- **Measured:** `bench_units.csv`: `stgradient_fit_one_gene` -> median
  1.02ms/call, 938 itr/sec on this fixture (whose data likely doesn't
  trigger the warning branch often, so this mostly reflects the single-call
  cost). `dev/profiling/out/rprof/STgradient.txt` (sampling.time 0.444s,
  the most substantial rprof sample size among the six targets) shows
  `.External2` at self.time 6.53%/**total.time 47.30%** and `lm.wfit` at
  self.time 6.98%/total.time 8.33% -- i.e. nearly half of sampled time in
  the whole `STgradient` step traces through `.External2`, which is
  consistent with `cor.test`'s C-level tie-handling/exact-p-value
  computation (both the `lm()`/`lm.wfit` and `cor.test` paths are called
  from `stgradient_fit_one_gene`, and this function runs once per gene per
  sample).
- **Assessment:** this is a code-duplication cleanup opportunity (extract a
  single `spearman_cor_safe(x, y)` helper used by all three branches) more
  than a raw performance win, since only one branch runs per call. The
  actual perf lever is different from what "duplicate cor.test" suggests at
  first read: it's the unavoidable-by-design probe/retry-on-warning
  pattern, which already only pays the 2x cost when ties are actually
  present. Given `.External2`'s outsized total.pct in the STgradient
  profile, this whole function (not just the duplication) is worth keeping
  in view for Phase 10, but the concrete, safe action here is the
  deduplication refactor, not a behavior change to the retry logic itself.

## Other hot spots observed (not in the original Phase 10 list)

- **`stdiff_spatial_one_cluster`** (`bench_units.csv`): median **20.8s**,
  **27.4GB** allocated, for a *single* spaMM mixed-model fit
  (`test_type='mm'`) on the synthetic fixture -- by far the largest cost of
  any unit measured, roughly 4 orders of magnitude beyond every other
  target in this file. This path is only reached under `STdiff(test_type=
  'mm')`, which the golden baseline never exercises (it uses `'t_test'`), so
  it doesn't show up in any rprof step file. Not one of the plan's named
  Phase 10 targets, but given the magnitude, worth flagging for anyone
  reaching for `test_type='mm'` on real data -- the memory footprint in
  particular (27GB for one cluster on a tiny fixture) suggests spaMM's
  internal matrices scale poorly and may need chunking/downsampling
  guidance in documentation even if the algorithm itself isn't changed.
- **The whole `spaMM::fitme` family is DEFERRED, not scheduled.** Both
  `spatial_de()`/`stdiff_spatial_one_cluster()` (above) and its non-spatial
  counterpart `non_spatial_de()` (`R/STdiff_fit.R:325`) run one `spaMM::fitme`
  per gene x cluster. The cost lives inside spaMM's own RcppEigen internals, not
  in spatialGE's R loops, so there is no local dplyr/vectorization lever; and
  reproducing spaMM's numeric output exactly is a research task rather than a
  refactor. Phase 4 item 2 of `plans/this-is-a-package-fluffy-marble.md` is
  marked DEFERRED for this reason -- the actionable part of that item is
  `stdiff_mean_test()` only. Revisit only with a dedicated investigation;
  chunking/downsampling guidance is the likely first step.
- **`sct_transf_one_sample`** (`bench_units.csv`): median 0.60s/call,
  32.9MB/call -- an order of magnitude slower than every other per-sample
  transform (`log_transf_one_sample` is 8.7ms). Expected, since SCTransform
  fits a regularized model per gene, but worth knowing it's the dominant
  cost of the transform stage specifically (not a bug, just a
  characteristic of the method).
- **`stdiff_nonspatial_one_sample`** (`bench_units.csv`): median 753ms/call,
  26.2MB/call, 8 GC cycles in 5 iterations -- notably GC-heavy relative to
  its peers, consistent with the full-matrix `expandSparse()` call inside
  it (target #2 above) triggering repeated large allocations.

## Scaling defects invisible at fixture size

Everything above was measured at 120 spots. These were not, and could not have
been -- they need `inst/testdata/make_scale_fixtures.R`, whose presets bracket
the range where quadratic terms dominate (`small` 120, `visium` 4096, `xenium`
50000).

### `create_listw_from_dist()` builds a complete graph

`R/utils.R:210`. Measured, single sample, `cores=1`:

| preset | spots | wall (s) | `listw` size |
|---|---:|---:|---:|
| small | 120 | 0.00 | 1.1 MB |
| visium | 4096 | 2.67 | **1216.8 MB** |

1.2 GB for one sample at 4,096 spots -- 1,100x the memory for 34x the spots.
Worse than quadratic because the loop attaches **character names** to every
neighbour vector (`names(neighbours[[id]]) = idx[...]`), and computes membership
via `rownames(adj) %in% rownames(adj)[id]`, an O(n) string comparison inside an
O(n) loop. The result is a graph in which every spot neighbours every other.

Three things make this the top priority rather than a curiosity:

1. **It is the default.** `SThet(k=NULL)` takes this path (`R/SThet.R:118`).
2. **The documentation steers users onto it.** Fidelity rule 13(b) recommends
   `SThet()` without `k` for order-insensitivity, which is correct about
   stability and unaware of the cost.
3. **Extrapolation rules out real data.** At the `xenium` preset (50,000 cells)
   the same formula gives roughly 181 GB, so it cannot run at all -- that preset
   is a post-fix target, not a before/after comparison.

`sthet_stat_one_sample()` (`R/SThet.R:245`) then loops genes in R, calling
`spdep::moran`/`geary` once per gene against that graph. Replacing both halves
with one compiled all-genes pass over the coordinates avoids materializing the
graph entirely. Before claiming a speedup, confirm whether `spdep::moran`/`geary`
are themselves compiled -- it determines whether the win is the construction
alone or both halves.

Note `create_listw_from_knn()` (`R/utils.R:192`) hardcodes `cores = 1`, so its
`sge_lapply` never forks regardless of what the caller requests.

## Reproducing

```sh
Rscript dev/profiling/profile_all.R      # writes dev/profiling/out/rprof/*.rds
Rscript dev/profiling/rprof_summary.R    # writes dev/profiling/out/rprof/*.txt
Rscript dev/profiling/bench_units.R      # writes dev/profiling/out/bench_units.csv
Rscript dev/profiling/coverage.R         # writes dev/profiling/out/coverage.txt
```

All four scripts are run from the package root and write into the
gitignored `dev/profiling/out/`. All four require
`options(spatialGE.cores=1)` or `SPATIALGE_CORES=1`, or every profile bottoms
out in `mclapply`.

The scale measurements above are reproduced with:

```sh
Rscript -e 'devtools::load_all(); source("inst/testdata/make_scale_fixtures.R")'
```

then building an STlist per preset via `fixture_scale_stlist(preset=)` and timing
the unit directly.
