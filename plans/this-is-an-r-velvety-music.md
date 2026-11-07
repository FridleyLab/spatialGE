# spatialGE refactor: make the package profilable without changing behavior

## Context

`spatialGE` 1.2.2 is on CRAN. The analysis code works, but it is effectively
un-profilable, and the repository has no safety net for changing it.

Three structural facts cause the profiling problem:

1. **Almost all compute runs in forked children.** ~25 `parallel::mclapply` call
   sites across 12 files, every one passing an *anonymous* `function(i){ ... }`
   closure containing the entire per-sample algorithm. `Rprof`/`profvis` in the
   parent cannot see inside a fork, so `profvis(STclust(x))` shows flat time in
   `mclapply` and nothing else.
2. **There is no way to force serial execution globally.** Each function
   re-derives its own worker count from one of three subtly different copies of
   the same boilerplate (13 copies), so going in-process means editing 13 places.
   The source shows this being done by hand already: `R/SThet.R:172-174` and
   `R/STdiff.R:493-495` keep commented-out serial / `bettermc` / `pbmcapply`
   variants beside the live `mclapply`.
3. **The units of work have no names.** Even serial, a profile of `STgradient`
   bottoms out in one 290-line anonymous frame (`R/STgradient.R:118-408`).
   `STdiff` is a 591-line function doing validation, VST, grid construction, two
   parallel compute phases and two assembly phases with no internal boundary.

Layered on top: the test suite has **one assertion**, runs only when GitHub is
reachable, and there are no committed fixtures (no `data/`, no `inst/extdata/`).
There is currently no way to prove a refactor preserved behavior.

**Intended outcome.** Bit-identical output and an unchanged user API, plus: one
global switch that makes the whole package run in-process; named, independently
callable units of work at every hot spot; deterministic offline fixtures; and
golden-output regression tests that prove equivalence. Genuine bugs are
documented and fixed *afterwards* in separate commits, so the harness can always
tell an intended change from an accidental one.

All work on branch `alex-refactor`. No pushes. No mention of Claude or AI in any
commit, comment, or file.

---

## Phase 0 — Finish the environment (blocking, no commit)

I verified the conda environment. It is about half built:

- `/opt/anaconda3/envs/spatialge-dev`: **R 4.5.3**, `clang` 22, `gfortran`.
  (System R outside conda is 4.6.1 — use the conda R.)
- The hard-to-build geospatial stack is present: `sf`, `sp`, `spdep`, `gstat`,
  `arrow`, `s2`, `units`, `classInt`, `deldir`, `stars`.
- **Only 56 R packages total.** Missing, and also absent from the renv cache for
  R 4.5: `Matrix` (a hard `import()`), `BiocParallel`, `ComplexHeatmap`,
  `DelayedArray`, `DelayedMatrixStats`, `dynamicTreeCut`, `EBImage`, `GSVA`,
  `ggplot2`, `ggforce`, `ggpolypath`, `ggrepel`, `hdf5r`, `jpeg`, `khroma`,
  `readr`, `readxl`, `sfsmisc`, `spaMM`, `wordspace`, `sctransform`, `uwot`,
  `dplyr`, `tibble`, `tidyr`, `stringr`, `scales`, `png`, `RColorBrewer`,
  `RcppEigen`. `src/` cannot compile without `RcppEigen`/`RcppProgress`.
- renv **is** activated (`.Rprofile` → `renv/activate.R`) but has **no
  `renv/library/` and no `renv.lock`** — only an empty `renv/staging/1`. The
  autoloader **hangs** under the sandbox because it reaches for the network
  (`RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` starts R instantly, confirming it).
- **`parallel::detectCores()` returns `NA` in this sandbox** (`sysctl` is
  blocked). So `count_cores()` returns `NA` and `mclapply` hard-errors
  `'mc.cores' must be >= 1`. **With its default `cores=NULL`, spatialGE cannot
  run here at all.** Consequence: the harness must pass `cores=` explicitly on
  every call — which is better for reproducibility regardless.

Work:

1. `dev/environment/spatialGE-dev.yml` — conda spec pinning R 4.5.3, the
   toolchain, and system libs (gdal, geos, proj, hdf5, fftw, libtiff, libjpeg,
   libpng).
2. `dev/environment/Makevars` — flags this env needs (clang 22 + gfortran paths,
   `-fopenmp` on arm64 macOS); committed, applied via `R_MAKEVARS_USER`.
3. `dev/environment/install_deps.R` — CRAN + Bioconductor for every
   `Imports`/`LinkingTo`/`Suggests` entry. **Pin `ggplot2`**: the 4.x S7 rewrite
   changes `ggplot_build` layer-data columns, which the plot comparisons depend
   on, and `R/STdiff_volcano.R:118` (`element_rect(fill=NULL)`) is a known 4.x
   risk. Record `sessionInfo()` + `extSoftVersion()` + the BLAS in
   `dev/environment/RECORDED.md`.
4. `dev/environment/README.md` — activation, and the
   `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` workaround for offline sessions.

Needs network, which the sandbox blocks. Either run it yourself with
`! Rscript dev/environment/install_deps.R`, or I run it and you approve the
un-sandboxed prompt. **Everything after this depends on it.**

**renv decision.** You asked for `renv/` and `renv.lock` in `.gitignore` and
`.Rbuildignore`. Then `.Rprofile` must be ignored too, or a fresh clone gets a
`.Rprofile` sourcing a missing `renv/activate.R` and every R session fails. So:
`.Rprofile`, `renv/`, `renv.lock` all ignored; reproducibility comes from the
conda yaml + `install_deps.R`. The tradeoff is no exact version pinning on
clone — say the word if you'd rather track `renv.lock` instead.

---

## Phase 1 — Build hygiene (zero `R/` changes)

Create `alex-refactor`. Then:

- `git rm --cached src/*.o src/spatialGE.so` — these 4 files are **tracked**
  despite `.gitignore` listing `src/*.o` (ignore rules don't apply to tracked
  files). A fresh clone gets `.o` newer than `.cpp`, so `make` may skip
  recompilation and link a stale arm64 `.so` — which would silently invalidate
  every baseline comparison.
- `.gitignore` — add `.claude/`, `.Rprofile`, `renv/`, `renv.lock`, `dev/data/`,
  `tests/golden/`, `dev/profiling/out/`, `*.Rcheck/`, `*.tar.gz`, `.Rhistory`,
  `.RData`, `.DS_Store`.
- `.Rbuildignore` — add `^\.claude$`, `^dev$`, `^plans$`, `^renv$`,
  `^renv\.lock$`, `^\.Rprofile$`, `^tests/golden$`. **Add `^dev$` and
  `^tests/golden$` now, before those directories exist**, so it can't be
  forgotten. Also fix the missing trailing newline.
- `DESCRIPTION` — fix `Copyright: ins/COPYRIGHT` → `inst/COPYRIGHT`.
- **Record the current `R CMD check` output as the reference NOTE/WARNING set.**
  It is not clean today: `LazyData: true` with no `data/` directory is a warning.
  That set is what later phases must not grow.

---

## Phase 2 — Fixtures and harness (zero `R/` changes)

### Synthetic fixtures — generator at `inst/testdata/make_fixtures.R`

Seeded, deterministic, offline. Generated into `tempdir()` by
`tests/testthat/helper-fixtures.R` rather than committed as binaries — and
**never into `data/`**, since `data/` objects need `\docType{data}` pages or they
raise an "Undocumented data sets" NOTE, forcing exactly the `man/` churn we're
avoiding.

Shape: **3 samples** (not 1 or 2 — see the fidelity rules; `mclapply`
short-circuits to `lapply` whenever `length(X) < 2`, so fewer than 3 samples
never exercises chunking), ~200 genes × ~120 spots, seeded Poisson counts, a
spatial gradient in ~10 genes and a 2-domain structure so `STclust`, `SThet`,
`STgradient` and `STdiff` produce non-degenerate results.

**Genes must stay below `sct_n_regr_genes` (3000)** so `sctransform::vst` skips
its gene sampling — otherwise SCT is irreproducible (see fidelity rules).

Deliberate degenerate content, because that's where the branching is: one
all-zero gene; one constant gene (hits the `'standard deviation is zero'`
`cor.test` branch); one **zero-library-size spot** (hits the `zero_size` path in
`log_transf`/`sct_transf` and the slot surgery at `R/transform_data.R:88-97`);
one sample with differently-ordered barcodes (`read_list_dfs`/`process_lists`
matching); one gene name with a shell metacharacter and one with a
syntactically-invalid-name character (proves `check.names=FALSE` is preserved).

Input layouts covered synthetically:

- csv/tsv matrices, list-of-data.frames, `samplesfile_matrices`.
- Samples file as **both csv and xlsx** — separate `detect_input` branches, and
  the xlsx fixture is what will later prove the `//.xlsx$` fix.
- **MEX**: `Matrix::writeMM` + gzipped `features.tsv.gz`/`barcodes.tsv.gz` under
  literally-named `filtered_feature_bc_matrix/` and `cell_feature_matrix/`, plus
  `spatial/tissue_positions_list.csv` and `scalefactors_json.json`. This is the
  only way to demonstrate the `import_visium` shadowing bug.
- **10x h5** via `hdf5r` (`matrix/{data,indices,indptr,shape,barcodes,features/*}`).
  Not committed — hdf5 isn't byte-reproducible.
- **`cells.parquet`** via `arrow::write_parquet`; **CosMx** `exprMat`/metadata as
  plain csv; images via `png::writePNG` of a small array.

**Dev/network tier only:** real 10x/NanoString file quirks, Seurat-object input,
realistic-scale GSVA/`spaMM`/kriging (also the only meaningful profiling
targets), anything over a few seconds.

### Test wiring

- `tests/testthat/test-tests.R`'s top-level `download.file` moves *inside*
  `test_that()` behind `skip_on_cran()` + `skip_if_offline()` (curl is already in
  Suggests). Today, when the download fails the whole block is skipped and
  **zero** tests run while check still passes.
- Golden-comparison tests additionally gate on
  `Sys.getenv("SPATIALGE_GOLDEN") == "true"`, with `.rds` under an Rbuildignored
  `tests/golden/` — they must not ship in the tarball, or CRAN runs them against
  a dependency set we don't control.
- All skips are the first lines *inside* `test_that()`, never at file scope.

---

## Phase 3 — Golden baseline (the pivot commit, zero `R/` changes)

`dev/baseline/capture_baseline.R` runs every exported function over the fixtures
across an argument matrix, at **`cores=1` and `cores=3`** separately, writing one
`tests/golden/<fn>__cores<N>.rds` per case so a failure names the culprit.

`set.seed()` immediately before **every** call, never once at script top — a
`set.seed` inside a fork doesn't propagate to the parent, so global RNG state
after `cores=1` and `cores=3` runs genuinely differs.

Record up front which outputs are **not** bit-comparable, so later diffs don't
chase them: SCT above 3000 genes, `calc_umap=TRUE`, and error paths.

### Comparison mechanics

`saveRDS(version = 3)`. **Never `identical()`** on round-tripped objects — it
compares attribute order and environment identity. Use
`waldo::compare(tolerance = NULL)` for exact structural comparison with a located
diff; `all.equal` is wrong here (returns a message, tolerates class changes on
lists).

One `normalize_for_baseline()` recurses and (a) strips `attr(,"call")`, `terms`,
`srcref`; (b) replaces environment/closure leaves with a sentinel. Then per type:

- **`dgCMatrix`** — `Matrix::drop0()` both sides first (the "makeSparse on an
  already-sparse matrix" paths differ in explicitly stored zeros), then compare
  `@i`, `@p`, `@x`, `@Dim`, `@Dimnames`, and assert `class()` (`dgC` vs `dge` vs
  `ngC` is a real regression).
- **`listw`/`nb`** — `$neighbours`, `$weights`, `$style`, with `region.id`/`sym`
  asserted separately. `attr(,"call")` **will** change when the call site moves,
  and these objects are persisted into `x@misc[['sthet']][['listws']]`
  (`R/SThet.R:123`), so stripping it is mandatory.
- **`sf`** — WKT via `st_as_text(st_geometry(x))` + `st_crs(x)$wkt` +
  `st_drop_geometry(x)`; raw `sfc` comparison fails on `bbox` float noise.
- **`sp::SpatialPointsDataFrame`** — `sp::coordinates()` and `@data`.
- **`spaMM` HLfit** — never serialize. `STdiff` already reduces models to
  p-values before returning, so the returned tibbles *are* the contract.
- **ggplot** — never compare the object. Per plot: `names()` of the returned list
  (the keys are API — `paste0(gene, "_", i)`),
  `lapply(ggplot_build(p)$data, as.data.frame)` (the real payload, and it forces
  evaluation), `p$labels`, scale `limits`/`breaks` and scale classes, and
  `vapply(p$layers, \(l) class(l$geom)[1], "")`.
- **`STlist`** — a `stlist_digest()` walking slots in fixed order, each
  normalized as above.

---

## Phase 4 — Delete dead code (before any restructuring)

Deliberately before the utils unification: `R/SThet_invdist_test.R` contains two
of the 13 cores-resolution variants and a `create_listw_from_dist` caller, so
unifying first means writing variants that are then thrown away. And there is no
point refactoring code that's about to be deleted.

- **`R/SThet_invdist_test.R`** — entire 230-line file, zero callers, a
  near-verbatim clone of `R/SThet.R`.
- `create_listw` (`R/utils.R:183-198`), `calculate_euclidean_distances` and
  `calculate_sum_distances` (`R/STenrich.R:314-321`, `:423-433`),
  `process_sample_filepaths` (`R/STList.R:1008-1035`),
  `process_sample_names_from_file` (`R/STList.R:1092-1138`), commented-out
  `process_meta_geomx` (`R/STList.R:1192-1219`).
- ~150 lines of commented-out alternatives interleaved with live code
  (`R/SThet.R:172-174`, `R/STdiff.R:493-496`/`:532-535`,
  `R/STenrich.R:147-159`/`:254-259`, `R/gene_interpolation.R:74-86`/`:143-155`).
- Unreferenced `error0008`–`error0015` rows in `inst/err.csv`.

**Grep `man/` and `NAMESPACE` for every deleted symbol** — a documented-but-absent
object fails check. Watch for newly-unused `@importFrom` creating new NOTEs,
remembering we cannot re-run `document()`.

---

## Phase 5 — `R/utils_parallel.R`

```r
resolve_cores(n, cores = NULL, coerce = c("as_integer", "ceiling"))
sge_lapply(X, FUN, ..., cores, preschedule)
sge_bpparam(cores)
```

`sge_lapply` **forwards straight to `parallel::mclapply`** rather than
reimplementing the dispatch. I read R's `mclapply` source to confirm this is safe
— it short-circuits at line 26–27 with

```r
if (cores < 2L && is.null(affinity.list)) return(lapply(X = X, FUN = FUN, ...))
```

and independently when `length(X) < 2`. So `cores = 1` already means "plain
`lapply`, in-process, visible to `profvis`". Letting `mclapply` decide means we
inherit its exact semantics — including `.check_ncores()`, `mc.reset.stream()`,
and the clamp of `cores` to `length(X)` — instead of replicating them.

The new capability is a **single global override** honored by `resolve_cores`:

```r
options(spatialGE.cores = 1)   # or env var SPATIALGE_CORES
```

That forces every call site serial and in-process with no API change and no source
edits, replacing the comment-out-the-mclapply workflow visible in the source
today. `resolve_cores` also absorbs the Windows clamp (8 copies) and is checked
*after* it.

`sge_bpparam` returns `BiocParallel::SerialParam()` when `cores < 2`, otherwise
`MulticoreParam(workers = cores)`. Without it, `R/STenrich.R:388` stays a second
source of invisible work — `MulticoreParam(workers = 1)` still forks.

---

## Phase 6 — `R/utils_validate.R` (separate commit from Phase 5)

`assert_stlist()`, `assert_transformed()`, `resolve_samples()`. Separate commit
because a sample-order regression is the single most likely silent breakage, and
it should be bisectable alone.

The ~15 existing copies are **not** uniform: they differ in slot read, numeric
index handling, whether substring `grep` matching is used, and exact error
strings. `resolve_samples()` takes explicit parameters so each site keeps its
current behavior verbatim — **including the buggy substring `grep`**, which
changes *which* samples get selected. Tests assert each site's **ordered** sample
vector (not `setequal`) and its exact error string with `expect_error(fixed=TRUE)`.

---

## Phase 7 — Extract named units of work (one commit per function family)

The primary profiling win. Every `mclapply` anonymous closure becomes a **named
top-level internal function taking explicit arguments** instead of closing over
the whole `STlist` — readable stack frames, and each unit directly callable for
microbenchmarking one sample.

Order cheapest-first so harness bugs surface on fast functions before `STdiff`:
`transform_data`, `SThet`, `STclust`, `gene_interpolation`, `STgradient`,
`STenrich`, `STdiff`, then the `STList.R` importers.

| Current | Becomes |
|---|---|
| `R/transform_data.R:153`, `:226`, `:102` | `log_transf_one_sample`, `sct_transf_one_sample`, `gene_meta_one_sample` |
| `R/STclust.R:146`, `:162` | `stclust_topgenes_one_sample`, `stclust_one_sample` |
| `R/SThet.R:175`, `:250` | `sthet_stat_one_sample` (the two helpers differ only in the `spdep` call and target column) |
| `R/STgradient.R:118-408` | `stgradient_one_sample` → `stgradient_dists_one_sample`, `stgradient_expr_one_sample`, `stgradient_fit_one_gene` |
| `R/STdiff.R:280`, `:493` | `stdiff_nonspatial_one_sample`, `stdiff_spatial_one_cluster` |
| `R/STenrich.R:224`, `:332` | `stenrich_perm_one_combo`, `stenrich_score_one_combo` |
| `R/gene_interpolation.R:210` | `interp_one_sample` → `interp_krige_one_gene` |
| `R/STList.R:639`, `:799`, `:880`, `:493`, `:504` | `read_one_visium_sample`, `read_one_xenium_sample`, `read_one_cosmx_slide`, `read_one_counts_file`, `read_one_coords_file` |

Also here, since every call site is already being touched: **internal
non-exported slot accessors**, the precondition for any future slot-layout change
(e.g. making `gene_krige` sample-major like every other slot). A `validity()`
method is *not* free — it can start rejecting objects — and belongs after the
bug-fix phase.

---

## Phase 8 — Decompose the monoliths

Split by `validate → prepare → compute → assemble`. Deliberately after Phase 7,
so a pure file move isn't entangled with signature changes; verify with
`git diff -M --stat` showing renames plus baseline equality.

- `R/STdiff.R` (930 lines; `STdiff` alone is 591) → `STdiff.R` (entry) +
  `STdiff_prepare.R` + `STdiff_fit.R` + `STdiff_assemble.R`. This is where "is
  the time in `spaMM::fitme` or in `bind_rows`?" is currently unanswerable.
- `R/STList.R` (1220 lines) → `STlist.R` + `stlist_readers.R` +
  `stlist_process.R` + `sparse_utils.R`.
- `R/detect_input.R` (265-line function) → one named predicate per input type.
  It maintains **two independent, non-overlapping input-type vocabularies** that
  `STlist()` consumes simultaneously (`R/STList.R:228` branches on one, `:247` on
  the other); keep both, make the split explicit and documented.
- `R/filter_data.R` (209 lines, 19 parameters) → 5 named filter stages.
- `R/STenrich.R`, `R/gene_interpolation.R` — same, smaller scope.

Plot functions return **lists of unevaluated `ggplot` objects**, so rendering
happens outside any `profvis(STplot(...))` block. Rather than change the return
contract, the harness forces rendering via a `force_plots()` helper calling
`ggplot_build()` on every element — used in **both** the profvis block and the
comparison. Without it, profiling and regression coverage of ~8 plot functions
are both near-empty.

---

## Fidelity rules — the things that will silently break output

These come out of reading `parallel`'s source and the specific call sites. Each
is a hard constraint on the phases above.

1. **Forked children discard the RNG.** `parallel:::mc.set.stream()` deletes
   `.Random.seed` in each child under the default Mersenne-Twister, so children
   re-seed from time+PID. Therefore:
   - `transform_data(method='sct')` is **not reproducible at all** at `cores≥2`
     with ≥2 samples — `sctransform::vst` (`R/transform_data.R:243`) samples genes
     when `nrow > n_genes`. Mitigated by keeping fixture genes under 3000. Do
     **not** add a `set.seed` during the refactor; that's a behavior change for
     the bug phase.
   - `pseudobulk_samples(calc_umap=TRUE)` (`R/pseudobulk_samples.R:145`) calls
     `uwot::umap` unseeded and multithreaded. Baseline `calc_umap=FALSE`; assert
     only dims/class for that path.
   - `STenrich` **is** safe: `set.seed(seed)` at `R/STenrich.R:253` sits inside the
     worker before the only RNG consumer, with the same seed for every combo row,
     so results are independent of core count and chunking. Confirmed
     `computeSubsampleSums` uses R's RNG (`Rcpp::RNGScope` in
     `src/RcppExports.cpp:42`, `Rcpp::sample` at `src/spatialge_helpers.cpp:33`) —
     no independent C++ generator.
2. **Error semantics differ between the forked and serial paths, already today.**
   The preschedule path does `sendMaster(try(lapply(...), silent=TRUE))`, so **one
   failing element poisons every element sharing that core**, with only a
   `warning()`; `lapply` aborts instead. `R/STdiff.R:718` reaches this — when
   `fitme` errors, the handler assigns the *message string* to `res_mod` and
   `summary.HLfit` on a character then errors out through the closure. Don't
   preserve this: choose fixture data where every fit converges, and log the
   asymmetry for the bug phase.
3. **Leave `R/STList.R:294`/`:302` alone.** They call `mclapply` with **no
   `mc.cores`**, so they always use `getOption("mc.cores", 2L)`. Routing them
   through `resolve_cores` changes the fork count, can push `nc > 2` into
   `.check_ncores`, and changes which samples share a chunk. Keep verbatim through
   the refactor; fix in the bug phase.
4. **Forward `mc.preschedule` per call site.** `preschedule=FALSE` has a third set
   of semantics (per-element `mcparallel`, errors confined to one element, plus
   `if (is.list(child.res)) child.res <- child.res[[1]]` unwrapping). Never
   standardize it.
5. **`resolve_cores` must preserve the *value and type*, not just the fork
   count.** Two live variants: `as.integer(cores)` + the
   `'Could not recognize number of cores requested'` check; and
   **`ceiling(cores)` with no NA check** (`R/STclust.R:139`, `R/STenrich.R:124`,
   `R/STdiff.R:276`), which returns a **double** and lets `STclust(cores='a')`
   error from `ceiling` rather than with the spatialGE message. `count_cores`
   returns `detectCores()/2` **unrounded**, and that value is forwarded to
   non-`mclapply` consumers — `calculate_vst(cores=)`, `create_listw_from_dist(x,
   cores=)`, `MulticoreParam(workers=)` — so non-integral doubles must survive.
   Also: `R/STenrich.R:201`/`:212`/`:221` re-derive `cores` three more times from
   different unit counts, gated on `user_cores`. `cores` is not one value in that
   function; do not collapse it.
6. **Never relocate a `<<-`.** At `R/STclust.R:196-211` the `<<-` resolves to
   `STclust`'s own frame. Move that body into a top-level helper and it walks past
   into the package namespace or `globalenv` — silent corruption, no error.
   Convert to a pure `x <- ...; return(x)`. Same at
   `R/pseudobulk_samples.R:80-83`. That loop also mutates **sequentially** and its
   drop-existing-columns test reads state from iteration `w-1`, so building all
   columns and joining once changes results on column collisions.
7. **Sample order is user-visible.** It feeds `mclapply`'s `X`, which feeds
   `do.call(rbind, ...)` row order (`R/STenrich.R:275`), and at `R/STenrich.R:293`
   ties in `adj_p_value` are broken by incoming row order. Enumerate each site's
   order *source* (`names(x@spatial_meta)` vs `names(x@counts)` vs
   `unique(combo[[1]])` vs the user's `samples`).
8. **Column position and type are API.** `add_column(.before='p_value')`
   (`R/STenrich.R:281`), `relocate(.after='exp_adj_p_val')` in `STdiff`. Any switch
   from incremental `bind_rows` accumulation to a single `bind_rows` must be
   `waldo`-checked including types and positions: `exp_p_val = exp_summ[[2,4]]`
   (`R/STdiff.R:583`) is logical `NA` in some rows and double in others, and the
   accumulator starts from an empty `tibble()` — early-row type promotion is the
   concrete failure mode.
9. **`expandSparse` hoisting is not trivially safe.** `R/plot_spatial_expression.R:82`
   already subsets with a **logical mask**, preserving matrix row order. Fixing
   the other sites with character indexing (`m[genes, ]`) would give the order of
   `genes` instead. Use the identical logical-mask expression. Also:
   `data.frame()` errors on duplicate row names, so full-matrix vs subset
   densification can differ in whether an error occurs; and `check.names = FALSE`
   is load-bearing.
10. **Timing output is inherently non-reproducible.** The 6 `cat` timers print
    `round(difftime(...), 2)`. Compare return values only; if stdout is compared
    at all, mask `\d+\.\d+ min`.
11. **`sge_worker_msg` must keep the `system(echo)` form verbatim** during the
    refactor. `system(sprintf('echo "%s"', msg))` shell-interpolates `"`, `$`,
    backticks and `\`, and real CosMx sample names can contain them. It also
    writes to *process* stdout, which `capture.output` cannot capture from a fork
    — so don't build the harness around stdout. (Note `R/STclust.R:79`'s
    `verbose = 1L` clobber makes the `verbose > 1L` branch in the STclust worker
    dead today, so refactoring it is free.)
12. **`_R_CHECK_LIMIT_CORES_`** makes `.check_ncores` `stop()` when `nc > 2`, so
    any `count_cores()` path errors under `R CMD check` on a many-core machine —
    a latent CRAN failure today, independent of this work. Baseline capture and
    comparison must run with it set identically or diffs are spurious.

---

## Phase 9 — Profiling harness

`dev/profiling/`:

- `profile_all.R` — `options(spatialGE.cores = 1)`, then `profvis::profvis()` per
  exported function against a fixture, HTML into `dev/profiling/out/`. Plot
  functions wrapped in `force_plots()`.
- `rprof_summary.R` — `Rprof` + `summaryRprof` text tables, diffable, for
  tracking regressions. Include `Rprof(memory.profiling = TRUE)`: full-matrix
  `expandSparse` is an allocation problem more than a time problem and won't show
  in a CPU profile.
- `bench_units.R` — times the individual named workers from Phase 7 directly, so
  one sample's cost is measurable without the surrounding orchestration.
- `coverage.R` — `covr::package_coverage(type="none", code=<harness>)`, the
  objective answer to "do the fixtures reach every branch".
- `dev/profiling/FINDINGS.md` — measured hot spots, so the optimization work
  starts from data rather than from reading.

Also add minimal CI (`.github/workflows/`, `r-lib/actions`) running the synthetic
tier. There is no CI today, so "check stays clean" is otherwise an assertion
rather than a verified property.

---

## Phase 10 — Output-neutral optimizations (one commit each)

These are the only changes that can move floating-point values, so each gets its
own commit and its own baseline run:

- Hoist `hclust` out of the k loop (`R/STclust.R:334`) — verify
  `get_hier_clusters_ks` doesn't reorder `weighted_dists` between k values.
- Subset before `expandSparse` (`R/STdiff.R:288`, `:480`, `R/STgradient.R:217`,
  `R/plot_spatial_geneset.R:57`) — logical mask only, per fidelity rule 9.
- Hoist the per-gene prediction grid (`R/gene_interpolation.R:231-235`). Safe
  because the `match()` at `:223` always returns all spots, so the grid is
  gene-invariant — but hoist the **grid**, not the `SpatialPointsDataFrame`
  (`sp::coordinates(df) = ~xpos+ypos` mutates in place).
- Memoize `inst/err.csv` (`R/utils.R:364`) **lazily on first call, not in
  `.onLoad`** — under `load_all()` `system.file()` resolves differently and the
  csv is edited during development. Preserve `header=FALSE, quote="'"` and the
  `suppressWarnings`.
- `Matrix::colSums` instead of `colSums(as.matrix(...))` twice per sample
  (`R/STList.R:975-976`).
- Collapse the duplicate `cor.test` (`R/STgradient.R:302-311`, `:332-340`,
  `:361-370`) via `withCallingHandlers` + `invokeRestart("muffleWarning")`
  keeping the first result. Output-neutral because `cor.test` already takes the
  asymptotic branch when ties exist, but the saving is smaller than it looks —
  `tryCatch(warning=)` already aborts the first call at the warning signal.

---

## Phase 11 — Bug fixes (one commit each)

Each documents how output changes and updates the affected baseline in the same
commit. Ordered by severity:

1. `R/import_xenium.R:19` defines a function named **`import_visium`**, shadowing
   the real one (bodies byte-identical, which hides it). `read_xenium_outs` calls
   `import_xenium_mex`, which **does not exist**, so Xenium MEX import always
   fails; check is silent only because `R/globals.R:5` lists the name. The copied
   body also reads Visium-format coordinates, so renaming is necessary but not
   sufficient.
2. `R/STclust.R:79` — `verbose = 1L` clobbers the user's argument.
3. `R/filter_data.R:170` — malformed 3-subscript index on a 2-D tibble, and uses
   `rownames()` (which is `"1","2",…`) instead of `[['gene']]`.
4. `R/detect_input.R:232` — `grepl("//.xlsx$", ...)` should be `"\\.xlsx$"`.
5. `R/STList.R:294`/`:302` — `mclapply` without `mc.cores`;
   `R/transform_data.R:101` hard-overwrites `cores = 1`.
6. `R/import_visium.R:112` / `R/import_xenium.R:108` — warning path references
   undefined `h5_file`, so the warning itself errors.
7. `R/pseudobulk_samples.R:361` — `mean/sd` sorted decreasing selects the *least*
   variable genes, opposite to the documentation.
8. `R/STdiff_volcano.R:24` — `is.numeric(x)` should be `is.numeric(samples)`.
9. `R/STList.R:948-950` — `arrange(!!array_col)` splices a character literal, so
   the sort is a silent no-op (needs `!!rlang::sym()`).
10. `R/STList.R:148-155` — `cores=NULL` with `rnacounts=NULL` errors "argument is
    of length zero".
11. `vignettes/basic_functions_vignette.Rmd:94` — undefined `thrane_tmp` (should
    be `visium_tmp`) throws inside the `tryCatch`, setting `net_check <- FALSE`, so
    **all 29 chunks of the only shipped vignette are currently skipped**.
12. The `cores=1`-vs-forked error asymmetry from fidelity rule 2, and the
    `_R_CHECK_LIMIT_CORES_` clamp from rule 12.

---

## Repository layout

**Deviation from your instruction, flagged for veto.** You asked for the plan and
`CLAUDE.md` in a git-ignored `.claude/`, but also said the layout should let
future clones continue development easily. Those conflict — an ignored `.claude/`
is lost on clone. Resolution:

- **`.claude/`** (gitignored + Rbuildignored, as instructed) —
  `settings.local.json` plus a short `CLAUDE.md` on how to work in this repo,
  which `@`-imports the tracked architecture doc so there's one source of truth
  and no drift.
- **`dev/`** (git-tracked, Rbuildignored) — `ARCHITECTURE.md` (the `STlist` slot
  contract, the `misc` key census, call graph, parallelism model, naming
  conventions, profiling switches, regression workflow, and the fidelity rules
  above), plus `environment/`, `baseline/`, `profiling/`. This is what makes a
  fresh clone productive.
- **`plans/`** (git-tracked, Rbuildignored) — matches your `spatialTIME`
  convention; this plan lives here.

**Do not run `devtools::document()`.** Installed roxygen2 is 8.1.0 against
`RoxygenNote: 7.3.2`; regenerating would churn all of `man/` and `NAMESPACE` and
bury the real diff. Every function added here is internal, so use `pkg::fn()` in
new code and no `@importFrom` is needed. If an `.Rd` must change, hand-edit it and
say so in the commit.

---

## Verification

Per commit:

```bash
Rscript -e 'devtools::load_all(); testthat::test_local()'
```

Per milestone:

1. `Rscript dev/baseline/compare_to_baseline.R` — every exported function over the
   fixture matrix, exact comparison against the Phase 3 goldens, at **`cores=1`
   and `cores=3`** (not 1 vs 2 — `mclapply` short-circuits below 2 elements, and
   3 samples at 3 cores is the smallest non-degenerate chunking).
2. `R CMD build` + `R CMD check --as-cran` — no growth over the Phase 1 reference
   NOTE/WARNING set, with `git diff --stat man NAMESPACE` empty.

Proof the profiling goal is met, not just claimed:

3. `options(spatialGE.cores = 1); profvis::profvis(STclust(fixture))` yields a
   profile whose stack contains `stclust_one_sample`, `calculate_dist_matrices`
   and `get_hier_clusters_dtc` as named frames — where the same run today shows
   only `mclapply` → `FUN`.
4. `dev/profiling/bench_units.R` times one sample's clustering without running
   `STclust`, which is impossible today.
5. `dev/profiling/coverage.R` reports fixture branch coverage, so gaps are known
   rather than assumed.
