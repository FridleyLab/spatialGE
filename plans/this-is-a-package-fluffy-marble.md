# spatialGE 2.0.0: scale out, interoperate, and ship to CRAN

## Context

The Phase 1–11 refactor made the package profilable and left every unit of
parallel work as a named, top-level, directly callable function. That work
deliberately preserved behavior. This plan spends that seam: it changes how data
moves through those units so the package works on datasets and platforms it
currently cannot, then adds the interop and method surface, then submits to CRAN.

Five structural facts drive the phase order:

1. **Workers receive the whole STlist, not their slice.** 23 of the 24
   `sge_lapply` call sites iterate an integer/name index and reach into the
   enclosing frame, so each worker's environment closes over every sample's
   `@counts` and `@tr_counts`. `R/STenrich.R:159` (`x = NULL  # Clear memory by
   offloading STlist`) is the only site that drops it. Under `fork()` this costs
   copy-on-write faults as the GC touches mark bits; under any non-forking
   backend it would serialize the entire object *per task*. The `*_one_sample`
   units already take explicit per-sample arguments, so this is mechanical.

2. **There is no parallelism on Windows at all.** `resolve_cores()`
   (`R/utils_parallel.R:39`) returns `1L` on Windows *before* consulting the
   override, and `sge_lapply` forwards only to `parallel::mclapply`. A CRAN
   release ships a serial package to most of its users.

3. **`SThet()`'s default path builds a complete graph.**
   `create_listw_from_dist()` (`R/utils.R:210-247`) materializes
   `as.matrix(dist(coords))`, then for every spot allocates two length-(n−1)
   vectors *with character names*, computing membership via
   `rownames(adj) %in% rownames(adj)[id]` — an O(n) string comparison inside an
   O(n) loop — producing a `listw` in which every spot neighbors every other. At
   4,000 Visium spots that is ~16M weights plus 16M names; at Xenium/CosMx scale
   it is infeasible. It has zero test coverage and zero benchmarks.
   `create_listw_from_knn()` (`R/utils.R:192`) hardcodes `cores = 1`.

4. **The numeric safety net is invisible outside this working tree.**
   `tests/golden/` is gitignored and `SPATIALGE_GOLDEN` is unset in
   `.github/workflows/test.yml`, so all three golden tests skip on CI and on any
   fresh clone. No direct numeric test exists for `SThet`, `STclust`,
   `STgradient`, `STenrich`, `STdiff`, or `gene_interpolation`. Line coverage is
   59.41%. The fixture is 3 samples × ~200 genes, which is why fact 3 was never
   observed — `dev/profiling/FINDINGS.md` says so itself, and its every
   `file:line` reference is now stale (all six of its Phase 10 targets shipped).

5. **The accessor seam is 3% adopted.** `R/utils_accessors.R` is 17 lines and 4
   one-line functions serving 9 call sites, against ~306 direct `@slot` reads and
   44 unabstracted writes. Worse, the accessors return *eagerly* and the
   `*_one_sample` signatures then take a bare matrix — so provenance is erased
   before any operation happens. Routing all reads through the accessors as they
   exist would not create a lazy seam.

**Intended outcome.** 2.0.0 goes to CRAN with: workers that receive only their
own sample; parallelism on all three platforms; an opt-in parquet-backed store
for data that exceeds memory; the O(n²) `listw` construction replaced;
bidirectional Seurat v5 interop; generics that dispatch on STlist vs Seurat; and
the functional cell–cell interaction method integrated.

## Standing constraints

**Commit, never push. This is a hard limit with no exceptions.** Every phase
below lands as one or more commits on `alex-refactor`. Nothing is ever pushed to
any remote — not `origin`, not a fork, not a new branch — regardless of how the
work is going or how convenient it would be. `git push` is not part of any step
in this plan. Publishing is yours alone to do, and that includes the CRAN
submission itself (Phase 9 prepares the tarball and gets `R CMD check --as-cran`
clean; it does not submit).

The rest, unchanged from the refactor:

- Commit when a meaningful, self-contained change is complete.
- No mention of AI tooling in any commit message, comment, or file.
- Do not run `devtools::document()` — installed roxygen2 is newer than
  `RoxygenNote`, and regenerating would churn all of `man/` and `NAMESPACE` and
  bury the real diff. Hand-edit `.Rd` and say so in the commit message.
- Any output-affecting change is its own commit that both makes the change and
  regenerates `tests/golden/` via `dev/baseline/capture_baseline.R`, so `git log`
  always has a single commit to bisect to for "when did this number change".

---

## Phase 0 — Make the verification real (blocking; no `R/` changes)

Nothing after this is verifiable without it, and Phase 4 is a numeric rewrite.

- **Commit the golden baseline.** Un-gitignore `tests/golden/`, keep it
  Rbuildignored. Set `SPATIALGE_GOLDEN=true` in `.github/workflows/test.yml`. The
  54 `.rds` files are the only numeric net that exists; a net that skips by
  default is not a net.
- **Add scale fixtures** to `inst/testdata/make_fixtures.R`: one Visium-like
  (~4,000 spots × ~2,000 genes) and one Xenium-like (~50,000 cells × ~400 genes),
  generated from a seed, materialized to `tempdir()`, and **not** part of the
  golden capture (too slow). They exist so `dev/profiling/bench_units.R` can see
  asymptotics. Fact 3 is invisible below ~1,000 spots.
- **Add direct numeric tests** for the six compute algorithms —
  `tests/testthat/test-numeric-<fn>.R` — asserting against literal expected
  values on the small seeded fixture at a stated tolerance. Prefer the
  order-insensitive variants per fidelity rule 13(b): `SThet()` without `k`,
  `STgradient(distsumm='min')`.
- **Extend `dev/baseline/compare_baseline.R`.** `bl_compare_stlist()` compares 2
  of the 10 `@misc` keys; add `scaled_pbulk_mtx` and `sp_images`, and add a
  `HeatmapList` branch reading `@ht_list[[1]]@matrix` for `pseudobulk_heatmap()`.
  Both gaps are documented in `dev/ARCHITECTURE.md`.
- **Refresh the stale docs in the same commit:** rewrite
  `dev/profiling/FINDINGS.md` against current line numbers, and fix
  `dev/ARCHITECTURE.md` fidelity rule 3 (all five `stlist_readers.R` reader sites
  now pass `cores=cores`) and its reference to `sge_worker_msg` (no such function
  exists; the 18 `system(sprintf('echo "%s"', msg))` calls are inline).

---

## Phase 1 — Turn the accessors into a real seam (output-neutral; one commit per slot family)

Every later phase depends on this. Widen `R/utils_accessors.R`:

```r
sample_counts(x, sample, genes=NULL, cells=NULL, as=c('dgCMatrix','arrow','manifest'))
sample_counts_lazy(x, sample)          # an arrow query or an in-memory shim
sample_counts_meta(x, sample)          # dim/dimnames only, zero I/O
`sample_counts<-`(x, sample, value)    # the 44 writes get a setter
```

Same four forms for `tr_counts`, `spatial_meta`, `gene_meta`. Add the missing
`sample_meta`/`misc`/`gene_krige` accessors. In-memory mode makes the lazy form a
trivial shim, so this phase changes no behavior and is fully golden-covered.

Then convert reads by category — the classification is complete and lands in
`dev/ARCHITECTURE.md` as a table:

| Category | Count | Conversion |
|---|---|---|
| (D) structural/metadata only | 76 | `sample_*_meta()`; zero I/O |
| (A) arrow-pushable aggregate | 24 | `group_by |> summarise |> collect` |
| (B) gene-subset materialization | 15 | `sample_*(genes=)` |
| (C) whole-matrix materialization | 16 | explicit `materialize()`; unavoidable |

Sequence the commits (D) → (B) → (A) → (C), smallest blast radius first.

### As executed — three scope corrections

**(D) is 22 sites, not 76.** 34 of the 76 are `names(x@counts)` /
`names(x@tr_counts)`, which need no change: the slot stays a named list in user
order with one manifest per sample, so `names()` keeps working while
`nrow(x@counts[[i]])` does not. Converting them would have been churn that buried
the 22 that actually break. Done in `3279e6c`.

**Setters dropped.** `sample_counts(x, s) <- v` expands to
``x <- `sample_counts<-`(x, s, v)``. R mutates in place only when an object's
reference count is 1, and routing a slot write through an extra function frame
risks defeating that — turning a free mutation into a full copy of a sample's
matrix. That is 44 sites of churn for a real performance risk, and Phase 5 does
not need it: writes only require interception in `filter_data` and
`transform_data`, which can call the store writer directly. A uniform setter
layer would solve a 2-site problem in 44 places.

**The `resolve_samples()` rewrite is unnecessary.** The original justification
was that Phase 7 could not dispatch through `methods::slot(x, slot)`. That holds
only for genuine dual implementations, which Phase 7 rejects in favour of
coercing at the boundary — under that design every internal only ever sees an
STlist, so the string-keyed slot lookup is never reached with a Seurat object.
It also survives Phase 5 untouched, since `names(x@counts)` works on a list of
manifests. Left alone, with its 20 tests pinning the deliberately-divergent
behaviour of its 14 callers.

### Also landed here (measured, not planned)

Densification removals found while converting the (A) sites. Each is its own
commit with a before/after measurement:

| change | scale | speedup | output |
|---|---|---|---|
| `filter_data` full densification removed | 200×50k | 27× | neutral (integer sums exact) |
| `apply(m, 1, sd)` → `sparse_row_sds()` | 200×50k | 42× | moved `gene_stdevs` at 1e-16 |
| log transform sparse end to end | 200×50k | 33× | bit-identical |

---

## Phase 2 — Give each worker only its own slice (output-neutral; one commit per function family)

At each of the 23 `sge_lapply` sites, build the payload *outside* the call and
pass it as `X`, so the closure captures nothing large. `R/STList.R:285,293` and
the five `stlist_readers.R` sites already do this and are the template:

```r
payloads = lapply(samples, function(s) list(name=s,
                                            counts=sample_counts(x, s),
                                            meta=sample_spatial_meta(x, s)))
res = sge_lapply(payloads, function(p) worker(p$counts, p$meta, ...), cores=cores)
```

Three sites need more than mechanics:

- `R/STclust.R:130` also captures `trcounts_df`, a list of *every* sample's
  top-gene matrix. Slice it per payload.
- `R/STdiff_fit.R:171` captures an already-densified `expr_df` plus
  `non_sp_models` for all samples.
- `R/SThet.R:160,204` capture `@misc$sthet$listws` — every sample's `listw`.
  After Phase 4 these are small; do Phase 4 first if ordering is convenient.

Preserve exactly: each site's order source (`names(x@spatial_meta)` vs
`names(x@counts)` vs `unique(combo[[1]])` vs the user's `samples=`) per fidelity
rule 7; the per-site `preschedule` value per rule 4; and the three `<<-` sites
(`R/STclust.R:145-148`, `R/STdiff_fit.R:364-367`) which must stay in their
enclosing frames per rule 6. `STenrich()` re-derives `cores` from three different
unit counts in one call — do not collapse it (rule 5).

### As executed — the stated deliverable was not achievable, and why

Measured with `dev/profiling/rss_peak.R`, `transform_data` over 4 samples,
before and after slicing `transform_data`'s three dispatch sites:

| preset | fork cost before | fork cost after |
|---|---|---|
| visium | +1619.5 MB | +1618.8 MB |
| xenium | +5022.8 MB | +5098.5 MB |

**Payload slicing does not reduce fork memory, and cannot.** Under `fork()` the
child inherits the parent's entire address space regardless of what the closure
captures; closure contents govern *serialization*, not forking. So the original
framing of this phase — "cuts COW pressure" — was wrong.

Worse, the metric was wrong too. Summed RSS charges copy-on-write shared pages to
every process, so the ~1.35 GB of loaded packages is counted once per worker:
that is the entire source of the apparent "+5 GB for 4 workers". Measuring the
largest single process instead:

| preset | cores | tree_sum (misleading) | worst single process |
|---|---|---|---|
| visium | 1 | 299.1 MB | 1651.9 MB |
| visium | 4 | 1950.3 MB | 1661.7 MB |
| xenium | 1 | 384.5 MB | 1737.3 MB |
| xenium | 4 | 5396.8 MB | 1769.9 MB |

No process grows meaningfully with core count. **There is no per-worker
environment duplication under fork on this workload** — the effect is real in
what a memory monitor displays, but the mechanism is shared-page accounting, not
copying.

Three consequences, which reshape the rest of the plan:

1. The genuine lever on memory is per-worker *allocation*, not closure contents.
   That is what the Phase 1 densification removals addressed (27x, 42x, 33x), and
   it is where further wins live.
2. Phase 2 is retained purely as the **prerequisite for Phase 3**. Its benefit is
   deferred, not present. Under a non-forking backend, `X` elements and the
   closure environment are both serialized per task, so an un-sliced call site
   would ship the entire STlist per worker — turning today's free fork into a real
   copy. Slicing first is what makes Phase 3 viable instead of a regression.
3. For Windows, where there is no fork at all, this is not optional: a
   non-forking backend is the only way to get parallelism, and slicing is what
   makes it affordable.

---

## Phase 3 — `future` backend, so Windows gets cores (**output-affecting**; own commit + baseline regen)

Add `future` and `future.apply` to `Imports`. Rewrite `sge_lapply` to dispatch on
a plan while keeping the signature:

- Default on Unix: `future::plan(multicore)` — preserves today's fork behavior.
- Default on Windows: `future::plan(multisession)` — PSOCK. **Delete the
  `.Platform$OS.type == 'windows'` short-circuit at `R/utils_parallel.R:39`**, and
  note that it currently sits *before* the override, so `SPATIALGE_CORES` is
  ignored on Windows today.
- `sge_cores_override()` maps to `plan(sequential)`, keeping `profvis` able to see
  through to the named units.
- Respect a user's own `plan()` if they set one; document that as the knob.

This changes output, deliberately and for the better: L'Ecuyer streams via
`future.seed` make `transform_data(method='sct')` reproducible at `cores>=2`,
which fidelity rule 1 currently documents as impossible. Amend rule 1 rather than
preserving it. `pseudobulk_samples(calc_umap=TRUE)` stays unseeded via
`uwot::umap`'s own threads — out of scope, note it.

Error granularity changes too (rule 2). Keep the `class(res_mod) == 'HLfit'` guard
at `R/STdiff_fit.R:381` and copy that pattern to every new failure point. Keep
`sge_bpparam()` for the one `GSVA::gsva(BPPARAM=)` consumer.

Verification is the golden gate at `cores=1` and `cores=3` under each plan, plus
the `cores1`-vs-`cores3` invariance test — which should now pass for `sct`, where
it previously could not.

---

## Phase 4 — Replace the O(n²) hot spots (mixed; one commit each)

Ordered by measured impact on the Phase 0 scale fixtures. Tiers 1a–1b are the
scaling blockers; the rest are allocation wins.

**1a. `create_listw_from_dist()` → C++ (`R/utils.R:210`).** The complete graph is
not a performance bug, it is a design bug: a distance-weighted `listw` where
every spot neighbors every other means `SThet()`'s default is O(n²) in both time
and memory with a character-names constant factor. Replace with a compiled
Moran/Geary that takes the coordinate matrix and the expression matrix and
returns all genes at once, never materializing the graph. This **changes
`SThet()` output** if the weight scheme changes — so first establish whether an
exact-equivalent formulation exists (same normalization, `1/d`, `style='B'`, no
truncation) and prove bit-equivalence against `main` on the small fixture. If it
does not, that is a documented, NEWS-worthy behavior change and needs its own
commit plus baseline regen. `sthet_stat_one_sample()` (`R/SThet.R:245-269`) loops
genes in R and calls `spdep::moran`/`geary` per gene against that graph; fold the
loop into C++. **Confirm first whether `spdep::moran`/`geary` are compiled** — it
determines whether the win is the graph construction alone or both halves.

**1b. `stgradient_dists_one_sample()` → C++ (`R/STgradient.R:256-327`).**
Materializes the full n×n `as.matrix(dist(coords))`, then wraps an n(n−1)/2 vector
in a `data.frame` to take a `min`, then `apply(..., 1, min|mean)` over a
data.frame. A C++ pass computing each non-reference spot's min/avg distance to
the reference set never allocates n×n. A fixed summation order also removes the
`distsumm='avg'` tie sensitivity that `dev/xver/REPORT.md` measures at 4.3e-3 —
that is an output change, so it is its own commit.

**2. `non_spatial_de()` / `stdiff_mean_test()` (`R/STdiff_fit.R:325-397,480-526`).**
**PARTIALLY DEFERRED** (see `dev/profiling/FINDINGS.md`, "Other hot spots
observed"). `non_spatial_de()` is a `spaMM::fitme` mixed-model fit per gene x
cluster, not trivial arithmetic; its cost is spaMM's own RcppEigen internals, so
there is no local vectorization lever and no bit-exact rewrite short of a
research task. The actionable part is `stdiff_mean_test()` (the non-default
`t_test`/`wilcoxon` path): `for(i in 1:nrow(combo))` with 4–5 separate dplyr
passes per combination, vectorizable over the whole gene × cluster matrix and
possibly achievable in pure R. The 753 ms / 8 GCs figure predates the
caller-side `expandSparse()` subsetting (`74e5935`); re-measure before investing.

**3. `calculate_dist_matrices()` / `calculate_weighted_dist()` (`R/STclust.R:258-300`).**
**MEASURED, NOT IMPLEMENTED** (see `dev/profiling/FINDINGS.md` and the note at
`calculate_weighted_dist()`). `Reduce('+', Map('*', dmxs, weight_ls))`
allocates three n x n matrices per weight (two scaled copies plus the sum) --
measured 384 MB per weight at 4,096 spots. Two corrections to the plan: the
"lower-triangle vector per weight" formulation is not possible, because
`get_hier_clusters_dtc()` passes the full matrix to
`dynamicTreeCut::cutreeDynamic(distM=)`; and a C++ elementwise fusion, while
3x lighter on allocation, is not bit-identical (clang contracts `wa*a + wb*b`
into an FMA, 1 ULP off) and was slower in the dev build. Revisit only with
`-ffp-contract=off` and an optimized build.

**4. `stenrich_score_one_combo()` (`R/stenrich_helpers.R:30`).**
`data.frame(as.list(pw_avg_exp))` builds one column *per cell* — a 100,000-column
data.frame per combination on a CosMx sample. Keep it a named numeric vector.
Trivial R fix, large real-data effect. `check.names=FALSE` is load-bearing here
(rule 9).

**5. `expandSparse()` (`R/sparse_utils.R:39-42`)** and its callers, which each do
`t(expr_df) %>% as.data.frame()` — three full dense copies. `Matrix::t()` before a
single `as.matrix()` halves it. Keep `expandSparse` as the one densification
chokepoint; callers depend on it returning a `data.frame` (see Phase 5 note 9.7).

**Explicitly not Rcpp targets:** `spatial_de()`/`stdiff_spatial_one_cluster()`
(20.8 s, 27.4 GB — that is `spaMM`'s own RcppEigen internals; the fix is chunking
guidance and documentation), `sctransform::vst`, `gstat` kriging,
`Seurat_FindVariableFeatures` (already RcppEigen), and `computeSubsampleSums`
(already Rcpp; any change to its sampling moves p-values per rule 13 row 1).
Also: `LinkingTo` declares `RcppProgress` and `src/seurat_helpers.cpp:24` declares
the dependency, but no `Progress` object is ever instantiated — either use it for
the new long-running kernels or drop it.

---

## Phase 5 — Opt-in parquet COO store (new API surface; several commits)

**Design.** Parquet is the store; `dgCMatrix` is a materialized view. One parquet
file per sample per slot, schema `gene_idx:int32, cell_idx:int32, value:double`,
**sorted by `(gene_idx, cell_idx)`**. Tabular slots and results tables are
parquet too. `arrow` is already in `Imports` and currently unused in `R/`, so the
dependency cost is paid.

Two non-obvious constraints on that sort key. Sorting by gene gives row-group
min/max statistics, so `filter(gene_idx %in% ...) |> collect()` prunes at the
row-group level — this is the `filter |> collect` win, and it avoids
hive-partitioning by gene, which would create tens of thousands of small files.
Including `cell_idx` in the key is what makes `group_by(gene_idx) |>
summarise(sum(value))` match `Matrix::rowMeans`' column-ascending walk
bit-for-bit; sorting by gene alone leaves intra-gene order unspecified and drifts
the last bits of `@gene_meta`, which is golden-covered.

**Layout.**

```
study.stlist/
  manifest.json                 # schema_ver, root, sample_order, portable
  sample_meta.parquet
  samples/<store_id>/
    counts.parquet   tr_counts.parquet
    spatial_meta.parquet   gene_meta.parquet
  results/{sthet,stgradient,stdiff,stenrich}.parquet
```

**Per-sample manifest** — held in the slot element so `sample_counts()` is the
only seam: `sample`, `store_id`, `path`, `generation`, `genes` (chr, authoritative
row order), `cells` (chr, authoritative column order), `n_genes`, `n_cells`,
`nnz`, `gene_dup`, `cell_dup`, `value_type`, `sort_key`, `explicit_zeros`,
`rowgroup_bounds`, `fingerprint`, `materialized`. Object-level in
`@misc$store`: `root`, `sample_order`, `schema_ver`, `portable`. Add the
`@misc$store` row to the key census in `dev/ARCHITECTURE.md` in the same commit.

**Never store an `arrow::Dataset`/`Table`/`RecordBatchReader` in the object** —
they are R6 wrappers around external pointers that survive neither `saveRDS` nor
`fork()`. Open lazily, per process, from `path`.

`validate(x)` reads parquet footers only, no data pages: paths exist and
fingerprints match; `identical(names(x@counts), store$sample_order)` (rule 7);
schema is exactly the three expected columns; `rowgroup_bounds` is non-decreasing
in `min_gene_idx`, i.e. the file really is gene-sorted (otherwise pushdown
silently degrades *and* `rowMeans` silently loses bit-fidelity);
`length(genes) == n_genes`; cross-slot `tr_counts$genes ⊆ counts$genes` in the
same relative order, which `R/filter_data_stages.R:216-217` assumes; and
`setequal(counts$cells, spatial_meta$libname)`. `validate(deep=TRUE)` checks index
bounds and duplicate `(gene_idx, cell_idx)` pairs.

**Order of commits.** Writer + manifest + `validate()` → materializer → the (A)
aggregate rewrites → `filter_data` → `transform_data`. Each is golden-gated in
memory mode, then re-run in backed mode.

### Eleven things that will otherwise break this

1. **`fork()` + arrow is a deadlock class, not a hypothesis.** Every site that
   benefits reads its slot *inside* `mclapply`. Arrow's C++ runtime holds a
   thread pool; calling into it from a child that inherited that state hangs.
   Mitigations: never let a handle cross the boundary (capture `path` + manifest,
   open in the child); `arrow::set_cpu_count(1)`. **This is why Phase 3 precedes
   Phase 5** — `multisession` sidesteps it entirely.
2. **`saveRDS` stops round-tripping, and that takes out the golden gate.**
   `dev/baseline/capture_baseline.R:22` saves full STlists to `tests/golden/` and
   `test-golden.R` reads them in a *later process*; a store under `tempdir()` is
   gone by then. `bl_compare_stlist()`'s `cmp_matrix_list`
   (`dev/baseline/compare_baseline.R:112`) branches on `is(la[[n]], 'dgCMatrix')`
   and would fall through to raw `waldo` on external pointers — the same noise
   that let the inverted-CV bug survive two releases. Fix in the same commit:
   `materialize()` both sides before comparing, and force in-memory mode in
   `bl_pipeline()`. This also hits any user who `saveRDS`es an STlist, which is
   the normal workflow after a multi-hour `STdiff`.
3. **Sample names are not valid filenames.** Rule 11 already documents `"`, `$`,
   backticks, `\` in real names; paths add `/`, `:`, leading `.`, macOS NFC/NFD,
   the 255-byte limit, and APFS case-insensitive collision (`Sample1` vs
   `sample1` are distinct samples, one file). Use
   `store_id = sprintf('%04d-%s', k, substr(hash(sample), 1, 12))`; keep the human
   name only in the manifest.
4. **Copy-on-modify.** `y <- filter_data(x, ...)` must leave `x` intact. An
   in-place parquet rewrite silently corrupts the caller. Every write creates a
   new `generation`. But `filter_data()` runs five stages in sequence
   (`R/filter_data.R:97-116`), so a naive implementation writes five full
   generations per call and never collects them. Fuse stages 1–4 (pure mask
   composition, so nothing observable changes) and write once — noting stage 5
   computes its statistics on stages 1–4's output, and stages 2–4 touch only a
   `samples` subset, so different samples can legitimately sit at different
   generations. Per-sample-per-slot manifests accommodate that.
5. **All-zero rows and columns do not exist in COO.** `filter_data()`'s defaults
   are all `>= 0` (`R/filter_data.R:80-85`), so all-zero genes and cells are
   retained and are real. Every per-cell aggregate must right-join
   `manifest$cells` filling `0`, or `min(colSums)` and the `spot_mask` at
   `R/filter_data_stages.R:169-171` are wrong. The manifest's `genes`/`cells` are
   the *only* authority for dim and order — never reconstruct them from the data
   (rule 13(c)).
6. **`sd`/`var` over a COO store are over nonzeros only** — semantically wrong,
   not merely imprecise. `R/filter_data_stages.R:223` and
   `R/transform_data.R:283` (`apply(m, 1, sd)`, which silently densifies today)
   need a two-pass mean-then-centered-sum with the structural zeros' contribution
   `mu^2 * (nc - nnz_g)` added back. A one-pass `sqrt((ss - s^2/n)/(n-1))` can
   produce `NaN` from a negative radicand on near-constant genes. `gene_stdevs`
   is golden-covered: budget one deliberately output-changing commit.
7. **Rule 9 arrives from the opposite direction.** `collect()` returns row-group
   order, which under a gene sort is *lexicographic* — exactly the `m[genes, ]`
   reordering rule 9 forbids. Every materializer must go through
   `match(gene_idx, sort(gidx))` against manifest order. Separately,
   `expandSparse`'s `data.frame()` currently *errors* on duplicate gene names
   (rule 9's second clause); a `gene_idx`-keyed COO makes duplicates
   representable and the error vanishes, so `validate()` must raise on
   `gene_dup` to preserve it.
8. **Do not use arrow's `grepl`.** `R/filter_data_stages.R:145` resolves
   `spot_pct_expr` (default `"^MT-"`, POSIX ERE) against rownames; arrow uses RE2.
   Resolve the regex in R against `manifest$genes`, then push down
   `gene_idx %in% gidx`.
9. **Arrow has its own CPU thread pool defaulting to `detectCores()`**, which
   R CMD check will see under `_R_CHECK_LIMIT_CORES_` (rule 12). Clamp it in
   `R/zzz.R` alongside the existing override, and hold it identical across
   `capture_baseline.R` and `bl_compare()` or diffs are spurious. Also set
   `options(arrow.use_threads=FALSE)` for any aggregate whose result is
   golden-covered, since arrow may otherwise chunk a `group_by`.
10. **Two paths cannot be converted and must require explicit `materialize()`.**
    `plot_counts()` (`R/plot_counts.R:83,94`) plots zero-inflation — the zeros
    *are* the data — and `sample(1:nrow(df_tmp), ...)` draws over positions in a
    frame whose `nrow` is `ng × nc` dense, which rule 13 row 4 names as
    position-dependent. `GSVA::gsvaParam(as.array(...))`
    (`R/stenrich_helpers.R:147`) needs a full dense `ng × nc` base array on the
    `score_type='gsva'` path; strictly worse than `expandSparse` and unavoidable.
    Document both.
11. **Type contracts.** `R/stlist_process.R:63`'s deliberate
    `as.double(as.vector(Matrix::colSums(m != 0)))` exists because `colSums` on an
    `lgCMatrix` returns integer; arrow's `n()` returns int64, so coerce or rule 8
    breaks. `nnz` can exceed 2^31 on a merged Xenium store — keep it `double`, and
    have `validate()` warn that `Matrix::sparseMatrix` itself caps at 2^31−1
    nonzeros before someone hits it three hours in.

Also repair the five hard `dgCMatrix` class assertions
(`R/sparse_utils.R:15`, `R/stlist_process.R:38`, `R/transform_data.R:84` — which
uses `class(...)[1]`, so any S4 class with a superclass chain fails to match —
`R/seurat_helpers.R:32,35`, and the `is(x, 'STlist')` dual-dispatch at
`R/utils.R:371`). Note `R/seurat_helpers.R:36` calls **`as.sparse()`, which does
not exist anywhere** and is masked only by `utils::globalVariables` in
`R/globals.R:1`; today neither branch fires because every caller passes a
`dgCMatrix`. It is a landmine sitting on the (C) path — fix it in Phase 1.

---

## Phase 6 — Seurat v5, both directions (new feature; two commits)

`read_seurat()` (`R/stlist_readers.R:44-66`) is the whole of the current interop
and is v4-only in three ways: `rnacounts@assays$Spatial@counts` (v5 uses `Assay5`
with `@layers[['counts']]`), a hardcoded `"Spatial"` assay name, and
`as.data.frame()` on the full matrix — densifying every gene × every cell across
all slices before any per-slice subsetting. `R/STList.R:212` then hardcodes
`platform = 'visium'` with a comment admitting it cannot.

- **Inbound:** support v4 and v5 by branching on `SeuratObject::Version()`; take
  the assay name as an argument defaulting to `DefaultAssay()`; go
  sparse-to-sparse via `SeuratObject::LayerData()`, never `as.data.frame`; infer
  `platform` from the image class (`VisiumV1`/`VisiumV2`/`FOV`) instead of
  assuming. `dev/xver/REPORT.md` records the Seurat path as never exercised, so
  add fixtures for both versions.
- **Outbound:** `as_Seurat(x, assay='Spatial')` — an `Assay5` per sample or one
  merged object with a slice per sample, carrying `spatial_meta` to
  `@meta.data`, `gene_meta` to `@assays[[a]]@meta.data`, images from
  `@misc$sp_images`, and STclust labels as identities.

Promote `SeuratObject` from `Suggests` to `Imports` only if outbound is
unconditional; otherwise keep it in `Suggests` behind `requireNamespace()`.

---

## Phase 7 — Generics (needs Phase 1; one commit)

All 24 exported functions are plain closures and all 23 that take an STlist name
it `x`, so `setGeneric(..., useAsDefault=)` is mechanical. There are zero
`setGeneric`, zero `setValidity`, and three `setMethod`s today (`show`,
`summary` — a byte-for-byte duplicate of `show`'s body — and `dim`).

**Recommendation: dispatch, then coerce at the boundary.** The Seurat method
converts to an STlist, calls the STlist method, and writes results back. Do not
write dual implementations: the bodies read slots directly at ~306 sites, so a
genuine Seurat method would mean maintaining two of everything. Coercion is
honest about the cost, is one code path to test, and can be replaced per-function
later where it pays. Document the round-trip cost in each method's `.Rd`.

Prerequisites: `assert_stlist()` (`R/utils_validate.R:15`) becomes a
per-method precondition rather than a universal gate, and `resolve_samples()`
must already be `sample_order`-based from Phase 1. Add `exportClasses(STlist)`
and `exportMethods()` to `NAMESPACE`, add the `setValidity` that has never
existed, and de-duplicate `show`/`summary`.

---

## Phase 8 — Functional cell–cell interaction from multiScaleAC (new feature; two commits)

`multiScaleAC` (repo `ACSoupir/KernelCommunication`, branch
`feature/s4-classes-modular-design`) computes bivariate Moran's I as a function of
Gaussian kernel bandwidth, giving a curve `I(h)` per gene pair, with RcppEigen +
optional OpenMP, closed-form null moments, and an `export_mxfda()` bridge for
FPCA.

- **Vendor, do not depend.** CRAN forbids a hard dependency on a GitHub-only
  package. Copy the needed `src/*.cpp` and `R/*.R` into spatialGE with a header
  comment on every file recording source repo, branch, commit SHA, retrieval
  date, and upstream license (GPL >= 2 — **confirm compatibility with
  spatialGE's MIT before copying**; this may force a license decision and is the
  one open item in this phase). Mirror the existing precedent in
  `inst/COPYRIGHT` and the `cph` entry in `DESCRIPTION` used for Seurat's
  `SparseRowVar2`. When multiScaleAC reaches CRAN, move it to `Imports`, delete
  the vendored files, and make the wrappers call `multiScaleAC::`.
- **Scope the vendored code to the new feature only.** Per your direction: the
  kernel is used *only* for the functional cell–cell interaction calculation and
  **not** for spatialGE's univariate Moran's I. If the Phase 4 work ever wants to
  reuse it, that requires first demonstrating bit-identical output against the
  installed `main` branch (1.2.2) using the `dev/xver/` harness — which already
  compares an installed `main` against `alex-refactor` across 17 scenarios. Until
  that benchmark exists, the two code paths stay separate.
- Expose `STinteract(x, pairs, bandwidths, ...)` returning a tidy curve table
  per sample and gene pair, an `STplot`-family plotting function for `I(h)`, and
  an `mxfda` bridge in `Suggests`. Carry through upstream's two documented
  traps: bandwidths are in coordinate units (not micrometres), and `pnorm()` on
  the `z_scores` slot is ~100× anti-conservative in the tail.

---

## Phase 9 — CRAN preparation (final; one commit — preparation only, no submission)

`R CMD check --as-cran` clean on macOS, Linux, and Windows (win-builder devel +
release). This phase ends with a checked tarball and a green NEWS; **submitting it
is yours to do**, per the standing constraints. Specific known items: confirm `src/*.o` and `src/spatialGE.so` are
untracked and Rbuildignored; add `R_registerRoutines`/`.registration=TRUE`, which
does not exist today; keep every example inside the 5-second budget and behind
`\donttest{}` where it downloads; verify `_R_CHECK_LIMIT_CORES_` behavior for both
`count_cores()` and arrow's thread pool (rule 12, Phase 5 note 9); resolve the
license question from Phase 8 before submission. Write the full 2.0.0 NEWS —
`NEWS.md` currently documents the refactor but nothing about performance,
parallelism, Seurat, the store, or `STinteract`.

---

## New fidelity rules to add to `dev/ARCHITECTURE.md`

Written in the existing style — bolded claim, mechanism with `file:line`,
consequence — as rules 14–17:

14. **Never let an arrow handle cross a `fork()`.** Capture `path` + manifest and
    open inside the worker. Arrow's thread pool inherited across a fork hangs.
15. **The manifest is the only authority for dim, dimnames, and order.** COO
    omits all-zero rows and columns; reconstructing `genes`/`cells` from the data
    silently shrinks and reorders the matrix, breaking rule 13(c).
16. **The COO sort key is `(gene_idx, cell_idx)`, not `gene_idx`.** The second
    component is what makes arrow's `group_by |> sum` match `Matrix::rowMeans`'
    summation order bit-for-bit on golden-covered `@gene_meta` values.
17. **Every parquet write creates a new generation.** An in-place rewrite
    violates R's copy-on-modify contract and silently corrupts the caller's
    object.

---

## Verification

**Per commit:**

```sh
Rscript -e 'devtools::load_all(); testthat::test_local()'
```

**Per phase:**

```sh
SPATIALGE_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_local()'
Rscript dev/profiling/bench_units.R      # scale fixtures, diff against prior phase
```

**Proof each goal was met, not just claimed:**

1. **Memory.** Peak RSS on the 50,000-cell fixture at `cores=4`, before Phase 2
   and after, with per-worker RSS. The claim is sublinear growth in core count.
2. **Windows.** `bench_units.R` on a Windows runner in
   `.github/workflows/test.yml` showing wall-clock scaling with cores where
   today it is flat.
3. **Scaling.** `SThet()` wall-clock and peak RSS at 1k / 4k / 50k spots before
   and after Phase 4a. The claim is that the O(n²) term is gone.
4. **Numeric fidelity.** The golden gate green at `cores=1` and `cores=3` under
   both `multicore` and `multisession`, plus the invariance test — which should
   pass for `sct` after Phase 3, where rule 1 says it cannot today.
5. **Store fidelity.** `bl_pipeline()` run twice, once in memory mode and once
   backed, compared with `bl_compare()`. Every difference must be a documented,
   NEWS-listed intentional change.
6. **Order sensitivity.** `dev/xver/30_order_probe.R` re-run after Phases 4 and 5.
   A diff in `@counts`, `@tr_counts`, or `vst.*` is real (rule 13(c)); a diff in
   the four order-sensitive paths is not a regression until the probe says so.
7. **Cross-version.** The full `dev/xver/` harness, installed `main` (1.2.2)
   against the finished branch, with every differing leaf path classified.
   Classification fails closed. Note the harness is currently untracked — track
   it, or this proof is not reproducible.
8. **Seurat round-trip.** `STlist -> as_Seurat -> STlist` preserves counts,
   coordinates, and metadata exactly, for both v4 and v5 inputs.
9. **CRAN.** `R CMD check --as-cran` with zero ERRORs and zero WARNINGs on all
   three platforms.

---

## Sequencing rationale

```
P0 verification debt  ──┬─> P1 accessor seam ──┬─> P2 payload slicing ──> P3 future backend ──┐
                        │                      │                                              │
                        └─> P4 O(n²) hot spots │                          P5 parquet store <───┘
                                               └─> P7 generics <── P6 Seurat v5
                                                                                   P8 STinteract
                                                                                   P9 CRAN
```

P0 is blocking because P4 is a numeric rewrite and the current tree cannot verify
one. P1 gates everything that touches data. P3 precedes P5 because
`multisession` avoids the arrow-across-fork deadlock. P4 precedes P5 because the
complete-graph fix is storage-independent and is the largest single scaling win.
P6 and P8 are independent and can land any time after P0.

**Flagged for veto:** Phase 8's GPL-vs-MIT license question may force either a
license change on spatialGE or waiting for multiScaleAC on CRAN. Phase 4a may not
have a bit-equivalent formulation, in which case `SThet()`'s default output
changes. Phase 5's `sd` fix changes `@gene_meta`. Each is called out in place.
