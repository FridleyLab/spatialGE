# Plan: prove `alex-refactor` (2.0.0) matches `main` (1.2.2) except where it fixes it

## Context

`alex-refactor` is a 54-commit refactor bumping `spatialGE` 1.2.2 → 2.0.0. Its
stated contract (`plans/this-is-an-r-velvety-music.md`) is *bit-identical output
with an unchanged user API*, with genuine bugs fixed afterwards in separate,
individually-documented commits.

That contract has so far only ever been checked **against itself**:
`tests/golden/` was captured from this same tree before the `R/` edits began, and
`tests/testthat/test-golden.R` diffs fresh runs against it. Nothing has ever
compared the branch to `main`. The golden harness sees only what its single
fixture and single input format reach, and is blind by construction to anything
the pre-edit capture already baked in.

The goal is to build the missing cross-version comparison and invert the burden
of proof: for every output that differs, either it is a documented fix **and we
assert the direction** of the change (main wrong, this branch right), or it is a
regression and we say so.

Nothing is committed or pushed. All new scripts and outputs land untracked.

### Established (do not re-derive)

- `origin/main` == local `main` == `849c174`; `origin/alex-refactor` == local
  `HEAD` == `98097bc`; 54 ahead, 0 behind, clean tree. **The GitHub comparison
  can be done entirely locally** — no fetch, no divergence.
- `NAMESPACE` byte-identical. `R/classDefinitions.R` byte-identical.
  `src/*.cpp`/`*.h` identical (only committed `.o`/`.so` removed) — the compiled
  layer is out of scope.
- All 24 exports have identical formals and defaults on both branches.
- `inst/testdata/make_fixtures.R` has no `spatialGE` dependency, so one fixture
  generator can feed both versions. No `:::` usage anywhere in the harness, so it
  works against an installed package via `library()`.

### Three findings that reshape the plan

**A. NEWS is wrong: the ordering fix does change values.** `STenrich`'s
permutation null is `set.seed(seed)` then
`computeSubsampleSums(coords_mat, n_subsample=…, n_samples=reps)`
(`R/stenrich_helpers.R:80-81`), and `computeSubsampleSums` draws
`Rcpp::sample(all_indices, …)` — **row positions**, not names
(`src/spatialge_helpers.cpp:26-33`). Reorder `coords_mat` and the same seeded
draws select different spots, so the null distribution and therefore
`p_value`/`adj_p_value` genuinely differ. The observed statistic is safe
(`coords_mat[high_spots_bc, ]` is character-indexed), which makes this the sneaky
case: the statistic looks stable, the p-value moves. Separately, every reduction
over spots (`colSums`, `var`, the Rcpp VST) changes summation order under a
column permutation — integer counts are exact, but `@tr_counts` and
`@gene_meta$vst.*` are not, and `arrange(desc(vst.variance.standardized))` +
`topgenes` means a 1e-16 change can flip a rank and select a different gene set.
So the plan needs a float-noise bucket, an order-induced-value bucket, and an
explicit check on the *derived discrete decisions*.

**B. `STclust` labels are order-sensitive, which corrupts three downstream
functions.** `hclust()` on the fixture's regular 12×10 integer grid has a
distance matrix saturated with exact ties, and tie merge order depends on input
row order. `bl_pipeline` then feeds `stclust_spw0.025_dsplFalse` and
`cluster_labels[1]` as `annot=`/`domain=`/`ref=` into `STenrich`, `STgradient`
and `STdiff` — so a label permutation means those three are asked *different
questions* on the two branches. Sorting rows does not help. Needs a label
canonicalizer, a membership gate, and a pinned-annotation scenario.

**C. Warnings are invisible today.** `bl_capture()` catches only `error`. Delta 8
(`ab76c49`) is a *warning* that errored — structurally undetectable by the
existing harness. Warnings and messages must be first-class payload.

### The 15 intentional deltas to reconcile against

Reachable by the existing 26-step baseline (2):

| # | Commit | Change |
|---|---|---|
| 1 | `f278bd0` | `process_lists()`'s coordinate sort was a silent no-op on main (`arrange(!!array_col)` splices a character literal). Sorts by the **third column of the coords frame**. Reorders `@counts`/`@tr_counts` columns, `@spatial_meta` rows, and downstream row order. |
| 2 | `3cd478f` | `STenrich(annot=, domain=)` always errored on main (`` `.data` must be a data frame without row names ``). A captured error on main becomes a real result. |

Not reachable by the existing baseline (13):

| # | Commit | Change | Why unreachable |
|---|---|---|---|
| 3 | `a60651a` | Xenium MEX import (duplicate `import_visium` shadowing; missing `import_xenium_mex`; parquet-vs-csv detection) | only the in-memory data-frame format is used |
| 4 | `b7e9c73` | `STclust()` clobbered its own `verbose` | `verbose` never varied |
| 5 | `d8a2e37` | `rm_genes_expr` did not filter `@gene_meta` (`rownames()` on a tibble returns `"1".."200"`, so `grepl` matches nothing — a silent no-op, not an error) | a later stage re-syncs it; deliberately not in NEWS |
| 6 | `777bb59` | xlsx samples-file detection (`grepl("//.xlsx$", …)`) | format not exercised |
| 7 | `fbe02cf` | `cores=` ignored in `STList.R`/`transform_data.R` (main has a literal `cores = 1` at `transform_data.R:101`) | affects fork count only, never return values |
| 8 | `ab76c49` | undefined `h5_file` in the import warning path | h5 formats not exercised; **and warnings aren't captured** |
| 9 | `e7917bc` | `pseudobulk_heatmap()` inverted CV (`mean/sd` for `sd/mean`) → wrong gene set | `pseudobulk_samples()` refuses <4 samples; fixture has 3, so it is a **captured error on both branches**. The most user-visible NEWS change, verified by inspection only. |
| 10 | `b8de052` | `STdiff_volcano()` checked `is.numeric(x)` not `is.numeric(samples)` | numeric `samples` never passed |
| 11 | `6954044` | `STlist(rnacounts=NULL)` gave "argument is of length zero" | counts always passed |
| 12 | `ede00c8` | `non_spatial_de()` passed a `conditionMessage()` string to `summary.HLfit()`, poisoning an `mclapply` preschedule chunk; `count_cores()` clamps under `_R_CHECK_LIMIT_CORES_` | only via `STdiff(test_type='mm')` — which is the **default**, and has zero coverage today |
| 13 | `f05c855` | `STgradient`: `pval_warn` initialised to `NA`, closing a state leak across gene-loop iterations. **Not in NEWS** | needs a gene with post-filter n ≤ 1 |
| 14 | `5e64688` | doubled roxygen marker broke all 14 examples | docs, not runtime |
| 15 | `d3484a9` | vignette `thrane_tmp`/`visium_tmp` | docs, not runtime |

---

## Design

### Process topology

Five processes. The comparator has **no `spatialGE` on its library path at all** —
asserted, not assumed, so a stray S4 object fails loudly instead of silently
deserializing against the wrong class definition.

```
P1  producer  lib-main  rep=1  -> payload_main_1.rds + stdout_main_1.txt
P2  producer  lib-main  rep=2  -> payload_main_2.rds
P3  producer  lib-alex  rep=1  -> payload_alex_1.rds + stdout_alex_1.txt
P4  producer  lib-alex  rep=2  -> payload_alex_2.rds
P5  comparator (spatialGE not loadable) -> report
```

Producers share the conda dependency library, so the ~366 dependencies are
byte-identical across both runs. Each producer asserts
`dirname(find.package('spatialGE'))` and `packageVersion()`, and records every
`Imports` version into the payload header; the comparator **refuses to compare**
if headers differ. This gate matters because `spatialGE` has never been
*installed* in this env (only `load_all()`), so `main` has never run against
these dependency versions and any diff could be dependency drift.

**Never `devtools::load_all()` in a producer** — it shadows the installed library
with the alex working tree and silently makes both children "alex". Use
`library(spatialGE)` plus the assertions above. Run children with `--vanilla`
(or `R_PROFILE_USER=/dev/null` and `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE`), since
a local `.Rprofile`/renv autoloader would hijack `.libPaths()`; CLAUDE.md says
those are gitignored but they may exist locally.

Producer settings, both sides identical: `options(mc.cores = 1)` (pins the one
knob main leaves to `getOption("mc.cores", 2L)` — fidelity rule 3),
`Sys.unsetenv("_R_CHECK_LIMIT_CORES_")` (rule 12: alex clamps, main doesn't),
`Sys.setlocale("LC_COLLATE", "C")` (canonical sorting must not be locale-dependent
— the fixture deliberately contains `Gene$special` and `1BadGene`), `pdf(NULL)`
(`pseudobulk_heatmap()` calls `ComplexHeatmap::draw()`, which otherwise writes
`Rplots.pdf` into cwd), and a dedicated temp root — never bare `tempdir()`, since
`test-tests.R:11` does `unlink(tempdir(), recursive=TRUE)`.

### Getting both sources: `git archive`, symmetrically

`git archive <ref> | tar -x` into a temp dir (verified working). Chosen over
`git worktree add`, which writes refs and admin files into `.git`, marks `main`
as in-use, and needs a cleanup step that can leave a half state on abort. The
extracted tree has no `.git`, so nothing the harness runs can touch the
repository — the whole argument for a read-only-by-contract exercise.

Extract **both** refs this way, not just main; installing alex from the live
working tree would be an asymmetry. Then, mandatorily: `rm src/*.o src/*.so`
(main has these committed, and `git archive` sets mtimes to commit time so `make`
may consider them current and relink a `.so` built against an unknown ABI), plus
`R CMD INSTALL --preclean --with-keep.source=FALSE -l <verlib> <srcdir>` as a
second belt. `--with-keep.source=FALSE` keeps srcrefs out of function bodies,
removing a noise source.

Minor known non-issue: main's `DESCRIPTION` has `Copyright: ins/COPYRIGHT`
pointing at a nonexistent path (the typo alex fixes). `R CMD INSTALL` does not
validate that field; it is a check NOTE only.

### Reduce payloads in the producing process

Forced, not preferred: an `STlist` is S4 and `@gene_krige` holds gstat objects,
so deserializing in the comparator would require loading *some* version of
`spatialGE` and defeat the isolation. Likewise a `ggplot` carries `plot_env` and
layers that resolve into the producing namespace, so `ggplot_build()` must run
where the plot was made.

Each producer runs `xv_assert_plain(payload)` before `saveRDS` — a recursive walk
erroring on any environment, closure, S4 object or `externalptr`. **This is the
single most important correctness property of the harness.** Payloads are nested
named lists of atomic vectors and plain data.frames, with class strings recorded
as *data* (`class(x)` as character) rather than attributes, so a diff names an
exact path like `STenrich/value/adj_p_value`. Raw objects are also saved
per-version for optional deep-dives, loaded only in a process that has that
version installed; the comparator never touches them. Each producer self-checks
that `reduce()` is deterministic (reduce twice, `waldo::compare` must be empty).

Reused, not reimplemented: `bl_compare_ggplot_signature()`
(`dev/baseline/compare_baseline.R:80`) — it is version-neutral (touches only
`waldo`/`Matrix`/`sf`/`sp`/`spdep`/`ggplot2`). Extend, don't duplicate, per
`ARCHITECTURE.md`.

**Two comparator gaps must be closed** (in a wrapper — `dev/baseline/*` is never
modified, as that would invalidate `tests/golden/`):

1. **`@misc` is almost entirely uncompared.** `bl_compare_stlist()`
   (`compare_baseline.R:100-153`) walks `@counts`, `@tr_counts`,
   `@spatial_meta`, `@gene_meta`, `@sample_meta`, `@misc$sthet$listws`, and
   `@gene_krige` — nothing else. The package writes **10** `@misc` keys; eight are
   never compared, including `scaled_pbulk_mtx` (the direct input to the delta-9
   CV fix) and `sp_images` (the entire output of `load_images()`). An inherited
   false-clean hole; the comparison must walk all of `@misc`.
2. **`pseudobulk_heatmap()` returns a `HeatmapList`, not a ggplot or a
   `Heatmap`** — it calls `ComplexHeatmap::draw()` on the way out
   (`R/pseudobulk_samples.R:373`), so the matrix is at
   `hm@ht_list[[1]]@matrix`, *not* `hm@matrix`. `bl_compare()` falls through to
   generic `waldo`, and `normalize_for_baseline()` recurses only into *lists*, not
   S4 slots, so grobs and environments get compared raw — guaranteed noise.
   Extract per `ht` in `@ht_list`: `list(name, class, matrix, row_order,
   column_order, row_title, column_title)` and nothing else.

**The reducer must fail loudly on an unrecognized class**, never fall through to
a generic path that yields `NULL` — `NULL == NULL` is a clean, so a tolerant
reducer manufactures false passes. Assert every payload node is non-empty. Named
offenders that will otherwise bite: `HeatmapList` (S4, grobs + environments),
`gstat`/krige objects in `@gene_krige` (formulas + environments + `call`),
`attr(listw$neighbours, "call")` (a call referencing the producing namespace),
`sp::SpatialPointsDataFrame` (S4 → reduce to `coordinates()` + `@data`), and
`dgCMatrix` (S4 from `Matrix`; `drop0()` must happen in the producer).

Also extend the ggplot signature with facet class, panel layout and coord class;
the existing signature captures built data, labels, scale limits/classes and
geoms, so a faceting or coord regression would read as clean. And scrub every
character leaf of `tempdir()`, lib paths, extract dirs and PIDs, or the two
versions differ on every stored path.

`@spatial_meta` row names must be recorded as their own payload field, not
dropped — main's non-contiguous `1, 5, 9, …` names *are* the mechanism of delta 2.

### Canonicalization: sort on identity, never on values

Each rule individually toggleable so the report can say which rule absorbed which
diff.

| Object | Rule |
|---|---|
| `sample_order` | **not** canonicalized — fidelity rule 7 says sample order is user-visible and is not expected to change; canonicalizing would hide a regression |
| `@counts`/`@tr_counts` | columns by `libname`, rows by gene name; original `Dimnames` kept as separate raw fields |
| `@spatial_meta` / `@gene_meta` | sort rows by `libname` / `gene`; **keep column order** (rule 8: column position is API) |
| result tibbles | sort by **declared identity columns only** (below) |
| ggplot built data | multiset compare; sort by a **rounded** key (below) |
| `@gene_krige` | sort grid rows by `format(x, digits=15)` |
| cluster labels | relabel each cluster by `min(libname)` in it; compare as a partition |
| ggplot `group`/`PANEL` codes | relabel in order of first appearance after canonical sort |

Three places where "just sort it" is wrong:

- **Never include value columns in a sort key.** Sorting by all columns would
  reorder rows *to make values line up*, masking the very diffs we hunt. Declare
  identity keys per result: `STenrich` → `(sample_name, gene_set)`; `STgradient`
  → `(sample_name, gene)`; `STdiff` → `(sample, gene, cluster_1, cluster_2)`
  (`+ spatialmod` for spatial); `get_gene_meta` → `(sample, gene)`. Then **assert
  the key is unique**; if not, report `AMBIGUOUS_KEY`, fall back to whole-row
  multiset comparison, and state in the report that a value diff inside a
  duplicated key is undetectable. This is the honest answer to the
  `adj_p_value`-tie concern: sorting on identity columns makes the shipped
  tie-break order irrelevant rather than something to reconcile.
- **Never sort on raw doubles.** One 1e-16 difference flips the sort and the diff
  becomes a huge apparent permutation. Round the *sort key* to 12 significant
  digits; compare the *unrounded* values.
- **Order-sensitive algorithms get a bucket, not a sort rule.** No sort makes
  `STenrich`'s permutation p-value comparable (finding A). Declare these fields
  up front and route them to `ORDER_INDUCED_VALUE`.

Before any sorting, assert **dimensions and name multisets** are equal. Otherwise
canonicalization could mask a genuinely dropped or duplicated spot or gene.

### Passes: four, plus a positive invariant

- **Pass 0 — self-reproducibility, per version.** rep 1 vs rep 2, raw, in
  **separate processes** (a same-process rep 2 inherits mutated RNG state, a
  populated temp root, and warm caches; also `img_dir` collides). Must be empty.
  Every non-empty path goes on a *printed* noise allowlist. Inspected before any
  cross-version comparison is attempted. This matters concretely: because of
  delta 7, main's `transform_data(method='sct')` forks even when asked for
  `cores=1`, and fidelity rule 1 says forked children discard the RNG.
- **Pass A — raw cross-version**, minus the allowlist. Deliberately noisy.
- **Pass B — canonical cross-version.** Diffs here are value diffs.
- **Pass B′ — canonical at relative `tolerance = 1e-12`.** Surviving diffs are
  *material*; diffs in B but not B′ are float noise from reduction-order changes.
- **Invariant C — keyed alignment.** This *replaces* "pass A minus pass B", which
  does not work: `waldo` reports are not set-algebraic, and pass A contains value
  diffs as well as order diffs. Instead, positively: assert
  `setequal(libnames)` and `setequal(genes)` per sample, then compare
  `counts[gene, libname]` and `tr_counts[gene, libname]` **matched by name**, and
  `spatial_meta`/`gene_meta` **joined on** `libname`/`gene`. If that passes while
  pass A is non-empty for the same object, NEWS's "identical once reordered" claim
  is *positively demonstrated* for that object. Report the max relative delta on
  `tr_counts` and `gene_meta` as its own number — that is the empirical
  measurement of how true the claim is at float precision.

Report as a first-class derived check: the **`topgenes` selection set** per sample
on each branch. Identical set ⇒ float noise was harmless. Different set ⇒
everything downstream is legitimately incomparable, and must be stated rather
than triaged row by row.

**Cores policy.** Primary comparison at `cores=1` only. `cores=3` is secondary
and explicitly lower-confidence: fork *topology* differs between versions at the
same nominal `cores`, so comparing there mixes two independent variables.

---

## Files to create (all untracked, under `dev/xver/`)

Nothing under `R/`, `tests/`, `inst/`, or `dev/baseline/` is modified — that is
what preserves `tests/golden/` and the `SPATIALGE_GOLDEN` gate. Add
`dev/xver/out/` to `.gitignore`.

| Path | Role |
|---|---|
| `00_gate.sh` | env preconditions: conda env, TLS gate, clean tree, both SHAs, `_R_CHECK_LIMIT_CORES_` unset, `detectCores()` sane |
| `10_build_libs.sh` | `git archive` both refs, strip `src/*.o`/`*.so`, `R CMD INSTALL --preclean` into two libs |
| `reduce.R` | version-neutral extractors; must not name a single `spatialGE` symbol |
| `canon.R` | canonicalization rules, on reduced payloads only |
| `scenarios.R` | the scenario matrix below; exported functions only |
| `produce.R` | sources `inst/testdata/make_fixtures.R`, `dev/baseline/pipeline.R`, `dev/baseline/compare_baseline.R` **by absolute path from the alex working tree**; applies the RNG shim; writes payload + stdout |
| `classify.R` | machine-readable expectation table |
| `compare.R` | the comparator (P5) |
| `run_all.sh`, `README.md`, `out/` | orchestration, docs, artifacts |

Critical rule: **both producers source the same copies of `make_fixtures.R` and
`pipeline.R`, from the alex working tree, never from the installed package.**
Main's tree has neither file, and `inst/testdata/` would install to a different
path anyway. One copy for both sides is the only way to guarantee an identical
call sequence, and it is safe because the fixture generator has no `spatialGE`
dependency and all 24 exports have identical formals.

Where a different fixture is needed, the driver **shadows** a generator in the
global environment after `source()` — this works because `fixture_base_data()`
resolves `fixture_sample_names()` lexically at call time and `bl_pipeline`'s
environment is `globalenv()` (assert that, don't assume it). `capture_baseline.R`
never applies these shadows, so the committed golden files stay valid.

### Fixing the `load_images` RNG false positive

`dev/baseline/pipeline.R:156-158` writes the PNGs via `fixture_write_image()`,
which calls `runif()` *outside* `bl_capture()`'s `set.seed()` — so the images
depend on RNG state left by the preceding 25 steps, which can differ between
branches. Useful corollary: since `bl_capture()` re-seeds before *every* step,
this loop is the **only** RNG-state-order dependency in the whole pipeline.

Shim `fixture_write_image()` in the producer to seed deterministically from the
output path. Belt and braces, independent of the shim: reduce a loaded image to
`(class, dim, range, digest)`, so a shim failure yields one compact flagged field
instead of a wall of pixel diffs.

### Scenario matrix

| ID | Scenario | Covers |
|---|---|---|
| S0 | `bl_pipeline(cores=1)` unmodified — the anchor | 1, 2 |
| S0b | `bl_pipeline(cores=3)` — secondary, lower confidence | — |
| S1 | `STlist()` on each of the 11 `fixture_materialize_all()` variants | **3, 6, 8** |
| S2 | `STlist` → `transform_data(log)` → `filter_data(rm_genes_expr='^Grad')` | **5** |
| S3 | `STclust(verbose=TRUE)` / `verbose=2L`, `cores=1` only | **4** |
| S4 | shadow `fixture_sample_names()` → 5 samples; pseudobulk chain to `pseudobulk_heatmap(hm_display_genes=5)` | **9** |
| S5 | `STdiff_volcano(samples=1)` (numeric) | **10** |
| S6 | `STlist(rnacounts=NULL)` and other NULL permutations | **11** |
| S7 | all 24 exports with no args / `x=NULL` / wrong type; compare `conditionMessage` only | new — catches the `assert_*`/`resolve_samples` refactor across six plotting functions. Zero coverage today, cheap, high yield |
| S8 | injected `annot`: 119 spots `'R'`, 1 spot `'X'`; `STgradient(ref='R')` | **13** |
| S9 | small driver-built STlist; `STdiff(test_type='mm', topgenes=10)` at `cores=1` | **12** partially — and `mm` is the *default*, currently unverified |
| S10 | `options(mc.cores=2)` vs `1`, both versions | **7** |
| S11 | pinned-annotation variants of `STenrich`/`STgradient`/`STdiff`, annot from `grid$domain` | isolates finding B's cascade — makes those diffs interpretable at all |
| S12 | real Thrane Visium data: `STlist` → `transform_data` → `SThet` → pseudobulk | 1 and 9 on real unsorted coordinate files; the fixture's shuffling is synthetic |
| S13 | **non-default argument enum sweep**: `STclust(ks=<numeric>)`, `STgradient(out_rm=TRUE, robust=FALSE, distsumm='avg', log_dist=TRUE)`, `SThet(method='geary', k=)`, each `filter_data` stage, `STplot(plot_meta=)` | new — the baseline uses exactly one parameter set per function, so whole branches are uncovered. `STclust(ks=numeric)` means `get_hier_clusters_ks()` is **never called**. Untaken branches in `STdiff_prepare/fit/assemble`, `filter_data_stages` and `stlist_readers` are precisely where refactor bugs hide |
| S14 | **`samples=` sweep**: numeric index, exact name, and a substring matching two samples | new — the refactor replaced inline guards with `resolve_samples(numeric_handling=, char_transform=, post_check=)`. A wrong flag at one call site still returns all 3 samples on this fixture, so this is a permanent false clean otherwise. The unanchored `grep(paste0(samples, collapse='|'))` substring semantics are deliberately preserved and load-bearing |

S3 note: `verbose` output is only observable via **process stdout** (fidelity
rule 11: `system(sprintf('echo …'))` bypasses `capture.output`). Producers
redirect whole-process stdout to a file; the comparator diffs transcripts after
scrubbing timing lines. `cores=1` only — from a fork that output is
unattributable.

S4 note: adding samples only *appends* `rpois`/`sample.int` draws, so
`alpha`/`bravo`/`charlie` stay bit-identical to the 3-sample fixture.
`calc_umap=FALSE` is mandatory (rule 1: `uwot::umap` is unseeded and
multithreaded).

### Directional oracles: difference is not evidence

A fixed function *differing* proves nothing about which side is correct. Each
expected fix gets an oracle computed independently of the package:

| Delta | Oracle |
|---|---|
| 1 | Assert alex's per-sample coords are non-decreasing in the third coords column and main's are not; then invariant C. |
| 2 | main errors with the row-names message; alex returns a tibble with documented columns and ≥1 row. |
| 9 | First assert `scaled_pbulk_mtx` is identical across versions (else the oracle is void). Then independently compute both statistics and assert `alex_rownames == top5_by(sd/mean)` **and** `main_rownames == top5_by(mean/sd)`. A falsifiable prediction, not a diff. |
| 6 | main errors; alex canonically equals the csv-samples-file variant built by the same version. |
| 3 | main errors; alex canonically equals its own `xenium_h5` variant. |
| 4, 10, 11, 13 | Assert the specific observable (captured stdout; volcano label; error message vs success; `pval_comment` column) rather than a whole-object diff. |

For delta 13 the prediction is concrete: main's `rm(list=grep(…))` does not match
`pval_warn`, and in the default `robust=TRUE` branch `pval_warn = NA_character_`
sits inside `if(nrow(df_gene_range) > 1)`. So with ≤1 non-reference spot, **main's
first gene should error with `object 'pval_warn' not found`** while alex returns
`pval_comment=NA`. If main does not error, delta 13 is unreachable and stays
unverified — the report says so.

---

## Triage: fail closed

`classify.R` holds one record per expected difference (id, commit SHA, scenario,
path regex, class, rationale). Classes evaluated in order; **first match wins; no
match means REGRESSION**:

| Class | Condition | Verdict |
|---|---|---|
| `NOISE` | path appears in Pass 0 for either version | IGNORED — but every such path is **printed**; silent allowlists are how harnesses rot |
| `NON_COMPARABLE` | on the declared list below | N/A |
| `INCOMPARABLE_CASCADE` | the per-sample cluster **partition** differs between versions, so `annot=`/`domain=`/`ref=` differ and the two branches were asked different questions | N/A for every downstream result — explicitly **not** `ORDER_ONLY`. This gate runs *before* any downstream verdict is believed. Same rule if the `topgenes` selection set differs, or if `bl_pipeline`'s `st1..st5` error-fallback chain took a different branch on the two versions |
| `EXPECTED_FIX` | matches a table record naming a commit **and its oracle passes** | PASS |
| `ORDER_ONLY` | in Pass A, absent from Pass B, **and** invariant C holds for the containing object | PASS, attributed to `f278bd0`. The conjunction matters: absence from B alone could mean a canonicalization rule is masking a real diff |
| `ORDER_INDUCED_VALUE` | in Pass B′ and on the declared order-sensitive list (`STenrich$p_value`/`adj_p_value`, `stclust_*` labels, anything downstream of a changed `topgenes` set) | PASS-WITH-CAVEAT + an explicit report line that **this contradicts NEWS and NEWS should be amended** |
| `FLOAT_NOISE` | in Pass B, not in Pass B′ | PASS-WITH-CAVEAT; report max relative delta and whether any derived discrete decision changed |
| `REGRESSION` | everything else, or a mapped delta whose **oracle fails** | FAIL |

Two structural requirements:

- **Exhaustiveness assertion**: `sum(counts per class) == total observed diff
  paths`. If the classifier cannot account for every path, it says so.
- **Vacuous-coverage metric**: report "N of M steps compared only an error
  string" and treat that count as coverage debt, not green. `bl_capture` stores
  only `conditionMessage`, and `pseudobulk_heatmap`/`pseudobulk_dim_plot`
  short-circuit to the literal `"pseudobulk_samples did not succeed"` on *both*
  branches — an `IDENTICAL` verdict that tests zero code. That is exactly how
  delta 9 stayed unverified.
- **`MISSING_EXPECTED`** — the inverse check. Every record in the expectation
  table must actually be *observed*. An expected fix producing no diff means
  either the scenario doesn't reach it or the fix doesn't work → FAIL, or
  `COVERAGE_GAP` if the scenario didn't run. This is the check that would have
  caught "delta 9 is entirely unverified" a month ago, and it mechanically
  enumerates the honest gap list instead of leaving it to memory.

---

## Verification

```sh
source /opt/anaconda3/etc/profile.d/conda.sh && conda activate spatialge-dev
export R_MAKEVARS_USER="$(pwd)/dev/environment/Makevars"
```

Each step gates the next:

1. **TLS gate — currently FAILING.** `curlGetHeaders()` from conda R returns
   `libcurl error code 60`. The intercepting root `CN=ca.moffitt.goskope.com` *is*
   in the System keychain (verified, 1 match), so the `dev/environment/README.md`
   recipe will work. Build the bundle, export `CURL_CA_BUNDLE`/`SSL_CERT_FILE`,
   then require `HTTP/2 302` from the melanoma_thrane URL. If it does not pass,
   S12 records a **`COVERAGE_GAP` failure, not a skip** — a `tryCatch` that
   reports "skipped" is how a coverage gap becomes a pass.
2. **Fixture dependencies present in the shared library**: `writexl`, `png`,
   `jsonlite`, `hdf5r`/`rhdf5`, `arrow` as applicable. `writexl` and `waldo` are
   in **alex's** Suggests only, so installing main will not surface a gap — and
   `fixture_write_samplesfile(format='xlsx')`, the fixture for delta 6, needs
   `writexl`.
3. **Both versions install and load**, each reporting its expected
   `packageVersion()` from its own library. **Read main's install log** — it has
   never been compiled in this conda toolchain.
4. **Smoke-run S0 against `main` alone**, and triage every error/warning not on
   the 15-delta list *before anything else*. Main (1.2.2) has never run against
   `ggplot2` 4.0.3 and current `dplyr`/`Matrix`. The plot code is functionally
   unchanged between branches, so dependency drift should hit both sides equally
   — but "should" is not "does", and an unexplained main-only failure here
   invalidates everything downstream.
5. **The existing golden gate still passes**, confirming nothing in the working
   tree drifted:
   ```sh
   Rscript -e 'devtools::load_all(); testthat::test_local()'
   SPATIALGE_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_local()'
   ```
6. **Pass 0 clean or fully explained** before any cross-version verdict.
7. **Self-check that the comparison can fail**: compare main against main with
   one value perturbed, and confirm it is reported. A comparison that has never
   failed is not evidence.
8. **Coverage assertion**: every one of the 15 deltas marked
   covered-and-verified or explicitly unverified with a reason. No silent gaps.

## Deliverable

`dev/xver/REPORT.md`: per-step verdict table, the order-only set as evidence for
NEWS's ordering claim, each oracle's result, the `MISSING_EXPECTED` list, and any
regression.

Then, only for what the comparison confirms, proposed wording in `NEWS.md` and
`dev/ARCHITECTURE.md`, left uncommitted in the working tree. Two candidates are
already identified:

- **NEWS's ordering entry is inaccurate** (finding A) — "no expression values
  changed" is false for `STenrich`'s permutation p-values, and arguably for
  `@tr_counts`/`vst.*` at float precision. The comparison will quantify it; the
  entry should be amended to say so.
- **Delta 13** (`f05c855`, the `STgradient` `pval_warn` state leak) is a real
  behavior change documented nowhere.

## Known to remain unverifiable, and why

- **Delta 7's actual effect** — fork count appears in no return value. S10
  verifies the *absence* of an output change; it cannot verify the presence of
  the parallelism change.
- **Delta 12's poisoning guard** — needs a non-converging `spaMM` fit sharing a
  preschedule chunk with converging ones. `ARCHITECTURE.md` fidelity rule 2
  explicitly says fixture policy is the opposite ("design fixture data so every
  fit converges"). If S9 produces no failing fit, the guard is unexercised.
  Delta 12's `count_cores()` clamp is fork-count-only and non-exported → verify
  single-version and say so.
- **`STdiff(test_type='mm')` at realistic size** — ~21 s, ~27 GB per cluster. S9
  runs a shrunken version; the released default on realistic data is uncovered.
- **Deltas 14 and 15** — not runtime outputs. Cheap side check via
  `tools::Rd2ex(commentDonttest = FALSE)`; worth running, not part of the
  comparison.
- **Seurat input path, geomx, and the other `detect_input_cases.R` branches** —
  no fixtures; `SeuratObject` is Suggests-only.
- **Pixel-level plot identity** — only built data, labels, scales, facets, coords.
  Two plots with identical built data can still render differently.
- **`pseudobulk_samples(calc_umap=TRUE)`** and **`sct` at `cores>=2`** —
  structurally non-reproducible (rule 1). Excluded, not "verified equal".
- **Windows** — `resolve_cores` forces 1 core; that branch is unreachable here.
- **Non-comparable by construction**: `cat`ted `difftime` timers (rule 10);
  environments, closures, `srcref`, `call`/`terms`; `conditionCall` (only
  `conditionMessage` is compared, since `assert_*` extraction legitimately moves
  the call frame); ComplexHeatmap grobs; `gstat` `call` slots; temp/lib paths and
  PIDs (scrubbed); wall clock, memory, fork counts.
- **The dependency dimension** — this holds *one* dependency set fixed. It proves
  nothing about how 1.2.2 behaved on the dependency versions users actually had.

Meta-caveat: this compares the two branches *as executed by identically-sourced
driver code*. It cannot detect a change both branches agree on — most importantly
any input format or argument combination no scenario reaches. The
`MISSING_EXPECTED`/`COVERAGE_GAP` verdict is what forces that set to be
enumerated rather than assumed empty.
