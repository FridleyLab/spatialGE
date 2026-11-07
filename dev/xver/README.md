# `dev/xver/` — cross-version comparison harness

Compares an installed **`main` (1.2.2)** against an installed
**`alex-refactor` (2.0.0)** and classifies every differing output, so that each
difference is either a documented fix whose *direction* is asserted, or a
regression.

Results: **`REPORT.md`**.

## Why this exists

`tests/golden/` can only compare this tree against *itself*. The baseline was
captured from the same working tree before the `R/` edits began, so it cannot see
anything already baked into that capture. It also reaches only one input format,
one parameter set per function, three samples, and no warnings — which is how the
`pseudobulk_heatmap` coefficient-of-variation fix, the most user-visible change in
`NEWS.md`, went two releases and an entire refactor without ever being executed.

## Design, in the order it matters

**Isolation.** Both versions are the same package name, so they cannot coexist in
one process. Each producer runs in its own `Rscript --vanilla`, asserts it loaded
the intended library *and* version, and never uses `devtools::load_all()` (which
would silently make both children `alex`). Results are reduced to plain data
*inside* the producer; the comparator runs with `spatialGE` absent from its library
path, so a leak fails loudly instead of deserializing an S4 object against the
wrong class table. `xv_assert_plain()` enforces that before anything is written.

**Fail loudly, not quietly.** The reducer errors on a class it has no extractor
for. A reducer that falls through and returns `NULL` manufactures clean
comparisons, because `NULL == NULL`.

**Two passes, because order changed.** The branch's ordering fix (`f278bd0`) means
almost every output differs in row or column order, which would bury everything
else. So each result is reduced twice: as produced, and canonicalized (spots sorted
by `libname`, genes by name, result tables by declared *identity* columns only —
never by value columns, which would reorder rows to make values line up and hide
the very differences being hunted). A third pass applies a scale-relative 1e-12
tolerance to separate float noise from material change.

**Difference is not evidence.** A fixed function differing says nothing about which
side is correct. `50_oracles.R` recomputes the right answer from the input data,
independently of both packages, and asserts which version matches.

**Fail closed.** Anything matching no expectation record is a `REGRESSION`, and the
classifier asserts it accounted for every observed path. The table is also used in
reverse: a claimed fix that produces *no* observed difference is
`MISSING_EXPECTED`, because that means either the scenario does not reach it or the
fix does not work.

## Files

| File | Role |
|---|---|
| `env.sh` | paths, CA bundle, `LC_COLLATE`, `mc.cores`. Source, don't run |
| `00_gate.sh` | preconditions: conda env, clean tree, both SHAs match `origin`, TLS, deps |
| `10_build_libs.sh` | `git archive` both refs, strip committed `src/*.o`/`.so`, install into two libraries |
| `common.R` | per-producer determinism, version assertions, the fixture-image RNG shim, `xv_capture()` |
| `reduce.R` | version-neutral extractors, canonicalization, `xv_assert_plain()` |
| `scenarios.R` | S0–S15 |
| `produce.R` | runs a scenario set against one version, writes a payload |
| `classify.R` | the expectation table, matched on path *and* value |
| `compare.R` | passes 0/A/B/B′, gates C/P/S, classification, report data |
| `run_all.sh` | orchestrator |
| `20_smoke.R`, `21_test_reduce.R` | single-version smoke and reducer self-test |
| `30_order_probe.R` | proves which differences are order-induced, using one version |
| `40_stdout_probe.sh` | the `verbose` fix, observable only in process stdout |
| `50_oracles.R` | directional oracles |

Artifacts (payloads, logs, installed libraries) go to `$XVER_ROOT`, default
`~/.spatialge-dev/xver` — outside the repo, since they are large and
machine-specific.

## Scenarios

`S0`/`S0b` are `bl_pipeline()` unmodified at `cores=1`/`3` — the anchor, and the
same code the committed golden gate runs. The rest exist because the anchor reaches
only 2 of the 15 deliberate behavior changes.

| | Covers |
|---|---|
| `S1` | all 11 input formats → Xenium MEX, xlsx samples file |
| `S2` | `rm_genes_expr` reaching `@gene_meta` |
| `S3` | `STclust(verbose=)` (process stdout only) |
| `S4` | 5-sample pseudobulk chain → the inverted-CV fix |
| `S5` | `STdiff_volcano(samples=<numeric>)` |
| `S6`, `S7` | invalid-input and error-message sweeps across all 24 exports |
| `S8` | `STgradient` `pval_warn` state leak |
| `S9` | `STdiff(test_type='mm')` — STdiff's *default*, otherwise uncovered |
| `S10` | `cores=` and `mc.cores` invariance |
| `S11` | pinned annotation, isolating `STclust`'s order sensitivity |
| `S12` | real Thrane Visium data |
| `S13` | non-default arguments (`ks=` numeric, `method='geary'`, `distsumm='avg'`, each filter stage) |
| `S14` | `samples=` as index, exact name, and substring |
| `S15` | multi-group `.h5`, the only way to reach the `h5_file` warning fix |

## Extending it

- New `@misc` key, or a function returning a new class → add an extractor to
  `reduce.R`. It will error until you do, which is the intent.
- New deliberate behavior change → add a record to `classify.R` naming its commit,
  and a scenario that reaches it. The `MISSING_EXPECTED` check will fail until the
  scenario actually produces the difference.
- A difference you believe is order-induced → prove it with `30_order_probe.R`
  before classifying it as such.
