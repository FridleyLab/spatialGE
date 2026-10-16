# A pkgdown-only vignette for on-disk STlists

## Context

`STlist_to_disk()` / `STlist_to_memory()` / `is_backed()` / `validate_store()`
shipped in 2.0.0 with one `NEWS.md` paragraph, no `@examples` on any of the four
reference pages, and no tutorial. It is the one feature in the package a user
cannot discover by reading `?`-help, because the thing that needs explaining is
a *workflow* (build in memory, move to disk, analyze, come back) rather than any
single argument. So: a new article on the pkgdown site, `.Rbuildignore`d like
three of the four existing vignettes, so it reaches readers without adding a
network-dependent vignette to the CRAN tarball.

Writing it surfaced three defects that have to be fixed first, or the article
would be a tour of workarounds. All three are reads of `@counts`/`@tr_counts`
that bypass `R/utils_accessors.R`, and `base::nrow()`/`rownames()` return `NULL`
on an `sge_manifest` rather than erroring — so each fails by *silently* reporting
nothing, which `.claude/CLAUDE.md` calls out as worse than an error:

| Call | On a backed STlist today |
|---|---|
| `x` / `print(x)` / `summary(x)` | prints `sample_117d ( ROIs\|spots\|cells x  genes)` — blank dimensions, no warning |
| `STplot(x, genes='NDRG1')` | warns "None of the requested genes are available" per sample, then stops with "Data was not found in the specified slot of this STList." |
| `distribution_plots(x, genes='NDRG1')` | messages "does not contain NDRG1" for every sample, returns plots that only fail when printed |

`STlist_to_disk()`'s own roxygen claims "All package functions continue to work,
because slot reads go through internal accessors." These three are why that is
currently false. Printing the object is also the first thing both existing
tutorials do after building one (`tnbc`, `melanoma`), so the vignette cannot
avoid it.

Scope confirmed with the user: fix those three, then write the article, then a
full `pkgdown::build_site()`.

## 0. Check the existing work in first

`R/classDefinitions.R` (+279) and `R/distribution_plots.R` (+4) — two of the
four files to edit — already carry uncommitted changes, so my commits cannot
cleanly separate from them. They belong to one coherent body of finished refactor
work already written up in `NEWS.md` (+126): 28 modified files (the `dim()`
base-contract rework, `stop()` → `raise_err()` conversions with 27 new
`inst/err.csv` codes, `STplot_interpolation`/`filter_data`/`transform_data`/
`utils_parallel` changes, hand-edited `man/`, a 107-line `.github/workflows/test.yml`
expansion) plus 7 new `tests/testthat/` files and `dev/ci/`.

So: one checkpoint commit of exactly that, before anything else. I will survey
the diff and **show the proposed message for approval before committing** — it
is not my work to characterize unprompted.

Excluded from it: `docs/` (the full rebuild in step 8 overwrites the partial
`build_reference()` sitting there, and the complete site lands in commit 6) and
`plans/golden-finding-boot.md` (this file).

## 1. Add the materializing slot accessor

`R/utils_accessors.R`, after `sample_tr_counts()` (line 31). Both plotting fixes
branch on `data_type` (`'tr'`/`'raw'`), which is a different vocabulary from the
accessors' `slot=` (`'counts'`/`'tr_counts'`); `sample_slot_mtx()` already does
that dispatch but returns the raw element, so it only serves the structural
accessors. The materializing counterpart, with `genes=` forwarded:

```r
sample_slot_counts = function(x, sample, slot=c('counts', 'tr_counts'), genes=NULL){
  slot = match.arg(slot)
  if(slot == 'counts') sample_counts(x, sample, genes=genes) else sample_tr_counts(x, sample, genes=genes)
}
```

Not exported, so no `NAMESPACE`/`man/` churn and no `devtools::document()`.
Name is unused in `R/`, `tests/` and `NAMESPACE`.

## 2. `show()` / `summary()` — report real dimensions

`R/classDefinitions.R` lines **50, 55, 90, 95** (four byte-identical lines, two
in `show`, two in `summary`, one per `>= 11`-samples branch). Substitute only
the two reads:

```r
# before
cat(paste0('\t', i, " (", ncol(object@counts[[i]]), ' ROIs|spots|cells x ', nrow(object@counts[[i]]), ' genes)\n'))
# after
cat(paste0('\t', i, " (", sample_n_cells(object, i, 'counts'), ' ROIs|spots|cells x ', sample_n_genes(object, i, 'counts'), ' genes)\n'))
```

Already the idiom in this file — the `dim()` method is
`as.integer(sample_dim(x, names(x@counts)[1], 'counts'))`.

## 3. `plot_spatial_expression()` — read the gene slice through the accessors

`R/plot_spatial_expression.R`. Keep `counts` as the accumulator (lines 79-89
depend on it: sample-dropping, the per-sample `warning()`, the
`rlang::is_empty(counts)` stop). Only carry the slot *name* alongside it and
read from `x`:

- lines 47-50: add `data_slot = 'tr_counts'` / `data_slot = 'counts'` to the two
  existing branches
- line 63: `if(any(rownames(counts[[i]]) %in% genes)){` →
  `if(any(sample_genes(x, i, data_slot) %in% genes)){`
- line 75: `expandSparse(counts[[i]][rownames(counts[[i]]) %in% genes, , drop=FALSE], transpose=TRUE)` →
  `expandSparse(sample_slot_counts(x, i, data_slot, genes=genes), transpose=TRUE)`

`plot_spatial_geneset.R:52` is the reference implementation that already does
this correctly, and is why `STplot(genes=list(...))` works on a backed object
while `STplot(genes='GENE')` does not.

## 4. `distribution_plots()` — same change in both helpers

`R/distribution_plots.R`, in `gene_plots()` (≈129-156) and
`cluster_gene_plots()` (≈226-254). Replace `expr_tmp = x@tr_counts` /
`expr_tmp = x@counts` with `data_slot = 'tr_counts'` / `data_slot = 'counts'`
(leaving the `p_title`/`ax_title` assignments alone), then in each loop:

- `if(!(gene %in% rownames(expr_tmp[[i]])))` → `if(!(gene %in% sample_genes(x, i, data_slot)))`
- add `mtx_tmp = sample_slot_counts(x, i, data_slot, genes=gene)` in the `else`
  branch and read `colnames(mtx_tmp)` / `mtx_tmp[gene, ]` in place of
  `colnames(expr_tmp[[i]])` / `expr_tmp[[i]][gene, ]`

`mtx_tmp` is this repo's existing name for a materialized one-sample matrix
(`summarize_STlist.R:53`, `filter_data_stages.R:147`, `STenrich.R:170`).
`expr_tmp` then has no remaining uses in either helper.

### Why steps 2-4 are output-neutral in memory

For a non-manifest element, each accessor reduces **by substitution** to the
expression it replaces: `sample_n_cells(x, i, 'counts')` → `base::ncol(x@counts[[i]])`,
`sample_genes(x, i, 'tr_counts')` → `base::rownames(x@tr_counts[[i]])`, and
`sample_tr_counts(x, i, genes=g)` → `m[rownames(m) %in% g, , drop=FALSE]` —
literally the same logical-mask expression already at both call sites, so
fidelity rule 9's matrix-row order is preserved by construction (it drives
`qp_list` key order downstream, so it is not cosmetic). Absent sample names,
`NA_character_` from an out-of-range numeric index, and an empty `@tr_counts`
all reduce to `NULL` through `[[` exactly as through a list slice, so the
"genuinely absent" branches are untouched.

**One behavior does move**, and it should be stated in the commit message: a
duplicated `samples=` (e.g. `STplot(x, genes='NDRG1', samples=c('a','a'))`,
never de-duplicated) currently has iteration 2 read the data.frame iteration 1
wrote back into `counts`, producing a false "None of the requested genes are
available" warning and, for a single sample, a false stop. Reading from `x`
fixes that. Not golden-covered (`bl_pipeline()` passes `samples=NULL`) and no
test exercises it.

Golden surfaces: `show`/`summary` are not in `bl_pipeline()` at all;
`STplot(st5, genes=c("Gene1","GradGene1"), data_type="tr")` and
`distribution_plots(st5, genes=c("Gene1","GradGene1"), data_type="tr")` are, and
both take the all-genes-present path. **So no commit here regenerates
`tests/golden/`.** The golden gate is run as proof — a diff there means the fix
is wrong, not that the baseline needs recapturing.

## 5. Tests

`tests/testthat/test-store_parquet.R`, whose premise is already parity — one
block per fix, using the file's existing `store_fixture()` / `backed_fixture(st, name)`:

- `show()`/`summary()`: `expect_identical(capture.output(methods::show(bk)), capture.output(methods::show(st)))`,
  plus one literal assertion that the real numbers appear (the current output is
  blanks, so an identity check alone would pass against another blank).
- `STplot(genes=)`: compare `names()` and per-panel rendered data both ways, for
  `data_type='tr'` and `'raw'`; keep assertions that an absent gene still warns
  per sample and then stops.
- `distribution_plots(genes=)`: same, plus the `plot_meta=` + `genes=` path
  through `cluster_gene_plots()` on an `STclust`ed fixture, and that an absent
  gene still messages.

Plots are compared through a local helper that `set.seed()`s and then
`ggplot2::ggplot_build()`s, because `ggforce::geom_sina` (the
`plot_type='violin'` default) draws fresh jitter inside the build —
`dev/baseline/compare_baseline.R:94-98` handles it the same way.

`test-graceful-errors.R` gets one test for the *silent* half of the rule (a
present gene must never be reported absent on a backed object). **Its `leaked`
vector is not touched** — none of these three surfaces as an internal error
string, and the backed failures that do leak internals are `transform_data()`
and `filter_data(verbose=TRUE)`, which stay documented-as-ordering below.

## 6. Correct the row-group pruning claim

`R/store_parquet.R:16-18`, `R/utils_accessors.R:16-17` and the store section of
`dev/ARCHITECTURE.md` all say a `genes=` read "touches only the row groups it
needs". `sge_store_materialize()` calls `arrow::read_parquet(m$path, as_data_frame=FALSE)`
and *then* `dplyr::filter()`s in Acero, so the whole file is read — pruning
would need `open_dataset()`. The saving is one full-size `dgCMatrix` allocation,
not I/O. Comment-only commit, fixed before the vignette because the article's
performance wording depends on getting this right, and because
`dev/ARCHITECTURE.md`'s own rule is to trust the code and fix the doc.

## 7. The vignette

New `vignettes/ondisk_stlists.Rmd`; add `^vignettes/ondisk_stlists\.Rmd$` to
`.Rbuildignore`. Sorts second in the auto-discovered articles index (after
`basic_functions_vignette`), which is the right slot for a follow-on topic.

House conventions to match exactly, from the existing four:

- YAML: the same five fields in the same order — `title`, ``date: '`r format(Sys.Date())`'``,
  `output: rmarkdown::html_vignette`, `bibliography: vignette_references.bib`,
  and the `vignette:` block with `%\VignetteIndexEntry{spatialge_tutorial_ondisk}`.
  No `description:` field — none of the four has one.
- First chunk `vig_settings` with `include=F`, `knitr::opts_chunk$set(collapse=TRUE, comment="#>", fig.width=8, fig.height=5)`
  and `options(rmarkdown.html_vignette.check_title=F)`.
- The `tryCatch` + `net_check <<- TRUE/FALSE` download guard from
  `basic_functions_vignette.Rmd:91-99` with `eval=net_check` on every
  downstream chunk. **Required** — the other two tutorials download with no
  guard, and a full `build_site()` aborts if a download fails.
- Code style: `<-` for objects, `=` with no surrounding spaces for arguments,
  single-quoted strings, `cores=2`, `warning=F`/`message=F` chunk options,
  snake_case chunk names.
- Prose style: "we"/"users", prose before each chunk ending in a colon, prose
  after interpreting the result, `***`spatialGE`***` at first mention, gene
  names italicized, `<details><summary>**...**</summary>` for asides, closing
  `## References` + `<div id="refs"></div>` and a Session Info `<details>`.

Data: `tnbc_bassiouni`, the same eight Visium samples as
`basic_functions_vignette.Rmd`, so the article reads as a continuation and the
`object.size()` numbers are worth quoting. `@bassiouni_2023` is already in
`vignette_references.bib`.

Outline:

1. What problem this solves — an STlist holds every sample's counts in memory at
   once; the store trades that for I/O. States up front that it buys capacity,
   not speed.
2. Building the STlist, then `filter_data()` and `transform_data()` — with the
   ordering stated as the deliberate workflow it is, since both rewrite the
   matrices and must run before the move.
3. `STlist_to_disk()`, `is_backed()`, printing the object, `summarize_STlist()`.
4. What is actually on disk: `list.files()`, the `print.sge_manifest` one-liner,
   `arrow::read_parquet()` of one file showing the `(gene_idx, cell_idx, value)`
   triplets, and why COO sorted by gene.
5. Memory, measured in the document with `object.size()` rather than quoting the
   fixture's 888 KB → 72 KB.
6. Analysis: `SThet()`, `STclust()`, `STplot()` for a gene and for clusters,
   `distribution_plots()` — with an explicit equality check against the
   in-memory result, and one honest timing comparison.
7. `validate_store()`, `validate_store(deep=TRUE)`, and `STlist_to_memory()`
   with an `identical()` round-trip check.
8. What to keep in mind — the store is a directory of files, not part of the R
   object, so `saveRDS()` records paths (`@misc$store$portable`, `FALSE` under
   `tempdir()`); transform and filter before moving; come back to memory before
   subsetting with `x[i]`, which leaves `@misc$store$sample_order` stale and
   makes `validate_store()` raise.

No `NEWS.md` entry — `NEWS.md` has never mentioned a vignette.

## 8. Rebuild the site

```sh
export CURL_CA_BUNDLE=/Users/4471592/.ssl/combined_cacert.pem
export SSL_CERT_FILE=/Users/4471592/.ssl/combined_cacert.pem
Rscript -e 'cat(curlGetHeaders("https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=")[1])'   # expect HTTP/2 302
Rscript -e 'pkgdown::build_site()'
```

The CA-bundle export is mandatory on this network (`dev/environment/README.md`)
— conda R's libcurl does not trust the intercepting root, and without it every
`download.file()` fails with libcurl error 60 while `/usr/bin/curl` succeeds.
Verified working above.

The committed site is stale (shows 1.2.0, built with pkgdown 2.1.1; 2.2.1 is
installed), so a full build churns all ~117 tracked HTML files and re-runs the
two unguarded tutorials. Expect 30-90 minutes. If an article fails, I report it
rather than committing a half-written `docs/`.

## 9. Commits

One per numbered step, staging only the paths that step touches — never
`git add -A`, given the dirty tree. No AI/tool mentions anywhere.

| # | Contents | Baseline regen |
|---|---|---|
| 0 | step 0 — the pre-existing work, message approved first | no |
| 1 | step 1 + step 2 + their test block | no |
| 2 | step 3 + its test block | no (commit body records the duplicated-`samples=` change) |
| 3 | step 4 + its test block + the `test-graceful-errors.R` test | no |
| 4 | step 6 — comments and `dev/ARCHITECTURE.md` | no |
| 5 | step 7 — `vignettes/ondisk_stlists.Rmd` + `.Rbuildignore` | no |
| 6 | step 8 — the rebuilt `docs/` | no |

Step 1 is one commit because committing `sample_slot_counts()` alone would add
an unused function.

## Verification

Per fix commit:

```sh
Rscript -e 'devtools::load_all(); testthat::test_local()'
```

Before the vignette, and again at the end:

```sh
SPATIALGE_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_local()'   # expect zero diffs
Rscript dev/ci/check_vignette_args.R                                             # statically checks the new .Rmd
```

`check_vignette_args.R` parses every `vignettes/*.Rmd` on every push and
validates named arguments against formals without evaluating anything, so the
new article is covered by CI from the moment it lands.

End-to-end, the thing the article actually claims — run before writing it, so
the prose describes observed behavior rather than intended behavior:

```r
devtools::load_all()
st <- STlist(rnacounts=visium_folders, samples=clin_file, cores=2)
st <- filter_data(st, spot_minreads=5000, spot_mingenes=1000, spot_maxreads=150000)
st <- transform_data(st, method='log', cores=2)
bk <- STlist_to_disk(st, file.path(tempdir(), 'tnbc_store'), overwrite=TRUE)

bk                                                    # real dimensions, not blanks
identical(summarize_STlist(bk), summarize_STlist(st))
validate_store(bk, deep=TRUE)
STplot(bk, genes='NDRG1', samples='sample_094c')      # errors today
distribution_plots(bk, genes='NDRG1')                 # silently empty today
identical(STlist_to_memory(bk)@counts, st@counts)
```

Then render the article alone before committing to a full site build:

```sh
Rscript -e 'pkgdown::build_article("ondisk_stlists")'
```

## Found but out of scope

Not touching these; listing them because they came out of the same read and are
each worth their own decision:

- **Duplicate gene names silently truncate a backed gene-subset read.**
  `sge_store_materialize()` resolves rows with `match()`, which keeps only the
  first position per name, while the in-memory mask keeps every duplicate row.
  `STlist_to_disk()` records `gene_dup` but does not refuse; only
  `validate_store()` raises (`error0067`). Reachable on real 10x data —
  duplicated symbols in `features.tsv` are routine, and
  `filter_data_stages.R` records that duplicates now pass through silently. The
  fix is to raise `error0067` at write time.
- **`distribution_plots(genes=)` matching nothing returns broken plots** rather
  than raising `error0090`; the ggplot2 internal `object 'samplename' not found`
  only appears at print time. Reachable in memory today, not specific to the
  store.
- **`plot_spatial_expression()` has no `assert_transformed()`** while
  `plot_spatial_geneset()` does, so the two `STplot` paths give different
  messages on an untransformed STlist.
- **`setValidity()` skips a backed STlist** — its guard requires a matrix, so a
  manifest bypasses the gene/spot-name and `spatial_meta`-row-count checks that
  `x[i]` would otherwise catch.
- **`@misc$store` has no branch in `bl_compare_stlist()`**, so adding a backed
  step to `bl_pipeline()` would trip the `misc_uncompared_keys` guard. It embeds
  tempdir paths and mtimes, so it needs a normalizer rather than a comparison.
