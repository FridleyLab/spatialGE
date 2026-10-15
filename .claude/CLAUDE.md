# Working in this repo

@../dev/ARCHITECTURE.md

## Environment

Use the `spatialge-dev` conda environment (R 4.5.3), not system R:

```sh
source /opt/anaconda3/etc/profile.d/conda.sh && conda activate spatialge-dev
export R_MAKEVARS_USER="$(pwd)/dev/environment/Makevars"
```

See `dev/environment/README.md` for env creation, dependency install, and
known gaps (e.g. `detectCores()` returning `NA` in sandboxed shells).

`renv/`, `renv.lock`, `.Rprofile` are gitignored and not the reproducibility
story here -- the conda yaml plus `dev/environment/install_deps.R` is.

### TLS: conda R cannot download anything until you fix the CA bundle

This network intercepts TLS with Netskope, and conda's libcurl does not trust
the intercepting root even though macOS (and therefore `/usr/bin/curl`) does. So
**every `download.file()` from conda R fails** with
`libcurl error code 60: self-signed certificate in certificate chain` until you
export a combined CA bundle. Full recipe, concrete CN, and the reason a green run
can be misleading: see the "Corporate network / TLS interception" section of
`dev/environment/README.md`.

Always run this gate before trusting anything network-dependent (tests,
examples, vignettes) -- expect `HTTP/2 302`, not `libcurl error code 60`:

```sh
Rscript -e 'cat(curlGetHeaders("https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=")[1])'
```

## Branch and commit discipline

- All work happens on `alex-refactor`. Never push -- pushes happen only when
  explicitly requested.
- Commit when a meaningful, self-contained change is complete. No mention of
  AI tools anywhere in commits, comments, or file contents.
- An output-affecting change gets its own commit that both makes the change
  and regenerates `tests/golden/` (via `dev/baseline/capture_baseline.R`) in
  the same commit -- see the "Regression / baseline workflow" section of
  `dev/ARCHITECTURE.md`.
- `devtools::document()` is now safe to run, and is the way to change an `.Rd`.
  It previously was not: several `.Rd` files had been hand-edited beyond their
  roxygen sources, so regenerating *deleted* documented behavior -- the
  `@details` and `@seealso` of `STlist_to_disk`, `STlist_to_memory`, `is_backed`,
  `validate_store` and `as_Seurat`. Those paragraphs now live in the roxygen
  blocks, so a regeneration round-trips. Do not hand-edit an `.Rd` any more:
  that is what created the divergence.

  Two things a `document()` run must not undo, both checked by CI:

  - **`importFrom(methods, .hasSlot)`.** `R/stlist_readers.R` calls it on a
    Seurat image slot. The call is conditional and `SeuratObject` is only in
    Suggests, so roxygen can drop the import with no other symptom. It is
    declared on `color_parse()` in `R/utils.R` for that reason.
  - **`S3method(base::print, sge_manifest)`.** A bare `@export` on
    `print.sge_manifest` makes roxygen2 >= 8 emit `export(print.sge_manifest)`
    instead, which leaves `print()` on a manifest undispatched unless the
    package happens to be attached. The tag is `@exportS3Method base::print`.

  roxygen2 8 also replaced `RoxygenNote` with `Config/roxygen2/version` in
  `DESCRIPTION`, and writes multi-symbol imports one per line. The second broke
  the first version of the `.hasSlot` CI guard, which grepped for the
  single-line form; it now asks `parseNamespaceFile()` instead.

## Failing gracefully (package-wide requirement)

**Every user-reachable failure must be an informative spatialGE error. An
internal error from `Matrix`, `dplyr`, `tibble`, `spdep`, `gstat` or base leaking
out to the caller is a bug, not an acceptable outcome.**

The bar: the message names the function and the offending argument, and tells the
caller what to do. These all failed it and were fixed:

| was | now |
|---|---|
| `object 'tr_results' not found` | `transform_data: Unrecognized method. Please use one of: log, sct.` |
| `attempt to select less than one element in get1index` | `STenrich: gene_sets must be a non-empty, named list of gene vectors.` |
| `object 'xpos' not found` | `gene_interpolation: None of the requested genes are present...` |
| `Can't extract columns past the end` | `SThet: Could not find normalized counts. Please use transform_data()...` |
| `knearneigh: fewer data points than k` | `SThet: k must be smaller than the number of spots/cells...` |

How to add one:

1. Append a row to `inst/err.csv` (`errorNNNN,'function: message.'`, next free
   code). Single-quoted, since `raise_err()` reads it with `quote="'"`.
2. Call `raise_err(err_code='errorNNNN')` as early as the condition is
   detectable -- validate at the entry point, not where the value is finally
   consumed. Most of the cases above failed many lines downstream of the bad
   input, which is why the messages named internals.
3. Add a case to `tests/testthat/test-graceful-errors.R`, including the specific
   leaked string in the `leaked` vector of the last test so the class cannot
   regress.
4. Prefer `assert_stlist()` / `assert_transformed()` / `resolve_samples()`
   (`R/utils_validate.R`) over a hand-rolled check where one fits.

Plain `stop('Please select a spatial weight between 0 and 1.')` is fine and does
not need converting -- it already names the argument and the fix. The rule is
about *internal* errors escaping, not about the `raise_err()` mechanism as such.

**Silently accepting invalid input also violates this**, and is worse than an
error because the caller believes it worked. The rule adopted is **error when
nothing matches, warn when only some does**:

| function | matched nothing | matched some |
|---|---|---|
| `SThet(genes=)` | `error0048` | message naming absent genes, proceed |
| `STgradient(annot=)` | `error0049` | message per sample, proceed on the rest |
| `STgradient(ref=)` | `error0050` | message per sample, proceed on the rest |
| `gene_interpolation(genes=)` | `error0046` | message naming absent genes, proceed |
| `STenrich(annot=`/`domain=)` | `error0051` | n/a -- both or neither |

`STenrich`'s case was different in kind and worse: domain-restricted mode tests
`!is.null(annot) & !is.null(domain)`, so supplying one without the other skipped
the entire block. The unrestricted analysis ran and returned results a caller
would reasonably read as domain-restricted, with nothing downstream able to
detect it because neither argument is consulted again.

## Verification

Per commit:

```sh
Rscript -e 'devtools::load_all(); testthat::test_local()'
```

Before merging, also run the golden gate (see `dev/ARCHITECTURE.md`):

```sh
SPATIALGE_GOLDEN=true Rscript -e 'devtools::load_all(); testthat::test_local()'
```

## Plans

Planning documents live in `plans/` (git-tracked, Rbuildignored), matching the
`spatialTIME` convention. The refactor this repo is mid-way through is
specified in `plans/this-is-an-r-velvety-music.md`.
