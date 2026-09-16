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
- Do not run `devtools::document()`. Installed roxygen2 is newer than
  `RoxygenNote` in `DESCRIPTION`; regenerating would churn all of `man/` and
  `NAMESPACE` and bury the real diff. If an `.Rd` genuinely needs to change,
  hand-edit it and say so in the commit message.

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
