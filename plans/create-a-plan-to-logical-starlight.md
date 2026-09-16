# Plan: isolated end-to-end validation of spatialGE (`alex-refactor`)

## Context

`alex-refactor` (HEAD `1b2663e`) is a large in-progress refactor: 71 files
changed vs `main`, +7290/−3342 lines, while carrying the **same
`Version: 1.2.2`** as `main`. All development so far has gone through
`devtools::load_all()` — **the package has never actually been installed** into
the `spatialge-dev` conda env (`packageVersion("spatialGE")` → MISSING). And it
has never been run against the real example data: the test suite uses only the
synthetic fixtures from `inst/testdata/make_fixtures.R`, and all 14 `@examples`
blocks are `\donttest{}`-wrapped *and* internally `tryCatch`ed, so they have
effectively never executed.

This plan validates the branch as a **package** rather than as a source tree:
install it, feed it the real `FridleyLab/spatialGE_Data` example data, run the
suite and every example, and diff its outputs against `main`.

**Nothing in the repository is modified. No git write operations of any kind**
(no commit/push/pull/fetch/checkout/stash/add). Every byte produced lands in
one folder that the user deletes afterward:

```
$T = /Users/4471592/Documents/share/lab_soupir/projects/spatialGE/DELETE_ME_spatialge_testrun_2026-09-16
$REPO = /Users/4471592/Documents/share/lab_soupir/projects/spatialGE
```

### Scope (confirmed with user)

Depth **Standard** — install, test suite (plain + golden), all 14 examples on
real data, light `R CMD check`; vignette knits deferred to an optional
appendix. Head-vs-`main` **output diff: yes**. ggplot2 pin: **leave the env
alone**. Temp folder at the repo root.

---

## Verified findings that shape this plan

Every item below was checked directly, and four of them overturned my initial
assumptions.

**F1 — BLOCKER: conda R cannot do HTTPS in this shell.**
```
$ /opt/anaconda3/envs/spatialge-dev/bin/Rscript -e 'curlGetHeaders("https://github.com/...")'
libcurl error code 60: self-signed certificate in certificate chain (19)
```
The TLS chain is corporate-intercepted (`O=Moffitt Cancer Center,
CN=ca.moffitt.goskope.com`). `/usr/bin/curl` works (macOS keychain trust);
conda's libcurl does not, and `CURL_CA_BUNDLE=/etc/ssl/cert.pem` alone does not
fix it. **Every `download.file()` in tests, examples, and vignettes fails until
a combined CA bundle exists.** `dev/environment/README.md` documents this fix;
it has not been applied in this shell. This gates every network tier.

**F2 — `tests/golden/` is NOT stale, so the golden tier should PASS.**
Golden files: `Sep 16 11:18:17`. `R/stlist_process.R` (the output-affecting
`f278bd0` sort fix): `11:17:32` — *before* the capture, and that commit's
message says it regenerated the baseline. Two files *are* newer
(`R/STList.R` 11:22, `R/utils.R` 11:32 — commits `6954044`/`ede00c8`), but both
are documented output-neutral (a cores-resolution tidy, plus an `HLfit`-class
guard on a never-taken branch and a `count_cores()` clamp active only under
`_R_CHECK_LIMIT_CORES_`). **Reversal of my earlier assumption: a golden failure
is a genuine finding, not expected staleness.**

**F3 — Do NOT install ggplot2 3.5.1; the pin is stale, not the env.**
ggplot2 4.0.3 was built `Sep 15 20:49:58`, *before* the `Sep 16 11:18` golden
capture — so **golden was captured under 4.0.3** and comparing under 4.0.3 is
apples-to-apples. Worse, `ggrepel` 0.9.8 declares `Depends: ggplot2 (>= 3.5.2)`,
so shadowing 3.5.1 in a temp lib would hard-break package loading. **Reversal:
plot-valued golden diffs are NOT expected**, and `dev/environment/install_deps.R:19`
is a doc-staleness finding to report, not act on. Also never *run*
`install_deps.R` — it does `sink("dev/environment/RECORDED.md")`, a repo write.

**F4 — Nothing needs installing.** All 44 `Imports`, `LinkingTo`, and every
`Suggests` are present, plus `rcmdcheck`, `pkgbuild`, `pkgload`, `remotes`,
`covr`, `msigdbr`. The temp libraries hold **only** spatialGE. (Drop my
"lib-deps" idea.)

**F5 — All 14 examples are broken, and it's pre-existing on `main`.**
`grep -l "#' zip_tmp" man/*.Rd` → **14 of 14**; `git show main:man/STlist.Rd`
has it too. The mis-escaped `#'   #' zip_tmp = list.files(...)` comments out the
assignment, so every example dies at `unzip(zipfile=zip_tmp)` with `object
'zip_tmp' not found` — swallowed into `message("Could not run example. Are you
connected to the internet?")`, which actively lies about the cause. **Reversal:
this is 14 examples, not just `R/STList.R:112`, and it is not a refactor
regression.**

**F6 — `git archive main` carries committed `src/*.o` and `src/spatialGE.so`**
(removed on `alex-refactor`). `git archive` gives all files identical mtimes, so
`make` may link stale objects from another toolchain. Must `rm -f` them after
extraction. Otherwise `git archive` is the right isolation mechanism:
tracked-only (skips the repo's own untracked `src/*.o`, `tests/golden/`,
`.claude/`), no working-tree touch, no `checkout`/`stash`.

**F7 — Four "escapes the temp folder" traps.**
- `R_LIBS_USER` = `~/Library/R/arm64/4.5/library`, which **does not exist**, so
  `.libPaths()` is a single entry and that entry **is** the shared conda
  library. A bare `install.packages()` writes straight into the shared env.
- `msigdbr` caches to `~/Library/Caches/org.R-project.R/R/msigdbr` → set
  `R_USER_CACHE_DIR`.
- `~/.R/Makevars` **exists** (system `/usr/bin/clang` + SDKROOT) and is R's
  default `R_MAKEVARS_USER`; `~/.Renviron` **exists** and exports
  `CC`/`CXX=/usr/bin/clang`. So the conda toolchain is *not* the default — set
  `R_MAKEVARS_USER` explicitly. `dev/environment/Makevars` uses `CONDA_PREFIX ?=`,
  so also export `CONDA_PREFIX` explicitly or it inherits `base`.
- `R CMD build` writes the tarball to cwd, and any base-graphics `plot()` drops
  `Rplots.pdf` in cwd → **never start R with cwd inside the repo.**

**F8 — `/usr/local/bin/R` is R 4.5.0 and is first on `PATH`.** Never invoke
bare `R`/`Rscript`; use absolute conda paths.

**F9 — Data sizes:** `melanoma_thrane.zip` 3.4 MB, `tnbc_bassiouni.zip`
**114 MB**, `nanostring_nsclc.zip` **51 MB**. `R CMD build` builds vignettes by
default, which alone pulls 114 MB and runs an 8-sample Visium pipeline → build
with `--no-build-vignettes`, check with `--no-vignettes`.

**F10 — `main` has essentially no harness**: its `tests/` is `testthat.R` +
`test-tests.R` only; `inst/testdata/make_fixtures.R`, `dev/baseline/*`,
`helper-fixtures.R`, `R/utils_parallel.R` are all absent, and its `DESCRIPTION`
lacks `waldo`/`writexl`. But `NAMESPACE` and all of `man/` are **byte-identical**
across the two branches, so HEAD's `dev/baseline/pipeline.R` is call-compatible
with `main` — the output-diff tier is viable by driving `main` with **HEAD's**
harness, sourced by absolute path.

**F11 — `main` has `Copyright: ins/COPYRIGHT`** (typo; `inst/` on HEAD) → the
`main`-side check gets one extra WARNING that HEAD does not. Useful as a sanity
check that the two checks aren't confused for each other.

**F12 — `pdflatex` present, `qpdf` absent** → manual builds; PDF-size NOTE only.

**F13 — `unlink(tempdir())` is contained but still disruptive.** `tempdir()` is
the per-session `$TMPDIR/RtmpXXXX`, not `$TMPDIR` — so it can't wipe sibling
dirs. But `test-tests.R` and all 14 examples call it, destroying materialized
fixtures mid-run. Keep downloaded data outside `TMPDIR`, and don't reorder the
suite with `filter=` and then blame the package.

**F14 — reportable harness gap:** `dev/baseline/pipeline.R:22-28` defines
`BL_NONCOMPARABLE` (`gene_interpolation`, `STplot_interpolation`, `STdiff`,
`STdiff_volcano`, `STenrich_gsva`) but `bl_compare()` has **no branch that
honors it** — those names are compared like any other despite being documented
as optimizer/KDE-nondeterministic.

---

## Layout

```
$T/
  .START_MARKER      # touched first; anchors the "nothing escaped" proof
  env.sh             # single sourceable preamble
  certs/             # netskope.pem, combined.pem  (F1)
  src-head/          # git archive HEAD  (+ seeded tests/golden/)
  src-head-recap/    # 2nd HEAD copy, for the re-capture only
  src-main/          # git archive main  (src/*.o and *.so deleted!)  (F6)
  lib-head/ lib-main/  # temp R libraries — spatialGE only  (F4)
  data/              # downloaded zips + unzipped trees (NOT under rtmp)
  rtmp/              # TMPDIR
  cache/             # R_USER_CACHE_DIR  (F7)
  check/{head,main}/ # R CMD check -o  (separate: both tarballs are 1.2.2)
  examples/{raw,run}/
  out/               # tarballs, result .rds, reports
  logs/              # one log per step
```

### `$T/env.sh`

```sh
export REPO=/Users/4471592/Documents/share/lab_soupir/projects/spatialGE
export T="$REPO/DELETE_ME_spatialge_testrun_2026-09-16"
export CONDA_PREFIX=/opt/anaconda3/envs/spatialge-dev          # F7
export PATH="$CONDA_PREFIX/bin:/Library/TeX/texbin:$PATH"
export RBIN="$CONDA_PREFIX/bin/R"; export RSCRIPT="$CONDA_PREFIX/bin/Rscript"  # F8
# containment
export TMPDIR="$T/rtmp"
export R_USER_CACHE_DIR="$T/cache"                             # F7
export R_LIBS_USER="$T/lib-head"; export R_LIBS="$T/lib-head"  # F7
export R_MAKEVARS_USER="$T/src-head/dev/environment/Makevars"   # F7
# TLS
export CURL_CA_BUNDLE="$T/certs/combined.pem"                  # F1
export SSL_CERT_FILE="$T/certs/combined.pem"
unset SPATIALGE_GOLDEN SPATIALGE_CORES _R_CHECK_LIMIT_CORES_
```

---

## Tiers

### T0 — Setup (~8 min, mostly the 170 MB download)

1. `mkdir -p` the layout; `touch "$T/.START_MARKER"`.
2. **Provenance manifests before anything runs** (read repo, write only to `$T`):
   `git status --porcelain`, `git status --porcelain --ignored`,
   `git rev-parse HEAD`, md5 **and** mtime of every `tests/golden/*.rds`, and
   `ls -1` of the shared conda library.
3. **CA bundle (F1)** — read-only on the keychain, writes only into `$T`:
   ```sh
   security find-certificate -a -c ca.moffitt.goskope.com -p /Library/Keychains/System.keychain > "$T/certs/netskope.pem"
   cat /etc/ssl/cert.pem "$T/certs/netskope.pem" > "$T/certs/combined.pem"
   ```
   **Gate:** `"$RSCRIPT" -e 'cat(curlGetHeaders("https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download=")[1])'`
   must print `HTTP/2 302`, not `libcurl error code 60`. If it still fails,
   append the login-keychain `goskope` certs and rebuild. Do not proceed to any
   network tier until this passes.
4. **Isolated sources:**
   ```sh
   git -C "$REPO" archive HEAD | tar -x -C "$T/src-head"
   git -C "$REPO" archive HEAD | tar -x -C "$T/src-head-recap"
   git -C "$REPO" archive main | tar -x -C "$T/src-main"
   rm -f "$T/src-main/src/"*.o "$T/src-main/src/"*.so          # F6 — mandatory
   mkdir -p "$T/src-head/tests/golden" && cp -Rp "$REPO/tests/golden/." "$T/src-head/tests/golden/"
   ```
   (52 files expected. The repo's copy is only ever read.)
5. **Example data** via `/usr/bin/curl` (keychain-trusted, works regardless of
   step 3) for all three zips, then `unzip -d "$T/data"`. Guard: any zip under
   100 KB is an unfollowed Git-LFS pointer (`git-lfs` is not installed locally,
   so `git clone` would yield only pointers) — stop and fix the URL.

### T1 — Install HEAD into the temp library (~3 min)

```sh
cd "$T/out" && "$RBIN" CMD INSTALL --library="$T/lib-head" --preclean --no-multiarch "$T/src-head"
"$RSCRIPT" -e 'library(spatialGE); print(packageVersion("spatialGE")); print(find.package("spatialGE"))'
```
**Pass:** exit 0, `DONE (spatialGE)`, all three `.cpp` compile,
`find.package()` → `$T/lib-head/spatialGE`.
**Watch:** if the log shows `/usr/bin/clang` rather than `$CONDA_PREFIX/bin/clang`,
`R_MAKEVARS_USER`/`CONDA_PREFIX` didn't take (F7).

### T2 — Plain test suite (~5 min)

```sh
cd "$T/src-head" && "$RSCRIPT" -e 'testthat::test_local(".", reporter="summary")'
```
Plus one run with `load_package="installed"` — that variant is what actually
proves `inst/testdata/` shipped, since `helper-fixtures.R` resolves
`system.file("testdata", ..., package="spatialGE")` through `$T/lib-head`.
`test_local()` uses `pkgload::load_all()`, so compiled artifacts land in
`$T/src-head/src/`, never in `$REPO/src/`.

**Expect:** 36 blocks — 11 `test-fixtures` + 1 `test-filter_data` + 20
`test-utils_validate` + 1 `test-tests` (network) + 3 `test-golden` (**SKIP**,
"golden comparison disabled"). **Pass = FAIL 0, SKIP 3.**
A `test-tests.R` skip reading "Could not download Thrane et al. test data"
means T0.3 didn't stick (F1) — an environment problem, not a package one.

### T3 — Golden regression tier (~12 min) — **run alone**

```sh
cd "$T/src-head" && SPATIALGE_GOLDEN=true "$RSCRIPT" -e 'testthat::test_local(".", reporter="summary")'
```
Leave `SPATIALGE_CORES` and `_R_CHECK_LIMIT_CORES_` unset — the 11:18 capture
ran with both unset and fidelity rule 12 requires both sides identical.

**Expect PASS** (F2/F3). Triage order on failure:
1. Which slot does the `waldo` diff name?
2. One of the five `BL_NONCOMPARABLE` names → known nondeterminism, *and*
   evidence for finding F14.
3. `@counts` / `@spatial_meta` row order → real regression in
   `process_lists()` territory.
4. A plot-valued result → unexpected per F3; investigate rather than
   dismiss as ggplot2 drift.

**T3b — re-capture in the *second* copy** so the seeded baseline stays pristine
(`capture_baseline.R` hardcodes `golden_dir <- "tests/golden"` relative to cwd):
```sh
cd "$T/src-head-recap" && "$RSCRIPT" dev/baseline/capture_baseline.R
```
then `bl_compare()` all 52 file pairs (`src-head/tests/golden` vs
`src-head-recap/tests/golden`) into `$T/logs/T3b_recapture_diff.log`. This is
also what makes the `cores1`-vs-`cores3` invariance test meaningful, since that
test compares *captured files*, not fresh runs.

### T4 — `R CMD build` + `R CMD check --run-donttest` (~20 min)

```sh
cd "$T/out"
"$RBIN" CMD build --no-build-vignettes "$T/src-head"                      # F9
tar -tzf "$T/out/spatialGE_1.2.2.tar.gz" | grep -E "^spatialGE/(dev|plans|tests/golden)/" \
  || echo "correctly excluded"
"$RBIN" CMD check --run-donttest --no-vignettes --no-multiarch \
  -l "$T/lib-head" -o "$T/check/head" "$T/out/spatialGE_1.2.2.tar.gz"
```
`-o` is mandatory or `spatialGE.Rcheck/` lands in cwd. Tests stay **enabled**:
`R CMD check` sets `_R_CHECK_LIMIT_CORES_`, so this is the only tier that
exercises the `count_cores()` clamp from `ede00c8`.

**T4b — same check on `main`**, purely to establish the pre-existing NOTE
baseline. Both tarballs are named `spatialGE_1.2.2.tar.gz`, so rename main's and
use `-o "$T/check/main"` and `-l "$T/lib-main"`. Then
`diff` the two `grep -E "^\* checking .*\.\.\. (NOTE|WARNING|ERROR)"` lists —
**anything on HEAD but not on `main` is new.** (Expect exactly one main-only
WARNING for the `ins/COPYRIGHT` typo, per F11, confirming the two runs weren't
crossed.)

**Pre-existing, expected:** `LazyData` specified without a `data/` directory;
`unable to run 'qpdf'` / PDF size NOTE; `Namespaces in Imports field not
imported from:` (44 `Imports` vs a thin `NAMESPACE`); installed-size NOTE.
**`checking examples ... OK` is a false positive** — the examples swallow their
own failures (F5). Examples are validated in T5, not here.

### T5 — Loud example driver (~25 min) — the tier that finds things

1. **Extract:** `tools::Rd2ex(rd, out, commentDontrun=TRUE, commentDonttest=FALSE)`
   over all `$T/src-head/man/*.Rd` → 14 scripts in `$T/examples/raw/`.
   `Rd2ex` writes only to the given `out` path.
2. **Generate two variants per example** in `$T/examples/run/` (patching copies
   only, never `$REPO`):
   - `__asis` — replace only the `download.file(lk, ...)` line with a
     `file.copy()` from `$T/data/melanoma_thrane.zip`, so the *only* possible
     failure is the package's own defect, not the network.
   - `__fixed` — additionally un-comment the mis-escaped `#' zip_tmp = list.files`
     assignment (F5).
3. **Run each in a fresh `Rscript --vanilla`** with `.libPaths("$T/lib-head")`
   and a narrow `globalenv()` shim that turns the swallow into a hard error —
   examples eval in `globalenv()`, so an unqualified name defined there wins:
   ```r
   message <- function(..., domain=NULL, appendLF=TRUE) {
     txt <- paste0(unlist(list(...)), collapse="")
     if (grepl("Could not run example|Could not download", txt))
       stop("SWALLOWED EXAMPLE FAILURE: ", txt, call.=FALSE)
     base::message(..., domain=domain, appendLF=appendLF)
   }
   ```
   Then a post-hoc assertion (`is(melanoma,"STlist")`, `length(@counts) > 0`)
   and `cat("EXAMPLE_OK")`. One process per example, so a failure can't poison
   the next; recreate `$TMPDIR` between runs (F13).
4. `$T/out/report_examples.md`: 28 rows → exit code, sentinel hit, verdict.

**Expect:** all 14 `__asis` **fail** with `object 'zip_tmp' not found`
(pre-existing, reproduces on `main`); all 14 `__fixed` **pass** with
`EXAMPLE_OK`. **A `__fixed` failure is a genuine, new finding** in that
function — report with traceback. An `__asis` *pass* means the patch didn't
apply.

`devtools::run_examples(".", run_donttest=TRUE, document=FALSE)` is an optional
smoke cross-check only — it reports OK because of the swallow, and
`document=FALSE` is mandatory or roxygen rewrites `man/`.

### T6 — HEAD vs `main` output diff (~20 min)

Drive both sides with **HEAD's** `make_fixtures.R` + `pipeline.R` by absolute
path, since `main` ships neither (F10):

```sh
"$RBIN" CMD INSTALL --library="$T/lib-main" --preclean --no-multiarch "$T/src-main"
# per side: .libPaths($T/lib-<side>); library(spatialGE);
#   source($T/src-head/inst/testdata/make_fixtures.R); source($T/src-head/dev/baseline/pipeline.R)
#   saveRDS(bl_pipeline(cores=1, verbose=FALSE), $T/out/pipeline_<side>_cores1.rds, version=3)
# then bl_compare() the 26 names pairwise -> $T/out/report_versions.md
```
`bl_capture()` already does `set.seed(20260115)` per call and `tryCatch`es
errors into the result, so a signature change on `main` is recorded as a
comparable error instead of aborting the run.

**Framing: a review artifact, not a pass/fail gate.** `STlist` and everything
downstream is *expected* to differ — `f278bd0` deliberately changed spot/cell
ordering, so a `DIFF` there **confirms the fix landed**. The five
`BL_NONCOMPARABLE` names differ for optimizer reasons. Also expect
error-message differences from the removed `inst/err.csv` rows. Anything else
that differs must be traced to a specific commit in `git log main..HEAD`.

---

## Verdict table — pre-existing vs new

| Symptom | Verdict |
|---|---|
| `libcurl error code 60` in any R network call | **Environment** (F1); fixed in T0.3 — if it persists, every network tier is invalid, not failing |
| All 14 `__asis` examples die on `object 'zip_tmp' not found` | **PRE-EXISTING**, `main` too (F5) |
| `R CMD check`: `checking examples ... OK` | **False positive** — trust T5 only |
| `LazyData` without `data/`; qpdf/PDF NOTE; `Namespaces in Imports` NOTE; installed-size NOTE | PRE-EXISTING |
| `Copyright field refers to nonexistent ins/COPYRIGHT` | PRE-EXISTING, **`main` only** (F11) |
| T3 golden diff on a `BL_NONCOMPARABLE` name | Known nondeterminism + evidence for F14 |
| `ggplot2` pin 3.5.1 in `install_deps.R:19` | Stale script (F3) — report, do not act |
| **T3 golden diff on `STlist`/`filter_data`/`transform_data_*`/`SThet`/`STclust`** | **NEW — real regression** (F2) |
| **Any T2 offline test failure** | **NEW** |
| **Any `__fixed` example failure in T5** | **NEW** |
| **Any HEAD-only NOTE/WARNING/ERROR in the T4-vs-T4b diff** | **NEW** |
| T6 `STlist` + downstream `DIFF` vs `main` | **Expected** — `f278bd0`'s intentional ordering fix |

Deliverable: `$T/out/SUMMARY.md` — verdict per tier plus the three sub-reports
(`report_golden`, `report_examples`, `report_versions`).

## Timing

Serial ≈ **1 h 15 m**. T3 must run alone (fork/timing sensitivity, fidelity
rule 12). T4/T4b are sequential w.r.t. each other only if `-o` collides — use
separate `-o` dirs and they can overlap with T5. 12 physical cores and each
call auto-detects `detectCores()/2 = 6` forks, so cap at **2 concurrent
tiers**. Do **not** use `SPATIALGE_CORES` to bound concurrency — it changes the
code path at every call site and would invalidate T3 and T6.

## Cleanup verification (before the user deletes `$T`)

1. **Tracked tree:** `git status --porcelain` diffs clean vs the T0 manifest;
   `git diff --stat HEAD` empty; `rev-parse HEAD` unchanged; `stash list`
   empty; `branch --show-current` still `alex-refactor`.
2. **Ignored set:** `git status --porcelain --ignored` diff shows *only* the new
   `DELETE_ME_…` line — proving `src/*.o`, `src/spatialGE.so`,
   `dev/profiling/out/` untouched.
3. **`tests/golden/` byte- AND mtime-identical** (md5 + `stat -f '%m %N'` vs the
   T0 manifests). Both must diff empty.
4. **Nothing else written in the repo:**
   `find "$REPO" -newer "$T/.START_MARKER" -not -path "$T/*" -not -path "$REPO/.git/*" -not -path "$REPO/.claude/*" -not -path "$REPO/plans/*"`
   → expect nothing but `$REPO` itself.
5. **Shared conda library untouched:** `ls -1` diff empty;
   `[ ! -d $CONDA_PREFIX/lib/R/library/spatialGE ]`.
6. **No `$HOME` writes:** nothing newer than the marker under `~/Library/R`,
   `~/Library/Caches/org.R-project.R`, `~/.R`; no new `~/.Rhistory`/`~/.RData`.
7. **System `TMPDIR` untouched:** no new `Rtmp*` under `/var/folders/.../T`.
8. **No strays at the repo root:** no `Rplots.pdf`, `*.tar.gz`, `*.Rcheck`.

Then the **user** runs: `rm -rf "$REPO/DELETE_ME_spatialge_testrun_2026-09-16"`.
I will not delete it.

## Out of scope

- **Vignette knits** (Standard depth). If revisited: knit *copies* in `$T`
  (knitr writes figures next to the `.Rmd`); `basic_functions_vignette.Rmd`
  gates ~30 chunks on `eval=net_check`, so **a failed download renders a green,
  empty document** — verification must grep the produced HTML for real output,
  not trust exit 0. `tnbc_bassiouni` is 114 MB; the SMI vignette will likely
  fail in `msigdbr()` because the `msigdbdb` data package isn't installed
  (an env gap, not a spatialGE bug); the STdiff vignette is the single slowest
  thing available (`bench_units.csv` measures `stdiff_spatial_one_cluster` at
  20.8 s/cluster).
- **Fixing anything.** At least three real defects will be documented and left
  alone: the 14 broken examples (F5), the stale ggplot2 pin (F3), and the
  unhonored `BL_NONCOMPARABLE` list (F14).
- `--as-cran`, `dev/profiling/`, `covr` coverage, and every `git` write.
