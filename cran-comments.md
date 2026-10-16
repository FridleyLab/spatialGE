# spatialGE 2.0.0

## Test environments

* macOS (local), R 4.5.3
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest, R 4.5.3

## R CMD check results

0 ERRORs, 0 WARNINGs, 3 NOTEs. All three are properties of the machine that ran
the check rather than of the package, and the reasoning is given for each so a
reviewer does not have to reconstruct it.

Four findings present in earlier drafts of this file have been fixed rather than
explained, and are recorded here only because a reviewer comparing submissions
may look for them:

* *Undocumented S4 methods* (`[`, `[[`, `length` on `STlist`) -- these are now
  documented, along with the `@param` entries for the generics' full argument
  lists.
* *Namespace in Imports field not imported from: 'R.utils'* -- `R.utils` is a
  **runtime** dependency of `data.table::fread()`, which loads it to decompress
  the gzipped `features`/`barcodes`/`matrix` files Cell Ranger always emits.
  `fread()` loads it internally rather than spatialGE calling it, so static
  analysis could not see the use. It is now named in an `@importFrom` beside
  those `fread()` calls. It is not moved to Suggests: a minimal install then
  cannot read standard 10x MEX output at all.
* *'qpdf' is needed for checks on size reduction of PDFs* -- `qpdf` is now
  installed on the machine that runs the check.
* *Possibly invalid URL `https://clustrmaps.com/site/1bwt8`* -- a visitor-counter
  badge, now commented out in README.md.

### NOTE: CRAN incoming feasibility

Two parts. The maintainer line is expected.

The URL reported as possibly invalid is
`https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/what-is-space-ranger`,
cited in the Space Ranger reference of `basic_functions_vignette.Rmd`. 10x
Genomics appears to have retired the `support.10xgenomics.com` host: it returns
403 to any client, including a browser user agent, and issues no redirect. The
current location is `https://www.10xgenomics.com/support/software/space-ranger`,
which could not be confirmed from the submitting network because 10x rate-limits
it (HTTP 429).

The incoming-feasibility check itself could not run locally, because the network
blocks the CRAN and Bioconductor package indexes.

### NOTE: checking for future file timestamps -- unable to verify current time

Environmental. The check contacts an external time service to compare file
timestamps, and the submitting network blocks it. No file in the tarball carries
a future timestamp.

### NOTE: checking HTML version of manual

Environmental: `HTML Tidy` on the machine that ran the check is older than the
version R requires for HTML validation, so that validation is skipped rather
than failed.

## Downstream dependencies

None known.

## Notes on this release

2.0.0 contains behavior changes that are documented in NEWS.md. Two are worth
drawing attention to because existing results may be affected:

* `STenrich()` previously accepted `annot` without `domain` and silently ran the
  **unrestricted** analysis, returning results a caller would reasonably read as
  domain-restricted. This is now an error. Results produced that way in 1.2.x are
  not domain-restricted.

* Requesting a gene, annotation or reference domain that matches nothing is now an
  error rather than a silent no-op, so a typo no longer looks like a run that
  legitimately found nothing.
