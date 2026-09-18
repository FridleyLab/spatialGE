# spatialGE 2.0.0

## Test environments

* macOS (local), R 4.5.3
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest, R 4.5.3

## R CMD check results

0 ERRORs. The NOTEs below are expected; the reasoning for each is given so a
reviewer does not have to reconstruct it.

### NOTE: Namespace in Imports field not imported from: 'R.utils'

This is a false positive and the dependency is load-bearing.

`R.utils` is a **runtime** dependency of `data.table::fread()`, which requires it
to decompress `.gz` input read directly from a path. `import_visium()` and
`import_xenium()` call `fread()` on the gzipped `features`/`barcodes`/`matrix`
files that Cell Ranger and Xenium Ranger always produce, so with `R.utils` only in
Suggests a minimal install cannot read standard 10x MEX output at all.

Because `fread()` loads it internally rather than spatialGE calling it by name,
static analysis cannot see the use. Moving it to Suggests would silence the NOTE
and reintroduce the failure.

### NOTE: CRAN incoming feasibility / possibly invalid URL

`https://clustrmaps.com/site/1bwt8` in README.md is a visitor-counter badge and
currently returns an error. It predates this release. It is not a documentation
or code link and can be removed if preferred.

The remainder of this NOTE is the maintainer line plus the incoming-feasibility
check, which could not run in the local environment because the network blocks
access to the CRAN and Bioconductor package indexes.

### WARNING: 'qpdf' is needed for checks on size reduction of PDFs

Environmental, not a package property: `qpdf` is not installed on the machine that
ran the check. It is present on CRAN's check machines, so this does not reproduce
there.

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
