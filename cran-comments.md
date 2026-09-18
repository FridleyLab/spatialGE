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

### NOTE: checking installed package size

The package ships compiled code (`src/`) for the spatial-autocorrelation and
permutation kernels.

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
