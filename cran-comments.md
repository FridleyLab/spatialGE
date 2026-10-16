# spatialGE 2.0.0

## Test environments

* macOS (local), R 4.5.3
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest, R 4.5.3

## R CMD check results

0 ERRORs, 0 WARNINGs, 1 NOTE.

### NOTE: checking for future file timestamps -- unable to verify current time

Environmental, and the only finding left. The check compares file timestamps
against an external time service, and both services R consults
(`worldclockapi.com` and `worldtimeapi.org`) were unreachable or erroring when
this check ran -- the first returned HTTP 503, the second did not respond. The
submitting network was otherwise working normally at the time, so this is an
outage of those services rather than a local block, and R has no fallback. No
file in the tarball carries a future timestamp.

### Findings fixed rather than explained

Earlier drafts of this file explained six other findings. All are now fixed, and
are listed here only because a reviewer comparing submissions may look for them:

* *Undocumented S4 methods* (`[`, `[[`, `length` on `STlist`) -- now documented,
  including the `@param` entries for the generics' full argument lists and an
  explicit alias for `[[`, whose registered signature loses its trailing `ANY`s.
* *Namespace in Imports field not imported from: 'R.utils'* -- `R.utils` is a
  **runtime** dependency of `data.table::fread()`, which loads it to decompress
  the gzipped `features`/`barcodes`/`matrix` files Cell Ranger always emits.
  `fread()` loads it internally rather than spatialGE calling it, so static
  analysis could not see the use; it is now named in an `@importFrom` beside
  those `fread()` calls. Deliberately not moved to Suggests, which would silence
  the NOTE by making a minimal install unable to read standard 10x MEX output.
* *'qpdf' is needed for checks on size reduction of PDFs* -- `qpdf` installed on
  the machine that runs the check.
* *checking HTML version of manual* -- the machine had Apple's 2006 build of
  HTML Tidy, so R skipped HTML validation; HTML Tidy 5.8.0 is now installed.
* *Possibly invalid URL `https://clustrmaps.com/site/1bwt8`* -- a visitor-counter
  badge, now commented out in README.md.
* *Possibly invalid URL `https://support.10xgenomics.com/.../what-is-space-ranger`*
  -- 10x Genomics appears to have retired that host (403 to any client,
  including a browser user agent, with no redirect), so the vignette now refers
  to Space Ranger without linking to it.

Four `doi.org` citations and one `satijalab.org` page were also reported as
possibly invalid at one point. That was local TLS interception, not a problem
with the links: the check reported `libcurl error code 60` together with
`(Status without verification: OK)`, i.e. it confirmed each URL resolves and
only certificate verification had failed. Pointing R's libcurl at a CA bundle
that includes the institution's root resolves it, and all five verify normally.

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
