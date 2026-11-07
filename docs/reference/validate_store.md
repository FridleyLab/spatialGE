# validate_store: Check an on-disk STlist against its manifests

Verifies that a backed STlist's files exist and still match what its
manifests record, reading parquet footers only.

## Usage

``` r
validate_store(x = NULL, deep = FALSE)
```

## Arguments

- x:

  a backed STlist

- deep:

  logical, whether to read data pages as well as metadata

## Value

invisibly TRUE, or stops with the first problem found

## Details

Checks sample order against the recorded `sample_order` (fidelity rule
7), that every file exists and its size/mtime fingerprint is unchanged,
that dimnames lengths agree with the recorded dimensions, and that
`@tr_counts` genes and cells are subsets of `@counts`' in the same
relative order – the invariant
[`filter_data()`](https://fridleylab.github.io/spatialGE/reference/filter_data.md)
assumes. With `deep=TRUE`, also reads the data to check index bounds and
duplicate coordinates.
