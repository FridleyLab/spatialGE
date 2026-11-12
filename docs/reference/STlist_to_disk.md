# STlist_to_disk: Move an STlist's count matrices to on-disk parquet

Converts an in-memory STlist into one backed by parquet files, so data
larger than memory can be held and read in slices.

## Usage

``` r
STlist_to_disk(x = NULL, root = NULL, overwrite = FALSE)
```

## Arguments

- x:

  an STlist

- root:

  directory to write the store into. Created if absent

- overwrite:

  logical, whether to overwrite an existing store at `root`

## Value

an STlist whose count slots are on-disk backed

## Details

Opt-in and reversible. `@counts` and `@tr_counts` elements are replaced
by manifests describing per-sample parquet files; every other slot is
untouched. All package functions continue to work, because slot reads go
through internal accessors that materialize on demand – so this trades
memory for I/O rather than changing results.

Matrices are stored in COO (long) form, one file per sample per slot,
sorted by gene then cell.
[`STlist_to_memory()`](https://fridleylab.github.io/spatialGE/reference/STlist_to_memory.md)
reverses the conversion.

Note the store is a directory of files, not part of the R object. An
STlist saved with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
records paths, so moving or deleting `root` breaks it;
[`validate_store()`](https://fridleylab.github.io/spatialGE/reference/validate_store.md)
reports that rather than failing obscurely later.
