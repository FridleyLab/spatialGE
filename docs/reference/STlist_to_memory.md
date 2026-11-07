# STlist_to_memory: Bring an on-disk STlist back into memory

Reverses
[`STlist_to_disk()`](https://fridleylab.github.io/spatialGE/reference/STlist_to_disk.md),
materializing every backed matrix.

## Usage

``` r
STlist_to_memory(x = NULL)
```

## Arguments

- x:

  an STlist, backed or not

## Value

an STlist with in-memory count matrices
