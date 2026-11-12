# STniche (mclust version): Detect functional spatial domains with model-based clustering

Unsupervised clustering of spots/cells in one or more ST samples using
gene-set GSVA scores -\> local Moran's I (Z) as features, then Gaussian
mixture modeling via mclust with BIC-based model/cluster selection.

## Usage

``` r
STniche(
  x = NULL,
  gene_sets = NULL,
  samples = NULL,
  ks = 2:10,
  knns = 6,
  min_genes = 1,
  out_i = FALSE,
  verbose = 1L,
  cores = NULL
)
```

## Arguments

- x:

  an STlist with normalized expression data in `x@tr_counts[[sample]]`
  (genes x spots)

- gene_sets:

  a *named* list of gene sets (each an atomic character vector of gene
  symbols)

- samples:

  vector of sample names or indices to run (default: all)

- ks:

  numeric vector of candidate cluster numbers to evaluate with BIC
  (default: 2:10)

- knns:

  integer k for kNN used in Local Moran's I (default: 6)

- min_genes:

  minimum overlapping genes required to compute a GSVA score (default:
  1)

- out_i:

  logical; if TRUE return list(STlist, Z.Ii matrix), else STlist only

- verbose:

  integer verbosity level (0=quiet, 1=basic, 2=chatty)

- cores:

  integer number of workers for parallel scoring (Unix only; default
  chooses sensibly)

## Value

If `out_i=TRUE`: `list(updated_STlist, Zi_matrix)` where `Zi_matrix` has
rows = gene sets, columns = "sample&&libname" spots. Otherwise, the
updated STlist.

## Details

Workflow:

1.  For each sample, compute GSVA scores for provided gene sets.

2.  For each gene set within each sample, compute Local Moran's I
    Z-score over a kNN graph (k = `knns`).

3.  Concatenate all samples into a single feature matrix (spots x
    gene-sets, entries are Z.Ii).

4.  Run mclust::Mclust with candidate cluster counts `ks` (e.g., 2:10).
    Best model & G selected by BIC.

5.  Write cluster labels back into `x@spatial_meta[[sample]]` under
    columns `stniche_dsplFalse` and `STniche_DTC` for plotting
    compatibility.
