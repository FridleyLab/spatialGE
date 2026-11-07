# STdiff: Differential gene expression analysis for spatial transcriptomics data

Tests for differentially expressed genes using linear models with or
without spatial covariance structures

## Usage

``` r
STdiff(
  x = NULL,
  samples = NULL,
  annot = NULL,
  w = NULL,
  k = NULL,
  deepSplit = NULL,
  topgenes = 5000,
  pval_thr = 0.05,
  pval_adj = "fdr",
  test_type = "mm",
  sp_topgenes = 0.2,
  clusters = NULL,
  pairwise = FALSE,
  verbose = 1L,
  cores = NULL
)
```

## Arguments

- x:

  an STlist

- samples:

  an integer indicating the spatial samples to be included in the DE
  tests. Numbers follow the order in `names(x@counts)`. Sample names are
  also allowed. If NULL, performs tests on all samples

- annot:

  a column name in `x@spatial_meta` containing the groups/clusters to be
  tested. Required if `k` and `w` are empty.

- w:

  the spatial weight used in STclust. Required if `annot` is empty.

- k:

  the k value used in STclust, or `dtc` for dynamicTreeCut clusters.
  Required if `annot` is empty.

- deepSplit:

  the deepSplit value if used in STclust. Required if `k='dtc'`.

- topgenes:

  an integer indicating the top variable genes to select from each
  sample based on variance (default=5000). If NULL, all genes are
  selected.

- pval_thr:

  cut-off of adjusted p-values to define differentially expressed genes
  from non-spatial linear models. A proportion of genes (`sp_topgenes`)
  under this cut-off will be applied the spatial models. Default=0.05

- pval_adj:

  Method to adjust p-values. Defaults to `FDR`. Other options as
  available from `p.adjust`

- test_type:

  one of `mm`, `t_test`, or `wilcoxon`. Specifies the type of test
  performed.

- sp_topgenes:

  Proportion of differentially expressed genes from non-spatial linear
  models (and controlled by `pval_thr`) to use in differential gene
  expression analysis with spatial linear models. If 0 (zero), no
  spatial models are fit. Default=0.2

- clusters:

  cluster name(s) to test DE genes, as opposed to all clusters.

- pairwise:

  whether or not to carry tests on a pairwise manner. The default is
  `pairwise=F`, meaning that DE genes are tested by comparing each
  cluster to the rest of the pooled cell/spots.

- verbose:

  output progress indicators. If `verbose=0`, no text is shown in
  console. Other values are 1 and 2 indicating increasing level of
  verbosity. Default is `verbose=1`

- cores:

  Number of cores to use in parallelization. If `NULL`, the number of
  cores to use is detected automatically

## Value

a list with one data frame per sample with results of differential gene
expression analysis

## Details

The method tests for differentially expressed genes between groups of
spots/cells (e.g., clusters) in a spatial transcriptomics sample.
Specifically, the function tests for genes with significantly higher or
lower gene expression in one group of spots/cells with respect to the
rest of spots/cells in the sample. The method first runs non-spatial
linear models on the genes to detect differentially expressed genes.
Then spatial linear models with exponential covariance structure are fit
on a subset of genes detected as differentially expressed by the
non-linear models (`sp_topgenes`). If running on clusters detected via
STclust, the user can specify the assignments using the same parameters
(`w`, `k`, `deepSplit`). Otherwise, the assignments are specified by
indicating one of the column names in `x@spatial_meta`. The function
uses [`spaMM::fitme`](https://rdrr.io/pkg/spaMM/man/fitme.html) and is
computationally expensive even on HPC environments. To run the STdiff
using the non-spatial approach (faster), set `sp_topgenes=0`.
