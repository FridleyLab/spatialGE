# STgradient: Tests of gene expression spatial gradients

Calculates Spearman's coefficients to detect genes showing expression
spatial gradients

## Usage

``` r
STgradient(
  x = NULL,
  samples = NULL,
  topgenes = 2000,
  annot = NULL,
  ref = NULL,
  exclude = NULL,
  out_rm = FALSE,
  limit = NULL,
  distsumm = "min",
  min_nb = 3,
  robust = TRUE,
  nb_dist_thr = NULL,
  log_dist = FALSE,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  an STlist with transformed gene expression

- samples:

  the samples on which the test should be executed

- topgenes:

  the number of high-variance genes to be tested. These genes are
  selected in descending order of variance as caclulated using Seurat's
  vst method

- annot:

  the name of a column in `@spatial_meta` containing the tissue domain
  assignments for each spot or cell. These assignments can be generated
  using the `STclust` function

- ref:

  one of the tissue domains in the column specified in `annot`,
  corresponding to the "reference" cluster or domain. Spearman's
  correlations will be calculated using spots assigned to domains other
  than this reference domain (or domains specified in `exclude`).

- exclude:

  optional, a cluster/domain to exclude from the analysis

- out_rm:

  logical (optional), remove gene expression outliers defined by the
  interquartile method. This option is only valid when `robust=F`

- limit:

  limite the analysis to spots/cells with distances to `ref` shorther
  than the value specified here. Useful when gradients might occur at
  smaller scales or when the domain in `ref` is scattered through the
  tissue. Caution must be used due to difficult interpretation of
  imposed limits. It is suggested to run analysis without restricted
  distances in addition for comparison.

- distsumm:

  the distance summary metric to use in correlations. One of `min` or
  `avg`

- min_nb:

  the minimum number of immediate neighbors a spot or cell has to have
  in order to be included in the analysis. This parameter seeks to
  reduce the effect of isolated `ref` spots on the correlation

- robust:

  logical, whether to use robust regression (`MASS` and `sfsmisc`
  packages)

- nb_dist_thr:

  a numeric vector of length two indicating the tolerance interval to
  assign spots/cells to neighborhoods. The wider the range of the
  interval, the more likely distinct neighbors to be considered. If
  NULL, `c(0.75, 1.25)` and `c(0.25, 3)` is assigned for Visium and
  CosMx respectively.

- log_dist:

  logical, whether to apply the natural logarithm to the spot/cell
  distances. It applies to all distances a constant (1e-200) to avoid
  log(0)

- cores:

  the number of cores used during parallelization. If NULL (default),
  the number of cores is defined automatically

- verbose:

  logical, whether to print text to console

## Value

a list of data frames with the results of the test

## Details

The `STgradient` function fits linear models and calculates Spearman
coefficients between the expression of a gene and the minimum or average
distance of spots or cells to a reference tissue domain. In other wordsm
the `STgradient` function can be used to investigate if a gene is
expressed higher in spots/cells closer to a specific reference tissue
domain, compared to spots/cells farther from the reference domain (or
viceversa as indicated by the Spearman's cofficient).
