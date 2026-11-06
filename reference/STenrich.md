# STenrich

Test for spatial enrichment of gene expression sets in ST data sets

## Usage

``` r
STenrich(
  x = NULL,
  samples = NULL,
  gene_sets = NULL,
  score_type = "avg",
  reps = 1000,
  annot = NULL,
  domain = NULL,
  num_sds = 1,
  min_units = 20,
  min_genes = 5,
  pval_adj_method = "BH",
  seed = 12345,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- x:

  an STlist with transformed gene expression

- samples:

  a vector with sample names or indexes to run analysis

- gene_sets:

  a named list of gene sets to test. The names of the list should
  identify the gene sets to be tested

- score_type:

  Controls how gene set expression is calculated. The options are the
  average expression among genes in a set ('avg'), or a GSEA score
  ('gsva'). The default is 'avg'

- reps:

  the number of random samples to be extracted. Default is 1000
  replicates

- annot:

  name of the annotation within `x@spatial_meta` containing the
  spot/cell categories. Needs to be used in conjunction with `domain`

- domain:

  the domain to restrict the analysis. Must exist within the spot/cell
  categories included in the selected annotation (i.e., `annot`)

- num_sds:

  the number of standard deviations to set the minimum gene set
  expression threshold. Default is one (1) standard deviation

- min_units:

  Minimum number of spots with high expression of a pathway for that
  gene set to be considered in the analysis. Defaults to 20 spots or
  cells

- min_genes:

  the minimum number of genes of a gene set present in the data set for
  that gene set to be included. Default is 5 genes

- pval_adj_method:

  the method for multiple comparison adjustment of p-values. Options are
  the same as that of `p.adjust`. Default is 'BH'

- seed:

  the seed number for the selection of random samples. Default is 12345

- cores:

  the number of cores used during parallelization. If NULL (default),
  the number of cores is defined automatically

- verbose:

  either logical or an integer (0, 1, or 2) to increase verbosity

## Value

a list of data frames with the results of the test

## Details

The function performs a randomization test to assess if the sum of
distances between cells/spots with high expression of a gene set is
lower than the sum of distances among randomly selected cells/spots. The
cells/spots are considered as having high gene set expression if the
average expression of genes in a set is higher than the average
expression plus `num_sds` times the standard deviation. Control over the
size of regions with high expression is provided by setting the minimum
number of cells/spots (`min_units`). This method is a modification of
the method devised by Hunter et al. 2021 (zebrafish melanoma study).
