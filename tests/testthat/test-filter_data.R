##
# filter_stage_genes_expr(), STEP 4 of filter_data()'s fixed filtering order.
# Regression test for the bug where the @gene_meta branch indexed with a
# malformed 3-subscript expression and rownames() (tibbles have no real
# rownames, so this returns row-number strings) instead of the `gene`
# column, making that branch a silent no-op. filter_data()'s own STEP 5
# (filter_stage_counts) happens to re-sync @gene_meta from @tr_counts
# afterward regardless, which masks the bug end-to-end -- so this test
# calls filter_stage_genes_expr() directly rather than through filter_data().
#

testthat::test_that("filter_stage_genes_expr removes matching genes from gene_meta, not just counts", {
  f <- fixtures$list_dfs
  st <- STlist(rnacounts = f$rnacounts, spotcoords = f$spotcoords, cores = 1)
  st <- transform_data(st, method = "log", cores = 1)

  target_gene <- "GradGene1"
  st_filtered <- filter_stage_genes_expr(st, names(st@counts), paste0("^", target_gene, "$"))

  for (s in names(st_filtered@gene_meta)) {
    testthat::expect_false(target_gene %in% st_filtered@gene_meta[[s]][["gene"]])
    testthat::expect_false(target_gene %in% rownames(st_filtered@counts[[s]]))
    testthat::expect_equal(nrow(st_filtered@gene_meta[[s]]), nrow(st@gene_meta[[s]]) - 1)
  }
})
