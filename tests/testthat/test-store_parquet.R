##
# Optional on-disk (parquet COO) backing for an STlist's count matrices.
#
# The property that matters is PARITY: a backed STlist must produce results
# identical to an in-memory one. It trades memory for I/O, and nothing else. So
# most of this file runs an analysis both ways and compares, rather than
# inspecting the store's internals.
#

store_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  transform_data(st, method = "log", cores = 1)
}

backed_fixture <- function(st, name = "store"){
  root <- file.path(tempdir(), paste0("sge-", name, "-", as.integer(runif(1, 1, 1e6))))
  STlist_to_disk(st, root, overwrite = TRUE)
}

testthat::test_that("STlist_to_disk writes one file per sample per populated slot", {
  st <- store_fixture()
  bk <- backed_fixture(st)

  testthat::expect_true(is_backed(bk))
  testthat::expect_false(is_backed(st))

  ## 3 samples x {counts, tr_counts}
  root <- bk@misc[["store"]][["root"]]
  testthat::expect_length(list.files(root, recursive = TRUE, pattern = "parquet$"), 6)

  ## sample_order is recorded explicitly so fidelity rule 7 is assertable rather
  ## than inferred from whatever order the filesystem enumerates.
  testthat::expect_identical(bk@misc[["store"]][["sample_order"]], names(st@counts))
  testthat::expect_identical(names(bk@counts), names(st@counts))
})

testthat::test_that("a round trip through disk is bit-identical", {
  st <- store_fixture()
  rt <- STlist_to_memory(backed_fixture(st))

  ## Not "close" -- identical. COO stores the exact doubles, and the manifest
  ## carries dimnames and dimensions, so nothing is inferred on the way back.
  for(s in names(st@counts)){
    testthat::expect_identical(rt@counts[[s]], st@counts[[s]],
                               label = paste("counts round trip", s))
    testthat::expect_identical(rt@tr_counts[[s]], st@tr_counts[[s]],
                               label = paste("tr_counts round trip", s))
  }
  testthat::expect_null(rt@misc[["store"]])
})

testthat::test_that("structural accessors read the manifest without touching data", {
  st <- store_fixture()
  bk <- backed_fixture(st)

  ## This is why the accessor layer was built before the store existed:
  ## base::nrow() returns NULL on a manifest, so an unrouted dim read would have
  ## silently produced NULL rather than failing loudly.
  for(s in names(st@counts)){
    testthat::expect_identical(sample_dim(bk, s, "counts"), sample_dim(st, s, "counts"))
    testthat::expect_identical(sample_genes(bk, s, "counts"), sample_genes(st, s, "counts"))
    testthat::expect_identical(sample_cells(bk, s, "counts"), sample_cells(st, s, "counts"))
    testthat::expect_identical(sample_n_genes(bk, s, "tr_counts"), sample_n_genes(st, s, "tr_counts"))
  }
})

testthat::test_that("a gene-subset read returns genes in manifest order, not requested order", {
  st <- store_fixture()
  bk <- backed_fixture(st)

  ## Fidelity rule 9 arriving from the opposite direction. collect() returns rows
  ## in row-group order, which under a gene sort is ascending gene_idx -- so the
  ## result must be mapped back through the MANIFEST's gene order, never through
  ## the order the caller asked for. Requesting them backwards must still return
  ## them in matrix order.
  genes <- c("DomainGeneA", "GradGene1")
  sub <- sample_tr_counts(bk, "alpha", genes = genes)

  all_genes <- sample_genes(st, "alpha", "tr_counts")
  expected_order <- all_genes[all_genes %in% genes]
  testthat::expect_identical(rownames(sub), expected_order)

  ## Values match the in-memory subset exactly.
  ref <- st@tr_counts[["alpha"]][expected_order, , drop = FALSE]
  testthat::expect_identical(as.matrix(sub), as.matrix(ref))
})

testthat::test_that("all-zero genes and cells survive the round trip", {
  ## COO omits structural zeros, so an all-zero gene has NO rows in the file. Its
  ## existence is known only from the manifest -- which is exactly why the manifest
  ## is the sole authority for dimnames. filter_data()'s defaults are all >= 0, so
  ## such genes are retained and real, and losing them would silently shrink the
  ## matrix (fidelity rule 13c).
  st <- store_fixture()
  testthat::expect_true("ZeroGene" %in% rownames(st@counts[["alpha"]]))
  testthat::expect_equal(sum(st@counts[["alpha"]]["ZeroGene", ]), 0)

  rt <- STlist_to_memory(backed_fixture(st))
  testthat::expect_true("ZeroGene" %in% rownames(rt@counts[["alpha"]]))
  testthat::expect_identical(dim(rt@counts[["alpha"]]), dim(st@counts[["alpha"]]))
  testthat::expect_equal(sum(rt@counts[["alpha"]]["ZeroGene", ]), 0)
})

testthat::test_that("analyses give identical results backed and in memory", {
  st <- store_fixture()
  bk <- backed_fixture(st)
  genes <- c("Gene1", "GradGene1", "DomainGeneA")

  ## The load-bearing property. If these diverge the store is worthless, however
  ## efficient it is.
  a <- SThet(st, genes = genes, method = "moran", cores = 1, verbose = FALSE)
  b <- SThet(bk, genes = genes, method = "moran", cores = 1, verbose = FALSE)
  testthat::expect_identical(get_gene_meta(b, sthet_only = TRUE)$moran_i,
                             get_gene_meta(a, sthet_only = TRUE)$moran_i)

  ca <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)
  cb <- STclust(bk, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)
  col <- "stclust_spw0.025_dsplFalse"
  for(s in names(ca@spatial_meta)){
    testthat::expect_identical(as.character(cb@spatial_meta[[s]][[col]]),
                               as.character(ca@spatial_meta[[s]][[col]]),
                               label = paste("STclust parity", s))
  }

  testthat::expect_identical(summarize_STlist(bk), summarize_STlist(st))
})

testthat::test_that("validate_store detects a missing or altered file", {
  st <- store_fixture()
  bk <- backed_fixture(st)
  testthat::expect_true(validate_store(bk))

  ## The failure mode a user will actually hit: the store directory is moved or
  ## deleted after the STlist was saved with saveRDS(). It must be reported as
  ## that, not as an obscure arrow error deep in an analysis.
  moved <- bk
  moved@counts[["alpha"]]$path <- file.path(tempdir(), "definitely-not-here.parquet")
  testthat::expect_error(validate_store(moved), "spatialGE exception.*missing")

  ## Sample order drifting from what was recorded is a fidelity-rule-7 violation.
  reordered <- bk
  reordered@counts <- reordered@counts[c("bravo", "alpha", "charlie")]
  testthat::expect_error(validate_store(reordered), "spatialGE exception.*order")
})

testthat::test_that("validate_store refuses a non-backed STlist, and deep mode reads data", {
  st <- store_fixture()
  testthat::expect_error(validate_store(st), "not backed")

  bk <- backed_fixture(st)
  testthat::expect_true(validate_store(bk, deep = TRUE))
})

testthat::test_that("STlist_to_disk refuses to clobber an existing store unless told to", {
  st <- store_fixture()
  root <- file.path(tempdir(), paste0("sge-clobber-", as.integer(runif(1, 1, 1e6))))

  bk <- STlist_to_disk(st, root)
  testthat::expect_error(STlist_to_disk(st, root), "already exists")
  testthat::expect_no_error(STlist_to_disk(st, root, overwrite = TRUE))

  testthat::expect_error(STlist_to_disk(st, root = NULL), "single directory path")
})

testthat::test_that("a store under tempdir is marked non-portable", {
  st <- store_fixture()
  bk <- backed_fixture(st)

  ## An STlist is not self-contained once backed: saveRDS() records paths, so a
  ## store under tempdir() cannot survive the session. Flagging it is what lets a
  ## caller know the object is session-scoped before they rely on it.
  testthat::expect_false(bk@misc[["store"]][["portable"]])
})
