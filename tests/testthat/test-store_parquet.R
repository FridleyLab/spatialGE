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

## Compare what a plot RENDERS, not the object graph -- environments and ggplot
## internals are noisy and non-substantive. Some geoms (ggforce::geom_sina, the
## plot_type='violin' default) draw fresh jitter inside ggplot_build() itself
## rather than when the object is constructed, so two independently built plots
## differ even when the upstream data is identical. Re-seeding before each build
## is how dev/baseline/compare_baseline.R handles the same problem.
built_data <- function(p){
  set.seed(20260115)
  lapply(ggplot2::ggplot_build(p)$data, as.data.frame)
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

testthat::test_that("show() and summary() report real dimensions when backed", {
  st <- store_fixture()
  bk <- backed_fixture(st)

  ## Printing the object is the first thing anyone does with it, and this was the
  ## silent-wrong-answer case rather than an error: nrow()/ncol() return NULL on a
  ## manifest and paste0() renders a zero-length argument as "", so an unrouted
  ## dim read printed "alpha ( ROIs|spots|cells x  genes)" with no warning at all.
  testthat::expect_identical(capture.output(methods::show(bk)),
                             capture.output(methods::show(st)))
  testthat::expect_identical(capture.output(summary(bk)),
                             capture.output(summary(st)))

  ## ...and the numbers are the real ones, which an identity check alone would not
  ## catch -- two blank renderings are also identical to each other.
  expected <- sprintf("\talpha (%d ROIs|spots|cells x %d genes)",
                      ncol(st@counts[["alpha"]]), nrow(st@counts[["alpha"]]))
  testthat::expect_true(any(capture.output(methods::show(bk)) == expected))
})

testthat::test_that("STplot(genes=) works on a backed STlist and matches in memory", {
  st <- store_fixture()
  bk <- backed_fixture(st)
  genes <- c("Gene1", "GradGene1")

  ## plot_spatial_expression() took a list slice of the slot and then tested
  ## rownames() on it, which is NULL on a manifest -- so every sample fell to the
  ## "none of the requested genes are available" branch and the function died
  ## naming the slot rather than the real cause. The gene-SET path
  ## (plot_spatial_geneset) already read through the accessors, so the two STplot
  ## entry points disagreed about the same object.
  a <- STplot(st, genes = genes, data_type = "tr")
  b <- STplot(bk, genes = genes, data_type = "tr")

  ## Panel keys are built in MATRIX gene order, not in the order the caller asked
  ## for (fidelity rule 9), so this also pins that the accessor's two gene-subset
  ## mechanisms -- a logical mask in memory, sorted match() positions on disk --
  ## agree on order and not merely on content.
  testthat::expect_identical(names(b), names(a))
  for(nm in names(a)){
    testthat::expect_identical(built_data(b[[nm]]), built_data(a[[nm]]),
                               label = paste("STplot panel", nm))
  }

  ## The raw slot takes the other branch of the data_type switch.
  ra <- STplot(st, genes = genes, data_type = "raw")
  rb <- STplot(bk, genes = genes, data_type = "raw")
  testthat::expect_identical(names(rb), names(ra))
  testthat::expect_identical(built_data(rb[[1]]), built_data(ra[[1]]))

  ## A gene that really is absent must still drop the sample, warn once per
  ## sample, and then stop. That branch is what made the bug look like an
  ## ordinary "gene not found", so it has to keep behaving identically.
  testthat::expect_error(suppressWarnings(STplot(bk, genes = "NoSuchGene")),
                         "Data was not found in the specified slot")
  w <- testthat::capture_warnings(try(STplot(bk, genes = "NoSuchGene"), silent = TRUE))
  testthat::expect_length(w, length(tissue_names(bk)))
  testthat::expect_true(all(grepl("None of the requested genes are available for sample", w)))
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
