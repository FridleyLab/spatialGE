##
# Seurat interoperability, both directions.
#
# dev/xver/REPORT.md records the Seurat input path as never exercised by any
# harness. It could not have been: STList.R tested `input_samples ==
# 'sample_names'` before the seurat branch, and detect_input_samples(NULL) returns
# a zero-length value, so calling STlist() on a Seurat object with samples=NULL --
# the normal way -- raised "argument is of length zero" before reaching the branch
# that would have handled it.
#

testthat::test_that("read_seurat handles a v5 object (Assay5 + FOV)", {
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(1)
  n <- 60; G <- 30
  cts <- Matrix::Matrix(matrix(rpois(G * n, 4), nrow = G,
                              dimnames = list(paste0("Gene", 1:G), paste0("cell", 1:n))),
                        sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = cts, assay = "Spatial")
  coords <- data.frame(x = rep(1:10, 6), y = rep(1:6, each = 10),
                       cell = paste0("cell", 1:n))
  obj[["slice1"]] <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(coords), type = "centroids", assay = "Spatial")

  ## v5 stores counts in @layers, so the old @assays$Spatial@counts access returned
  ## NULL -- the whole reason this did not work. And v5's FOV has no @coordinates
  ## data.frame, so coordinates come from GetTissueCoordinates().
  testthat::expect_s4_class(obj@assays$Spatial, "Assay5")
  testthat::expect_equal(class(obj@images[["slice1"]])[1], "FOV")

  res <- read_seurat(obj)
  testthat::expect_named(res$counts, "slice1")
  testthat::expect_equal(colnames(res$coords[["slice1"]])[1:3], c("libname", "ypos", "xpos"))

  st <- suppressWarnings(STlist(rnacounts = obj, verbose = FALSE))
  testthat::expect_s4_class(st, "STlist")
  testthat::expect_equal(dim(st@counts[["slice1"]]), c(G, n))
  testthat::expect_s4_class(st@counts[["slice1"]], "dgCMatrix")
})

testthat::test_that("platform is inferred from the image class, not assumed to be visium", {
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(2)
  cts <- Matrix::Matrix(matrix(rpois(20 * 30, 4), nrow = 20,
                              dimnames = list(paste0("Gene", 1:20), paste0("cell", 1:30))),
                        sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = cts, assay = "Spatial")
  obj[["fov1"]] <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(data.frame(x = rep(1:6, 5), y = rep(1:5, each = 6),
                                             cell = paste0("cell", 1:30))),
    type = "centroids", assay = "Spatial")

  ## @misc$platform is not cosmetic: it picks STgradient's default nb_dist_thr
  ## (c(0.75,1.25) for visium vs c(0.25,3) otherwise) and branches in plot_image()
  ## and STplot_interpolation(). Hardcoding 'visium' for every Seurat input, as
  ## before, changed results rather than just a label.
  testthat::expect_equal(seurat_platform(obj), "xenium")

  st <- suppressWarnings(STlist(rnacounts = obj, verbose = FALSE))
  testthat::expect_equal(st@misc[["platform"]], "xenium")
})

testthat::test_that("read_seurat honors the assay argument and fails clearly on a bad one", {
  testthat::skip_if_not_installed("SeuratObject")

  cts <- Matrix::Matrix(matrix(rpois(200, 4), nrow = 10,
                              dimnames = list(paste0("Gene", 1:10), paste0("cell", 1:20))),
                        sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = cts, assay = "MyPanel")
  obj[["s1"]] <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(data.frame(x = rep(1:5, 4), y = rep(1:4, each = 5),
                                             cell = paste0("cell", 1:20))),
    type = "centroids", assay = "MyPanel")

  ## The assay name was hardcoded to 'Spatial', so an object using any other name
  ## -- 'SCT' after SCTransform, 'RNA', a custom panel name -- silently produced
  ## NULL counts. DefaultAssay() is now the default.
  testthat::expect_no_error(read_seurat(obj))
  testthat::expect_error(read_seurat(obj, assay = "NoSuchAssay"),
                         "spatialGE exception.*assay")
})

testthat::test_that("as_Seurat round-trips counts and coordinates", {
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  st <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)

  so <- suppressWarnings(as_Seurat(st))
  testthat::expect_s4_class(so, "Seurat")

  ## One FOV per sample, so each keeps its own coordinate frame instead of all
  ## samples being overlaid in one space.
  testthat::expect_setequal(names(so@images), names(st@counts))

  ## Transformed counts land in the data layer when transform_data() has run.
  testthat::expect_true(all(c("counts", "data") %in% SeuratObject::Layers(so[["Spatial"]])))

  ## STclust assignments survive into @meta.data -- the point of the conversion is
  ## to take spatialGE results into Seurat workflows.
  testthat::expect_true(any(grepl("^stclust", colnames(so@meta.data))))

  ## Cell names are sample-prefixed, since spot barcodes are unique only within a
  ## sample and a merged object needs globally unique names.
  testthat::expect_true(all(grepl("^(alpha|bravo|charlie)_", colnames(so))))

  back <- suppressWarnings(STlist(rnacounts = so, verbose = FALSE))
  testthat::expect_setequal(names(back@counts), names(st@counts))

  a <- st@counts[["alpha"]]
  b <- back@counts[["alpha"]]
  shared <- intersect(rownames(a), rownames(b))
  cn <- paste0("alpha_", colnames(a))
  testthat::expect_identical(unname(as.matrix(a[shared, ])),
                             unname(as.matrix(b[shared, cn[cn %in% colnames(b)]])))
  testthat::expect_equal(sort(st@spatial_meta[["alpha"]]$xpos),
                         sort(back@spatial_meta[["alpha"]]$xpos))
})

testthat::test_that("as_Seurat output is usable by the rest of the package", {
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)

  back <- suppressWarnings(STlist(rnacounts = as_Seurat(st), verbose = FALSE))
  back <- transform_data(back, method = "log", cores = 1)

  ## A round trip that produces an object no downstream function accepts would be
  ## worthless, so assert an analysis actually runs on it.
  sh <- SThet(back, genes = c("Gene1", "GradGene1"), method = "moran",
              cores = 1, verbose = FALSE)
  gm <- get_gene_meta(sh, sthet_only = TRUE)
  testthat::expect_true(any(is.finite(gm$moran_i)))
})

testthat::test_that("every STlist entry point accepts a Seurat object", {
  testthat::skip_if_not_installed("SeuratObject")

  set.seed(3)
  n <- 60; G <- 40
  cts <- Matrix::Matrix(matrix(rpois(G * n, 4), nrow = G,
                              dimnames = list(paste0("Gene", 1:G), paste0("c", 1:n))),
                        sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = cts, assay = "Spatial")
  obj[["s1"]] <- SeuratObject::CreateFOV(
    SeuratObject::CreateCentroids(data.frame(x = rep(1:10, 6), y = rep(1:6, each = 10),
                                             cell = paste0("c", 1:n))),
    type = "centroids", assay = "Spatial")

  ## Coercion happens in assert_stlist(), so it reaches every function that
  ## validates -- which, after this change, is all of them.
  testthat::expect_no_error(suppressWarnings(tissue_names(obj)))
  testthat::expect_no_error(suppressWarnings(spatial_metadata(obj)))
  testthat::expect_no_error(suppressWarnings(summarize_STlist(obj)))

  ## And a real chain works starting from a Seurat object rather than an STlist.
  tr <- suppressWarnings(transform_data(obj, method = "log", cores = 1))
  testthat::expect_s4_class(tr, "STlist")
  sh <- SThet(tr, genes = c("Gene1", "Gene2"), method = "moran", cores = 1, verbose = FALSE)
  testthat::expect_true(any(is.finite(get_gene_meta(sh, sthet_only = TRUE)$moran_i)))
})

testthat::test_that("functions that previously skipped validation now reject bad input", {
  ## Eleven exported functions had no STlist check at all, so a wrong argument
  ## surfaced as whatever internal error came first. They validate now.
  for(fn in c("STgradient", "STenrich", "STdiff", "gene_interpolation", "STplot",
              "plot_counts", "plot_image", "distribution_plots", "spatial_metadata",
              "tissue_names")){
    f <- get(fn, envir = asNamespace("spatialGE"))
    testthat::expect_error(f(42), "STlist", label = paste(fn, "rejects a non-STlist"))
  }
})
