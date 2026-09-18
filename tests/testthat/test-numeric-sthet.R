##
# Numeric assertions for SThet()'s spatial-autocorrelation statistics.
#
# This is the platform-portable tier: unlike the golden gate (which is
# bit-exact and therefore same-machine only -- see dev/baseline/compare_baseline.R),
# these assert literal values at an explicit tolerance and run by default,
# including on CI. Before this tier existed, none of the six compute-heavy
# algorithms had a single direct numeric test; they were covered only
# transitively by a gate that skips on every fresh clone.
#
# TOLERANCE. 1e-6 throughout. Moran's I and Geary's C here are sums and
# products over a distance matrix, so cross-platform agreement is nearer 1e-12;
# 1e-6 leaves room for BLAS and dependency-version differences while staying
# far tighter than any real algorithmic change, which moves these in the second
# or third decimal (dev/xver/REPORT.md measured the kNN variant's order
# sensitivity at 1.8e-2).
#
# ORDER SENSITIVITY. Deliberately uses SThet() WITHOUT k. Per fidelity rule 13,
# SThet(k=) breaks tied k-th distances by spot position and moves by ~2e-2
# under reordering, which would make a pinned value meaningless. The default
# distance-based weighting agrees to ~1.5e-15 under reordering.
#
# Note the tension rule 13(b) creates: the order-insensitive default is also
# create_listw_from_dist()'s complete graph, which is the package's worst
# scaling path. These tests pin its *values* so a replacement can be shown
# equivalent.
#

sthet_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  ## bravo carries a deliberate zero-count spot, which filter_data warns about.
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  transform_data(st, method = "log", cores = 1)
}

SAMPLES <- c("alpha", "bravo", "charlie")

testthat::test_that("SThet Moran's I ranks the fixture's genes by designed spatial structure", {
  st <- sthet_fixture()
  genes <- c("Gene1", "GradGene1", "DomainGeneA")
  sh <- SThet(st, genes = genes, method = "moran", cores = 1, verbose = FALSE)

  gm <- get_gene_meta(sh, sthet_only = TRUE)
  mi <- stats::setNames(gm$moran_i, paste(gm$sample, gm$gene, sep = "|"))

  ## The fixture builds three regimes on purpose: Gene1 is Poisson(4) with no
  ## spatial term, GradGene1 has a linear gradient in x, DomainGeneA is a
  ## two-domain step. Moran's I must order them accordingly in EVERY sample.
  ## This is the assertion that survives any implementation change; the pinned
  ## values below are the regression detector.
  for(s in SAMPLES){
    testthat::expect_gt(mi[[paste0(s, "|DomainGeneA")]], mi[[paste0(s, "|GradGene1")]])
    testthat::expect_gt(mi[[paste0(s, "|GradGene1")]], mi[[paste0(s, "|Gene1")]])
    ## A gene with no spatial structure sits near zero, not merely below the
    ## others.
    testthat::expect_lt(abs(mi[[paste0(s, "|Gene1")]]), 0.05)
  }
})

testthat::test_that("SThet Moran's I matches pinned values", {
  st <- sthet_fixture()
  genes <- c("Gene1", "GradGene1", "DomainGeneA")
  sh <- SThet(st, genes = genes, method = "moran", cores = 1, verbose = FALSE)

  gm <- get_gene_meta(sh, sthet_only = TRUE)
  mi <- stats::setNames(gm$moran_i, paste(gm$sample, gm$gene, sep = "|"))

  expected <- c(
    "alpha|Gene1"         = -0.019182504865,
    "alpha|GradGene1"     =  0.114399515879,
    "alpha|DomainGeneA"   =  0.209595137645,
    "bravo|Gene1"         = -0.005225036527,
    "bravo|GradGene1"     =  0.103622707754,
    "bravo|DomainGeneA"   =  0.249129040457,
    "charlie|Gene1"       = -0.022435305149,
    "charlie|GradGene1"   =  0.086576657061,
    "charlie|DomainGeneA" =  0.229314160187
  )

  for(k in names(expected)){
    testthat::expect_equal(mi[[k]], expected[[k]], tolerance = 1e-6,
                           label = paste("moran_i", k))
  }
})

testthat::test_that("SThet Geary's C matches pinned values and inverts Moran's I", {
  st <- sthet_fixture()
  genes <- c("Gene1", "DomainGeneA")
  sg <- SThet(st, genes = genes, method = "geary", cores = 1, verbose = FALSE)

  gm <- get_gene_meta(sg, sthet_only = TRUE)
  gc_ <- stats::setNames(gm$geary_c, paste(gm$sample, gm$gene, sep = "|"))

  expected <- c(
    "alpha|Gene1"         = 1.0193410705,
    "alpha|DomainGeneA"   = 0.7917956706,
    "bravo|Gene1"         = 1.0114828890,
    "bravo|DomainGeneA"   = 0.7208013213,
    "charlie|Gene1"       = 0.9886193414,
    "charlie|DomainGeneA" = 0.7804281369
  )
  for(k in names(expected)){
    testthat::expect_equal(gc_[[k]], expected[[k]], tolerance = 1e-6,
                           label = paste("geary_c", k))
  }

  ## Geary's C runs opposite to Moran's I: ~1 means no autocorrelation, below 1
  ## means positive. The spatially structured gene must land well below 1 and
  ## the unstructured one near it. Catches a sign or normalization inversion
  ## that pinned values alone would not explain.
  for(s in SAMPLES){
    testthat::expect_lt(gc_[[paste0(s, "|DomainGeneA")]], 0.85)
    testthat::expect_equal(gc_[[paste0(s, "|Gene1")]], 1, tolerance = 0.05)
  }
})

testthat::test_that("SThet writes only the requested statistic and leaves the other NA", {
  st <- sthet_fixture()
  sg <- SThet(st, genes = c("Gene1"), method = "geary", cores = 1, verbose = FALSE)

  gm <- sg@gene_meta[["alpha"]]
  testthat::expect_true(all(c("moran_i", "geary_c") %in% colnames(gm)))

  ## Both columns are created, but a geary run must not populate moran_i --
  ## column presence and type are API (fidelity rule 8).
  row <- gm[gm$gene == "Gene1", ]
  testthat::expect_true(is.na(row$moran_i))
  testthat::expect_false(is.na(row$geary_c))
})

testthat::test_that("SThet tolerates non-contiguous spatial_meta row names", {
  st <- sthet_fixture()

  ## @spatial_meta is a plain data.frame, not a tibble, so row-subsetting it
  ## with base `[` leaves row names like 1, 5, 9, ... behind, and
  ## tibble::column_to_rownames() then aborts with "`.data` must be a data frame
  ## without row names". Both listw builders in R/utils.R hit that, so any
  ## STlist whose spatial_meta had been row-subset by base `[` -- an ordinary
  ## thing for a caller to do -- made SThet() fail with an error naming neither
  ## SThet nor the real cause. This is the same footgun that made STenrich()'s
  ## domain= path unreachable through 1.2.2.
  st@spatial_meta[["alpha"]] <- st@spatial_meta[["alpha"]][
    seq(1, nrow(st@spatial_meta[["alpha"]]), by = 2), , drop = FALSE]
  keep <- st@spatial_meta[["alpha"]][["libname"]]
  st@counts[["alpha"]] <- st@counts[["alpha"]][, keep, drop = FALSE]
  st@tr_counts[["alpha"]] <- st@tr_counts[["alpha"]][, keep, drop = FALSE]

  testthat::expect_false(identical(rownames(st@spatial_meta[["alpha"]]),
                                   as.character(seq_along(keep))))

  testthat::expect_no_error(
    sh <- SThet(st, genes = c("GradGene1"), method = "moran", cores = 1, verbose = FALSE)
  )
  gm <- get_gene_meta(sh, sthet_only = TRUE)
  testthat::expect_true(any(!is.na(gm$moran_i)))

  ## The kNN builder shares the defect and the fix.
  testthat::expect_no_error(
    SThet(st, genes = c("GradGene1"), method = "moran", k = 4, cores = 1, verbose = FALSE)
  )
})

testthat::test_that("SThet is invariant to spot order for the default weighting", {
  st <- sthet_fixture()
  genes <- c("GradGene1", "DomainGeneA")

  base <- SThet(st, genes = genes, method = "moran", cores = 1, verbose = FALSE)

  ## Permute spots within one sample and re-run. Fidelity rule 13 establishes
  ## that the default distance-based listw is order-insensitive to ~1.5e-15
  ## while SThet(k=) is not, so this is the property that makes the pinned
  ## values above legitimate rather than incidental.
  st2 <- st
  perm <- rev(seq_len(ncol(st2@tr_counts[["alpha"]])))
  st2@counts[["alpha"]] <- st2@counts[["alpha"]][, perm, drop = FALSE]
  st2@tr_counts[["alpha"]] <- st2@tr_counts[["alpha"]][, perm, drop = FALSE]
  sm <- st2@spatial_meta[["alpha"]][perm, , drop = FALSE]
  ## Row-subsetting a plain data.frame with base `[` leaves non-contiguous row
  ## names behind, which tibble::column_to_rownames() then rejects -- the
  ## footgun documented in dev/ARCHITECTURE.md's slot-contract section. Reset
  ## them so this test exercises order-invariance rather than that separate
  ## defect (which create_listw_from_dist() is also subject to; see the
  ## companion commit adding remove_rownames() there).
  rownames(sm) <- NULL
  st2@spatial_meta[["alpha"]] <- sm

  permuted <- SThet(st2, genes = genes, method = "moran", cores = 1, verbose = FALSE)

  gm_b <- get_gene_meta(base, sthet_only = TRUE)
  gm_p <- get_gene_meta(permuted, sthet_only = TRUE)
  key_b <- paste(gm_b$sample, gm_b$gene, sep = "|")
  key_p <- paste(gm_p$sample, gm_p$gene, sep = "|")

  for(k in key_b[grepl("^alpha", key_b)]){
    testthat::expect_equal(gm_p$moran_i[key_p == k], gm_b$moran_i[key_b == k],
                           tolerance = 1e-10, label = paste("order-invariant", k))
  }
})

testthat::test_that("computeMoranGeary reproduces spdep over the complete graph", {
  ## The C++ path replaced create_listw_from_dist() + per-gene spdep calls. This
  ## builds the old listw explicitly and asserts agreement, so the equivalence is
  ## pinned in the repository rather than resting on a one-off check.
  ##
  ## Agreement is ~1e-15 relative, not bit-identical: S0 IS bit-identical (the C++
  ## accumulates it in one continuous long double in spdep's own visit order), and
  ## the residue is summation order inside lag.listw's dot product. 1e-15 is four
  ## orders of magnitude below the order sensitivities fidelity rule 13 documents
  ## for this statistic family.
  set.seed(42)
  n <- 200
  side <- ceiling(sqrt(n))
  g <- expand.grid(xpos = seq_len(side), ypos = seq_len(side))[seq_len(n), ]
  coords <- as.matrix(g[, c("ypos", "xpos")])
  rownames(coords) <- paste0("s", seq_len(n))

  G <- 4
  expr <- matrix(as.double(rpois(G * n, 4)), nrow = G,
                 dimnames = list(paste0("g", 1:G), rownames(coords)))

  ## Rebuild the old complete-graph listw exactly as create_listw_from_dist did.
  adj <- as.matrix(stats::dist(coords))
  adj <- adj / max(adj)
  diag(adj) <- NA
  dimnames(adj) <- list(rownames(coords), rownames(coords))
  nb <- list(); wt <- list()
  for(id in seq_len(nrow(adj))){
    idx <- seq_len(nrow(adj))
    nb[[id]] <- idx[-id]; names(nb[[id]]) <- idx[-id]
    wt[[id]] <- 1 / (as.vector(adj[id, ])[-id])
  }
  class(nb) <- "nb"
  lw <- list(style = "B", neighbours = nb, weights = wt)
  class(lw) <- c("listw", "nb")

  ref_i <- vapply(seq_len(G), function(k)
    spdep::moran(expr[k, ], lw, n = n, S0 = spdep::Szero(lw))$I, numeric(1))
  ref_c <- vapply(seq_len(G), function(k)
    spdep::geary(expr[k, ], lw, n = n, n1 = n - 1, S0 = spdep::Szero(lw))$C, numeric(1))

  zt <- t(sweep(expr, 1, rowMeans(expr), "-"))
  res <- computeMoranGeary(coords, zt, TRUE, TRUE)

  testthat::expect_equal(res$S0, spdep::Szero(lw), tolerance = 0)
  testthat::expect_equal(res$moran, ref_i, tolerance = 1e-12)
  testthat::expect_equal(res$geary, ref_c, tolerance = 1e-12)
})

testthat::test_that("computeMoranGeary behaves on random and degenerate input", {
  set.seed(7)
  n <- 400
  side <- ceiling(sqrt(n))
  g <- expand.grid(xpos = seq_len(side), ypos = seq_len(side))[seq_len(n), ]
  coords <- as.matrix(g[, c("ypos", "xpos")])

  ## Spatially unstructured data: Moran's I near 0 and Geary's C near 1 are the
  ## theoretical expectations, and are what caught nothing being transposed.
  zt <- matrix(rnorm(n * 3), nrow = n)
  zt <- sweep(zt, 2, colMeans(zt), "-")
  res <- computeMoranGeary(coords, zt, TRUE, TRUE)
  testthat::expect_true(all(abs(res$moran) < 0.05))
  testthat::expect_true(all(abs(res$geary - 1) < 0.05))

  ## A constant gene has z == 0 throughout, so zz == 0 and both statistics are
  ## 0/0. spdep returns NaN there; so must this, or SThet would silently write a
  ## fabricated number into @gene_meta.
  zt0 <- cbind(zt[, 1], rep(0, n))
  res0 <- computeMoranGeary(coords, zt0, TRUE, TRUE)
  testthat::expect_true(is.nan(res0$moran[2]))
  testthat::expect_true(is.nan(res0$geary[2]))
  testthat::expect_false(is.nan(res0$moran[1]))

  ## Guard the input contract rather than letting a transposed matrix through.
  testthat::expect_error(computeMoranGeary(coords, t(zt), TRUE, TRUE), "one row per spot")
})
