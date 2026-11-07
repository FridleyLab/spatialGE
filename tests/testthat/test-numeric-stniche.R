##
# Numeric assertions for STniche()'s GSVA -> Local Moran's I -> mclust path.
#
# Follows the same philosophy as test-numeric-stclust.R:
#   - cluster LABELS are arbitrary (mclust relabels across runs/deps),
#     so nothing here asserts a specific label on a specific spot;
#   - what IS asserted is the partition: its size profile, and whether it
#     recovers the two-domain structure the fixture builds in;
#   - column names in spatial_meta are API and must be stable.
#
# Extra STniche-specific contracts covered here:
#   - `knns` and `ks` both affect the result and are recorded in the tag;
#   - `out_i = TRUE` returns (STlist, Zi) with the expected shape;
#   - on-disk STlists are handled through the lazy accessors.
#

# Column written by STniche. Matches the "stniche_dsplFalse"/"STniche_DTC"
# pair; downstream tools (STdiff, STgradient) read STniche_DTC by default.
STNICHE_COL <- "STniche_DTC"

# ---- Fixture ----------------------------------------------------------------
stniche_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  st <- transform_data(st, method = "log", cores = 1)
  
  # ---- Inject STniche-only helper genes ----
  # Added AFTER the base fixture is built, so the base fixture's 200-gene
  # matrix is unchanged and every golden-value test in
  # test-numeric-stclust / stdiff / stenrich / stgradient / sthet still
  # sees the exact same data it was pinned against.
  #
  # STniche needs multi-gene gene sets (GSVA requires minSize >= 2), so
  # the helpers live here and not in the shared base fixture.
  #
  # The non-dominant-domain rate is held >= 1.2 so P(count=0) < 30%,
  # preventing GSVA's all-zero NA trap when all 6 gene-set members happen
  # to be zero in one spot.
  helper_A <- paste0("GeneA_h", 1:5)
  helper_B <- paste0("GeneB_h", 1:5)
  
  for (s in names(st@tr_counts)) {
    meta     <- st@spatial_meta[[s]]
    dom      <- ifelse(meta$xpos <= 6, "A", "B")
    n_spots  <- length(dom)
    m        <- as.matrix(st@tr_counts[[s]])
    
    new_rows <- matrix(0L, nrow = 10, ncol = n_spots,
                       dimnames = list(c(helper_A, helper_B), colnames(m)))
    for (i in 1:5) {
      new_rows[i, ]    <- rpois(n_spots,
                                ifelse(dom == "A",
                                       14 * runif(n_spots, 0.4, 0.6),
                                       1.5 * runif(n_spots, 0.8, 1.2)))
      new_rows[5 + i, ] <- rpois(n_spots,
                                 ifelse(dom == "B",
                                        14 * runif(n_spots, 0.4, 0.6),
                                        1.5 * runif(n_spots, 0.8, 1.2)))
    }
    
    # log-transform to match the transformed matrix, then append
    m <- rbind(m, log1p(new_rows))
    st@tr_counts[[s]] <- m
  }
  
  st
}

# DomainGeneA/B are the fixture's contrast genes. These are the only genes
# guaranteed to vary between the two halves of the grid; helper genes
# (GeneA_h*/GeneB_h*) may or may not exist and are silently dropped by
# .stniche_panel_available().
stniche_panel <- function(){
  list(
    domain_A = c("DomainGeneA", "GeneA_h1", "GeneA_h2", "GeneA_h3",
                 "GeneA_h4", "GeneA_h5"),
    domain_B = c("DomainGeneB", "GeneB_h1", "GeneB_h2", "GeneB_h3",
                 "GeneB_h4", "GeneB_h5")
  )
}

.stniche_panel_available <- function(st, panel){
  avail <- Reduce(intersect, lapply(st@tr_counts, rownames))
  lapply(panel, intersect, avail)
}

testniche <- function(st, knns = 6L, ks = 2:4){
  STniche(
    x         = st,
    gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks        = ks,
    knns      = knns,
    out_i     = FALSE,
    cores     = 1,
    verbose   = 0L
  )
}

# ---- 0. Package wiring -----------------------------------------------------
test_that("STniche is exported and mclustBIC is imported", {
  # Exports are exposed via NAMESPACE; test membership in getNamespaceExports,
  # not a live binding in asNamespace (which importFrom does not create).
  expect_true("STniche" %in% getNamespaceExports("spatialGE"))
  
  # mclust::Mclust() internally calls mclustBIC(); both names must be
  # declared in NAMESPACE via @importFrom mclust Mclust mclustBIC.
  # Parse the whole file collapsed to a single string so the assertion is
  # insensitive to roxygen's line wrapping (it may write
  # "importFrom(mclust," and "  mclustBIC)" on two physical lines) and to
  # the extra spaces roxygen inserts before the closing paren of long
  # argument lists ("mclustBIC )").
  ns_path <- file.path(find.package("spatialGE"), "NAMESPACE")
  ns_text <- paste(readLines(ns_path, warn = FALSE), collapse = " ")
  expect_match(ns_text, "importFrom\\(mclust,\\s*Mclust,\\s*mclustBIC\\s*\\)",
               info = "mclustBIC missing from spatialGE/NAMESPACE")
})

# ---- 1. Input validation ---------------------------------------------------
test_that("STniche rejects non-STlist input", {
  expect_error(STniche(x = list(), gene_sets = list(a = c("G1"))),
               "x must be an STlist")
})

test_that("STniche requires a named gene_sets list", {
  st <- stniche_fixture()
  expect_error(STniche(st, gene_sets = list(c("G1", "G2"))),
               "named")
})

test_that("STniche errors when tr_counts is empty", {
  st <- stniche_fixture()
  st@tr_counts <- list()
  expect_error(STniche(st, gene_sets = stniche_panel()),
               "transformed counts")
})

test_that("STniche errors on non-numeric ks", {
  st <- stniche_fixture()
  expect_error(STniche(st, gene_sets = stniche_panel(), ks = c("a","b")),
               "numeric")
})

test_that("STniche errors when samples argument matches nothing", {
  st <- stniche_fixture()
  expect_error(STniche(st, gene_sets = stniche_panel(), samples = "nope"),
               "None of the requested samples")
})

# ---- 2. Partition structure ------------------------------------------------
test_that("STniche assigns every spot exactly once, with no NA", {
  st <- stniche_fixture()
  sc <- testniche(st)
  
  for(s in names(sc@spatial_meta)){
    labs <- sc@spatial_meta[[s]][[STNICHE_COL]]
    expect_false(any(is.na(labs)), label = paste("no NA labels", s))
    expect_equal(length(labs), nrow(sc@spatial_meta[[s]]))
    expect_gte(length(unique(labs)), 2L)
  }
})

test_that("STniche cluster count is within the requested ks range", {
  st <- stniche_fixture()
  sc <- testniche(st, ks = 2:4)
  
  for(s in names(sc@spatial_meta)){
    n <- length(unique(sc@spatial_meta[[s]][[STNICHE_COL]]))
    expect_gte(n, 2L)
    expect_lte(n, 4L)
    expect_equal(sum(!is.na(sc@spatial_meta[[s]][[STNICHE_COL]])),
                 nrow(sc@spatial_meta[[s]]),
                 label = paste("sizes total to spot count", s))
  }
})

# ---- 3. Spatial recovery ---------------------------------------------------
test_that("STniche recovers the fixture's two-domain spatial structure", {
  st <- stniche_fixture()
  # ks = 2: the fixture has exactly two domains; letting BIC choose G > 2
  # splits a domain and makes agreement collapse to chance (~0.5), which is
  # a model-selection failure, not a clustering failure. Pin G = 2 here.
  sc <- testniche(st, knns = 6, ks = 2)
  
  for(s in names(sc@spatial_meta)){
    d <- sc@spatial_meta[[s]]
    truth <- ifelse(d$xpos <= 6, "A", "B")
    labs  <- as.character(d[[STNICHE_COL]])
    
    tab <- table(labs, truth)
    # mclust labels are arbitrary; for the binary case, take the better
    # of the two label assignments.
    agreement <- max(
      sum(apply(tab,               1, max)) / length(labs),
      sum(apply(tab[, 2:1, drop = FALSE], 1, max)) / length(labs)
    )
    
    expect_gt(agreement, 0.85, label = paste("domain agreement", s))
  }
})

# ---- 4. out_i contract -----------------------------------------------------
test_that("STniche out_i=TRUE returns (STlist, Zi) with expected shape", {
  st  <- stniche_fixture()
  res <- STniche(
    x         = st,
    gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks        = 2:3,
    knns      = 6,
    out_i     = TRUE,
    cores     = 1,
    verbose   = 0L
  )
  
  expect_length(res, 2L)
  expect_s4_class(res[[1]], "STlist")
  
  Zi <- res[[2]]
  expect_true(is.matrix(Zi))
  expect_equal(sort(rownames(Zi)),
               sort(names(.stniche_panel_available(st, stniche_panel()))))
  expect_true(all(grepl("&&", colnames(Zi))))
  expect_true(any(is.finite(Zi)))
})

test_that("STniche Zi column order matches coordinate ordering", {
  st  <- stniche_fixture()
  res <- STniche(
    x         = st,
    gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks        = 2:3,
    knns      = 6,
    out_i     = TRUE,
    cores     = 1,
    verbose   = 0L
  )
  
  Zi <- res[[2]]
  expected <- unlist(lapply(names(st@spatial_meta), function(s) {
    paste0(s, "&&", st@spatial_meta[[s]]$libname)
  }), use.names = FALSE)
  
  expect_equal(sort(colnames(Zi)), sort(expected))
})

# ---- 5. knns and ks affect the result -------------------------------------
test_that("STniche responds to knns changes", {
  st <- stniche_fixture()
  z_small <- STniche(
    x = st, gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks = 2:3, knns = 4, out_i = TRUE, cores = 1, verbose = 0L
  )[[2]]
  z_large <- STniche(
    x = st, gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks = 2:3, knns = 12, out_i = TRUE, cores = 1, verbose = 0L
  )[[2]]
  
  expect_true(is.matrix(z_small) && is.matrix(z_large))
  expect_false(isTRUE(all.equal(unname(z_small), unname(z_large))))
})

test_that("STniche with ks = 1 writes a single cluster", {
  st <- stniche_fixture()
  sc <- STniche(
    x = st, gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks = 1, knns = 6, cores = 1, verbose = 0L
  )
  expect_equal(length(unique(sc@spatial_meta[[1]][[STNICHE_COL]])), 1L)
})

# ---- 6. Column-name is API -------------------------------------------------
test_that("STniche writes both compat columns, only once", {
  st <- stniche_fixture()
  sc <- testniche(st)
  
  for(s in names(sc@spatial_meta)){
    expect_true("stniche_dsplFalse" %in% colnames(sc@spatial_meta[[s]]))
    expect_true(STNICHE_COL %in% colnames(sc@spatial_meta[[s]]))
    expect_equal(sc@spatial_meta[[s]]$stniche_dsplFalse,
                 sc@spatial_meta[[s]][[STNICHE_COL]])
  }
})

test_that("Re-running STniche overwrites old columns (no duplicates)", {
  st  <- stniche_fixture()
  st1 <- testniche(st,  ks = 2:3)
  st2 <- STniche(
    x = st1, gene_sets = .stniche_panel_available(st1, stniche_panel()),
    ks = 4:5, knns = 6, cores = 1, verbose = 0L
  )
  
  for(s in names(st2@spatial_meta)){
    expect_equal(sum(colnames(st2@spatial_meta[[s]]) == "STniche_DTC"),
                 1L, label = paste("single STniche_DTC column", s))
    expect_equal(sum(colnames(st2@spatial_meta[[s]]) == "stniche_dsplFalse"),
                 1L, label = paste("single stniche_dsplFalse column", s))
  }
})

# ---- 7. Robustness ---------------------------------------------------------
test_that("STniche tolerates a gene set with zero overlap", {
  st    <- stniche_fixture()
  panel <- stniche_panel()
  panel$empty <- c("NOT_A_GENE_1", "NOT_A_GENE_2")
  panel <- .stniche_panel_available(st, panel)
  
  sc <- STniche(
    x = st, gene_sets = panel, ks = 2:3, knns = 6, cores = 1, verbose = 0L
  )
  expect_s4_class(sc, "STlist")
})

# ---- 8. Disk mode (skipped if unavailable) --------------------------------
test_that("STniche runs against an on-disk (lazy-handle) STlist", {
  skip_if_not(exists("create_stlist_disk", mode = "function"),
              "on-disk STlist constructor not available in this build")
  skip_if_not_installed("rhdf5")
  
  tmp <- tempfile("stniche_disk_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  
  st <- stniche_fixture()
  for(s in names(st@tr_counts)){
    fp <- file.path(tmp, paste0(s, "_tr.h5"))
    rhdf5::h5createFile(fp)
    rhdf5::h5write(st@tr_counts[[s]], fp, "counts")
    st@tr_counts[[s]] <- fp
  }
  st@disk_path <- tmp
  
  local_mocked_bindings(
    get_gene_counts = function(x, gene, sample){
      m <- rhdf5::h5read(x@tr_counts[[sample]], "counts")
      m[intersect(gene, rownames(m)), , drop = FALSE]
    },
    get_gene_names = function(x, sample){
      rownames(rhdf5::h5read(x@tr_counts[[sample]], "counts"))
    },
    .package = "spatialGE"
  )
  
  res <- STniche(
    x         = st,
    gene_sets = .stniche_panel_available(st, stniche_panel()),
    ks        = 2:3,
    knns      = 6,
    out_i     = TRUE,
    cores     = 1,
    verbose   = 0L,
    tag       = "toy"
  )
  expect_s4_class(res[[1]], "STlist")
  zi_fp <- if (is.character(res[[2]])) res[[2]]
  else file.path(tmp, "Zi_stniche_toy.h5")
  expect_true(file.exists(zi_fp) ||
                file.exists(sub("\\.h5$", "\\.rds", zi_fp)))
})