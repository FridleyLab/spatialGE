##
# Tests for plot_stniche().
#
# plot_stniche() produces three plot types from STniche() output:
#   - "clusters" : quilt map coloured by STniche_DTC
#   - "Zi"       : quilt map coloured by a gene set's Local Moran Z
#   - "ranking"  : bar chart of mean Zi per cluster with SD error bars
#
# Philosophy (same as test-numeric-stniche.R):
#   - cluster LABELS are arbitrary, so we do not assert specific labels;
#   - column names in spatial_meta are API and must be stable;
#   - the function must return a named list of ggplots, one per sample.
#
# The fixture and helper functions are the ones from test-numeric-stniche.R
# (same directory, loaded automatically by testthat).
#

# ---- Fixture helpers (reuse from test-numeric-stniche.R if present) -------
if (!exists("stniche_fixture", mode = "function")) {
    stniche_fixture <- function(){
        set.seed(20260115)
        fx <- fixture_write_list_dfs(fixture_base_data())
        st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
                     cores = 1, verbose = FALSE)
        st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
        transform_data(st, method = "log", cores = 1)
    }
}
if (!exists("stniche_panel", mode = "function")) {
    stniche_panel <- function(){
        list(
            domain_A = c("DomainGeneA", "GeneA_h1", "GeneA_h2",
                         "GeneA_h3", "GeneA_h4", "GeneA_h5"),
            domain_B = c("DomainGeneB", "GeneB_h1", "GeneB_h2",
                         "GeneB_h3", "GeneB_h4", "GeneB_h5")
        )
    }
}
if (!exists(".stniche_panel_available", mode = "function")) {
    .stniche_panel_available <- function(st, panel){
        avail <- Reduce(intersect, lapply(st@tr_counts, rownames))
        lapply(panel, intersect, avail)
    }
}

# ---- Helper: run STniche once and cache (STniche is slow) ------------------
.stniche_cached <- local({
    cache <- NULL
    function(st = NULL) {
        if (!is.null(cache)) return(cache)
        st <- stniche_fixture()
        res <- STniche(
            x         = st,
            gene_sets = .stniche_panel_available(st, stniche_panel()),
            ks        = 2:3,
            knns      = 6,
            out_i     = TRUE,
            cores     = 1,
            verbose   = 0L
        )
        cache <<- list(st = res[[1]], Zi = res[[2]])
        cache
    }
})

# ---- 0. Package wiring -----------------------------------------------------
test_that("plot_stniche is exported from spatialGE", {
    expect_true("plot_stniche" %in% getNamespaceExports("spatialGE"))
})

# ---- 1. Input validation ---------------------------------------------------
test_that("plot_stniche rejects non-STlist input", {
    expect_error(plot_stniche(x = list(), type = "clusters"),
                 "x must be an STlist")
})

test_that("plot_stniche errors when Zi-requiring type has no zi", {
    c  <- .stniche_cached()
    st <- c$st
    expect_error(plot_stniche(st, type = "Zi"),
                 "zi")
    expect_error(plot_stniche(st, type = "ranking"),
                 "zi")
})

test_that("plot_stniche errors when Zi-requiring type has no gene_set", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    expect_error(plot_stniche(st, type = "Zi", zi = Zi),
                 "gene_set")
    expect_error(plot_stniche(st, type = "ranking", zi = Zi),
                 "gene_set")
})

test_that("plot_stniche errors on unknown gene_set", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    expect_error(
        plot_stniche(st, type = "Zi", zi = Zi, gene_set = "NOT_A_REAL_SET"),
        "not found in rownames"
    )
})

test_that("plot_stniche errors when samples matches nothing", {
    c  <- .stniche_cached()
    st <- c$st
    expect_error(
        plot_stniche(st, samples = "nope", type = "clusters"),
        "No matching samples"
    )
})

# ---- 2. type = "clusters" --------------------------------------------------
test_that("plot_stniche(type='clusters') returns a named list of ggplots", {
    c  <- .stniche_cached()
    st <- c$st
    p_list <- plot_stniche(st, type = "clusters", ptsize = 0.8)
    
    expect_type(p_list, "list")
    expect_named(p_list, names(st@spatial_meta), ignore.order = TRUE)
    for (s in names(p_list)) {
        p <- p_list[[s]]
        expect_s3_class(p, "ggplot")
        expect_match(p$labels$title, s, fixed = TRUE)
    }
})

test_that("plot_stniche(type='clusters') handles a single sample", {
    c  <- .stniche_cached()
    st <- c$st
    s1 <- names(st@spatial_meta)[1]
    p_list <- plot_stniche(st, samples = s1, type = "clusters")
    expect_named(p_list, s1)
    expect_s3_class(p_list[[s1]], "ggplot")
})

test_that("plot_stniche(type='clusters') accepts a custom palette", {
    c  <- .stniche_cached()
    st <- c$st
    s1 <- names(st@spatial_meta)[1]
    lv <- sort(unique(as.character(st@spatial_meta[[s1]]$STniche_DTC)))
    my_cols <- grDevices::rainbow(length(lv))
    names(my_cols) <- lv
    
    p_list <- plot_stniche(st, samples = s1, type = "clusters",
                           palette = my_cols)
    expect_s3_class(p_list[[s1]], "ggplot")
    fills <- p_list[[s1]]$scales$scales[[1]]$palette(length(lv))
    expect_setequal(unname(fills), unname(my_cols))
})

# ---- 3. type = "Zi" --------------------------------------------------------
test_that("plot_stniche(type='Zi') returns one ggplot per sample", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    gs <- rownames(Zi)[1]
    p_list <- plot_stniche(st, type = "Zi", zi = Zi, gene_set = gs)
    
    expect_type(p_list, "list")
    expect_named(p_list, names(st@spatial_meta), ignore.order = TRUE)
    for (s in names(p_list)) {
        p <- p_list[[s]]
        expect_s3_class(p, "ggplot")
        expect_match(p$labels$title, gs, fixed = TRUE)
        expect_match(p$labels$title, s,  fixed = TRUE)
    }
})

test_that("plot_stniche(type='Zi') works for every gene set in Zi", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    for (gs in rownames(Zi)) {
        p_list <- plot_stniche(st, type = "Zi", zi = Zi, gene_set = gs)
        expect_type(p_list, "list")
        expect_true(all(vapply(p_list, inherits, logical(1), "ggplot")))
    }
})

test_that("plot_stniche(type='Zi') subsets by sample correctly", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    gs <- rownames(Zi)[1]
    all_samples <- names(st@spatial_meta)
    if (length(all_samples) < 2) skip("need >= 2 samples to test subsetting")
    s1 <- all_samples[1]
    p_list <- plot_stniche(st, samples = s1, type = "Zi",
                           zi = Zi, gene_set = gs)
    expect_named(p_list, s1)
})

# ---- 4. type = "ranking" ---------------------------------------------------
test_that("plot_stniche(type='ranking') returns a bar chart", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    s1 <- names(st@spatial_meta)[1]
    gs <- rownames(Zi)[1]
    p_list <- plot_stniche(st, samples = s1, type = "ranking",
                           zi = Zi, gene_set = gs)
    
    expect_named(p_list, s1)
    p <- p_list[[s1]]
    expect_s3_class(p, "ggplot")
    geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
    expect_true("GeomCol" %in% geoms)
    expect_true("GeomErrorbar" %in% geoms)
})

test_that("plot_stniche(type='ranking') has one bar per cluster", {
    c  <- .stniche_cached()
    st <- c$st; Zi <- c$Zi
    s1 <- names(st@spatial_meta)[1]
    gs <- rownames(Zi)[1]
    
    p_list <- plot_stniche(st, samples = s1, type = "ranking",
                           zi = Zi, gene_set = gs)
    p  <- p_list[[s1]]
    df <- p$data
    n_clusters <- length(unique(as.character(st@spatial_meta[[s1]]$STniche_DTC)))
    
    expect_equal(nrow(df), n_clusters)
    expect_equal(length(unique(as.character(df$STniche_DTC))), n_clusters)
})

# ---- 5. General behaviour --------------------------------------------------
test_that("plot_stniche warns when STniche_DTC is missing", {
    c  <- .stniche_cached()
    st <- c$st
    st_nostic <- st
    for (s in names(st_nostic@spatial_meta)) {
        st_nostic@spatial_meta[[s]]$STniche_DTC <- NULL
    }
    # plot_stniche warns once per sample missing STniche_DTC. capture_warnings()
    # collects all of them so nothing leaks as an uncaught warning; then assert
    # that at least one names the missing column. This replaces the deprecated
    # `all = TRUE` form of expect_warning() (removed in testthat 3.x).
    w <- capture_warnings(plot_stniche(st_nostic, type = "clusters"))
    expect_gte(length(w), 1L)
    expect_true(any(grepl("no STniche_DTC", w, fixed = TRUE)))
})

test_that("plot_stniche respects ptsize for cluster maps", {
    c  <- .stniche_cached()
    st <- c$st
    s1 <- names(st@spatial_meta)[1]
    p_small <- plot_stniche(st, samples = s1, type = "clusters", ptsize = 0.5)
    p_large <- plot_stniche(st, samples = s1, type = "clusters", ptsize = 2.0)
    
    size_small <- p_small[[s1]]$layers[[1]]$aes_params$size
    size_large <- p_large[[s1]]$layers[[1]]$aes_params$size
    expect_equal(size_small, 0.5)
    expect_equal(size_large, 2.0)
})

test_that("plot_stniche returns an empty list when all samples are bad", {
    c  <- .stniche_cached()
    st <- c$st
    st_bad <- st
    for (s in names(st_bad@spatial_meta)) {
        st_bad@spatial_meta[[s]]$STniche_DTC <- NULL
    }
    res <- suppressWarnings(plot_stniche(st_bad, type = "clusters"))
    expect_type(res, "list")
    expect_length(res, 0L)
})

# ---- 6. Composition with STniche output ------------------------------------
test_that("plot_stniche works end-to-end from STniche(out_i = TRUE)", {
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
    st2 <- res[[1]]
    Zi  <- res[[2]]
    
    p_clust <- plot_stniche(st2, type = "clusters")
    p_zi    <- plot_stniche(st2, type = "Zi",
                            zi = Zi, gene_set = rownames(Zi)[1])
    p_rank  <- plot_stniche(st2, type = "ranking",
                            zi = Zi, gene_set = rownames(Zi)[1])
    
    expect_type(p_clust, "list")
    expect_type(p_zi, "list")
    expect_type(p_rank, "list")
    expect_true(all(vapply(p_clust, inherits, logical(1), "ggplot")))
    expect_true(all(vapply(p_zi,    inherits, logical(1), "ggplot")))
    expect_true(all(vapply(p_rank,  inherits, logical(1), "ggplot")))
})