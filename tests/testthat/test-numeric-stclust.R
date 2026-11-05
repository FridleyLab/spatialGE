##
# Numeric assertions for STclust()'s spatially-weighted clustering.
#
# See test-numeric-sthet.R's header for why this tier exists and how it differs
# from the golden gate. Cluster LABELS are order-sensitive (dev/xver/REPORT.md
# found 24 of 24 clusterings had identical partitions but labels could permute),
# so nothing here asserts a specific label on a specific spot. What is asserted
# is the partition: its size profile, and whether it recovers the two-domain
# structure the fixture builds in.
#

stclust_fixture <- function(){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  transform_data(st, method = "log", cores = 1)
}

STCLUST_COL <- "stclust_spw0.025_dsplFalse"

testclust <- function(st){
  STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
          cores = 1, verbose = FALSE)
}

testthat::test_that("STclust produces the pinned partition sizes", {
  sc <- testclust(stclust_fixture())

  ## Sorted counts, not per-label counts: the label a cluster receives is
  ## order-dependent (dev/xver/REPORT.md), its size is not. bravo splits into 3
  ## rather than 2 -- it is the sample carrying the deliberate zero-count spot.
  expected <- list(alpha = c(52L, 68L), bravo = c(20L, 35L, 64L), charlie = c(45L, 75L))
  for(s in names(expected)){
    tab <- sort(as.integer(table(sc@spatial_meta[[s]][[STCLUST_COL]])))
    testthat::expect_equal(tab, expected[[s]], label = paste("cluster sizes", s))
    testthat::expect_equal(sum(tab), nrow(sc@spatial_meta[[s]]),
                           label = paste("sizes total to spot count", s))
  }
})

testthat::test_that("STclust assigns every spot exactly once, with no NA", {
  sc <- testclust(stclust_fixture())

  for(s in names(sc@spatial_meta)){
    labs <- sc@spatial_meta[[s]][[STCLUST_COL]]
    testthat::expect_false(any(is.na(labs)), label = paste("no NA labels", s))
    testthat::expect_equal(length(labs), nrow(sc@spatial_meta[[s]]))
    testthat::expect_gte(length(unique(labs)), 2)
  }
})

testthat::test_that("STclust recovers the fixture's two-domain spatial structure", {
  sc <- testclust(stclust_fixture())

  ## The fixture splits the 12x10 grid at xpos <= 6, with DomainGeneA/B
  ## expressed 14x vs 1x across that boundary. A spatially-weighted clustering
  ## that works must align with it. Measured as the best achievable agreement
  ## over label assignments, since labels themselves are arbitrary.
  for(s in names(sc@spatial_meta)){
    d <- sc@spatial_meta[[s]]
    truth <- ifelse(d$xpos <= 6, "A", "B")
    labs <- as.character(d[[STCLUST_COL]])

    tab <- table(labs, truth)
    ## For each cluster take its majority domain, then total the majorities:
    ## that is the number of spots a best-case label mapping gets right.
    agreement <- sum(apply(tab, 1, max)) / length(labs)

    ## Measured 0.875 (charlie) to 0.925 (alpha). The bound is set below the
    ## worst observed rather than just under the best, so a BLAS or
    ## dependency-version difference cannot trip it -- but it stays far above
    ## the ~0.5 a label-shuffled clustering would score, which is what makes it
    ## evidence that the spatial weighting works at all.
    testthat::expect_gt(agreement, 0.85, label = paste("domain agreement", s))
  }
})

testthat::test_that("STclust writes a column whose name encodes its parameters", {
  sc <- testclust(stclust_fixture())

  ## The column name is API -- downstream calls (STdiff, STenrich, STgradient)
  ## take it as `annot`, so a change to the encoding silently breaks them.
  for(s in names(sc@spatial_meta)){
    testthat::expect_true(STCLUST_COL %in% colnames(sc@spatial_meta[[s]]),
                          label = paste("column present", s))
  }

  ## Two weights produce two independently-named columns rather than
  ## overwriting one.
  sc2 <- STclust(stclust_fixture(), ws = c(0.025, 0.2), ks = "dtc",
                 deepSplit = FALSE, topgenes = 100, cores = 1, verbose = FALSE)
  cols <- grep("^stclust", colnames(sc2@spatial_meta[[1]]), value = TRUE)
  testthat::expect_length(cols, 2)
  testthat::expect_true(any(grepl("0\\.025", cols)))
  testthat::expect_true(any(grepl("0\\.2", cols)))
})
