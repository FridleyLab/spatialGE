##
# Every user-reachable failure must be an informative spatialGE error, never an
# internal one leaking out of Matrix, dplyr, tibble, spdep, gstat or base.
#
# This file pins that requirement for the cases found by auditing the main entry
# points against plausible user mistakes. See the "Failing gracefully" section of
# .claude/CLAUDE.md for the rule and how to add to it.
#
# The bar is: the message names the function and the offending argument, and
# tells the caller what to do. "object 'tr_results' not found" and "attempt to
# select less than one element in get1index" both fail that bar -- they name
# internal variables the caller has no way to act on.
#

ge_fixture <- function(transformed = TRUE){
  set.seed(20260115)
  fx <- fixture_write_list_dfs(fixture_base_data())
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords,
               cores = 1, verbose = FALSE)
  st <- suppressWarnings(filter_data(st, gene_minreads = 0, spot_minreads = 0))
  if(!transformed) return(st)
  transform_data(st, method = "log", cores = 1)
}

testthat::test_that("transform_data rejects an unrecognized method", {
  st <- ge_fixture(transformed = FALSE)

  ## Previously fell through the if/else without assigning tr_results and failed
  ## ten lines later with "object 'tr_results' not found".
  testthat::expect_error(transform_data(st, method = "bogus", cores = 1),
                         "spatialGE exception.*transform_data.*Unrecognized method")
  testthat::expect_error(transform_data(st, method = c("log", "sct"), cores = 1),
                         "Unrecognized method")

  ## Both documented methods still work.
  testthat::expect_no_error(transform_data(st, method = "log", cores = 1))
})

testthat::test_that("SThet requires transformed counts", {
  st <- ge_fixture(transformed = FALSE)

  ## Previously reached resolve_samples() with an empty slot and failed inside
  ## tibble with "Can't extract columns past the end", naming neither SThet nor
  ## transform_data.
  testthat::expect_error(SThet(st, genes = "Gene1", cores = 1, verbose = FALSE),
                         "spatialGE exception.*SThet.*transform_data")
})

testthat::test_that("SThet rejects a k that spdep cannot satisfy", {
  st <- ge_fixture()
  n <- ncol(st@tr_counts[[1]])

  ## spdep::knearneigh() needs k < the number of points. Previously surfaced as
  ## "knearneigh: fewer data points than k" from inside a worker.
  testthat::expect_error(SThet(st, genes = "Gene1", k = n + 50, cores = 1, verbose = FALSE),
                         "spatialGE exception.*SThet.*k must be smaller")
  testthat::expect_error(SThet(st, genes = "Gene1", k = n, cores = 1, verbose = FALSE),
                         "k must be smaller")
  testthat::expect_error(SThet(st, genes = "Gene1", k = 0, cores = 1, verbose = FALSE),
                         "k must be smaller")

  ## A workable k still runs.
  testthat::expect_no_error(SThet(st, genes = "Gene1", k = 4, cores = 1, verbose = FALSE))
})

testthat::test_that("STenrich rejects an empty or unnamed gene_sets", {
  st <- ge_fixture()
  sc <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)

  ## An empty list gave a zero-row combo and then "attempt to select less than
  ## one element in get1index".
  testthat::expect_error(
    STenrich(sc, gene_sets = list(), min_units = 1, min_genes = 1, reps = 10,
             cores = 1, verbose = FALSE),
    "spatialGE exception.*STenrich.*gene_sets")

  ## An unnamed list is equally unusable -- it produces gene_set=NA rows that
  ## silently test nothing, which is worse than an error.
  testthat::expect_error(
    STenrich(sc, gene_sets = list(c("Gene1", "Gene2")), min_units = 1, min_genes = 1,
             reps = 10, cores = 1, verbose = FALSE),
    "gene_sets")

  testthat::expect_error(
    STenrich(sc, gene_sets = "not a list", min_units = 1, min_genes = 1, reps = 10,
             cores = 1, verbose = FALSE),
    "gene_sets")
})

testthat::test_that("gene_interpolation reports absent genes and fails clearly when none remain", {
  st <- ge_fixture()

  ## Previously built empty kriging entries and then failed deep inside the gstat
  ## formula with "object 'xpos' not found".
  testthat::expect_error(
    gene_interpolation(st, genes = "NoSuchGene", cores = 1, verbose = FALSE),
    "spatialGE exception.*gene_interpolation.*None of the requested genes")

  ## A partially-valid request proceeds on the genes that exist, and says which
  ## it dropped rather than dropping them silently.
  testthat::expect_message(
    res <- gene_interpolation(st, genes = c("GradGene1", "NoSuchGene"),
                               cores = 1, verbose = FALSE),
    "NoSuchGene.*Not present")
  testthat::expect_true("GradGene1" %in% names(res@gene_krige))
  testthat::expect_false("NoSuchGene" %in% names(res@gene_krige))
})

testthat::test_that("no user-facing failure leaks an internal variable name", {
  st <- ge_fixture()
  stu <- ge_fixture(transformed = FALSE)

  ## A guard against regressions of the whole class: these are the specific
  ## internal names that used to escape. None may reappear in any message.
  leaked <- c("tr_results", "get1index", "xpos", "knearneigh",
              "Can't extract columns past the end")

  attempts <- list(
    function() transform_data(stu, method = "bogus", cores = 1),
    function() SThet(stu, genes = "Gene1", cores = 1, verbose = FALSE),
    function() SThet(st, genes = "Gene1", k = 500, cores = 1, verbose = FALSE),
    function() gene_interpolation(st, genes = "NoSuchGene", cores = 1, verbose = FALSE)
  )

  for(f in attempts){
    msg <- tryCatch({ f(); NA_character_ }, error = function(e) conditionMessage(e))
    testthat::expect_false(is.na(msg))
    for(bad in leaked){
      testthat::expect_false(grepl(bad, msg, fixed = TRUE),
                             label = sprintf("message leaks '%s': %s", bad, msg))
    }
    testthat::expect_match(msg, "spatialGE exception")
  }
})
