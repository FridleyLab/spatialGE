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

testthat::test_that("SThet errors when no requested gene exists, warns when some do", {
  st <- ge_fixture()

  ## Used to return the STlist untouched, so a typo'd gene name looked like a run
  ## that legitimately found nothing.
  testthat::expect_error(SThet(st, genes = "NoSuchGene", cores = 1, verbose = FALSE),
                         "spatialGE exception.*SThet.*None of the requested genes")

  ## Partially valid: proceeds, naming what was dropped.
  testthat::expect_message(
    res <- SThet(st, genes = c("Gene1", "NoSuchGene"), cores = 1, verbose = FALSE),
    "NoSuchGene.*Not present")
  gm <- get_gene_meta(res, sthet_only = TRUE)
  testthat::expect_true(any(gm$gene == "Gene1"))
})

testthat::test_that("STgradient errors on an annot or ref that matches nothing", {
  st <- ge_fixture()
  sc <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)
  annot <- "stclust_spw0.025_dsplFalse"

  testthat::expect_error(
    STgradient(sc, annot = "NoSuchCol", ref = "1", cores = 1, verbose = FALSE),
    "spatialGE exception.*STgradient.*annot")

  testthat::expect_error(
    STgradient(sc, annot = annot, ref = "NoSuchRef", cores = 1, verbose = FALSE),
    "spatialGE exception.*STgradient.*ref")

  ## The valid combination still runs.
  testthat::expect_no_error(
    STgradient(sc, annot = annot, ref = "1", topgenes = 10, distsumm = "min",
               cores = 1, verbose = FALSE))
})

testthat::test_that("STenrich requires annot and domain together", {
  st <- ge_fixture()
  sc <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)
  sets <- list(gradient = c("GradGene1", "GradGene2", "GradGene3"))
  annot <- "stclust_spw0.025_dsplFalse"

  ## Domain-restricted mode requires both. Supplying one silently skipped the
  ## whole block, so the unrestricted analysis ran and the caller would
  ## reasonably read the result as domain-restricted. Nothing downstream could
  ## detect it, since neither argument is consulted again.
  testthat::expect_error(
    STenrich(sc, gene_sets = sets, annot = annot, min_units = 1, min_genes = 1,
             reps = 10, cores = 1, verbose = FALSE),
    "spatialGE exception.*STenrich.*annot and domain")

  testthat::expect_error(
    STenrich(sc, gene_sets = sets, domain = "1", min_units = 1, min_genes = 1,
             reps = 10, cores = 1, verbose = FALSE),
    "annot and domain")

  ## Both supplied, and neither supplied, are both legitimate.
  testthat::expect_no_error(
    STenrich(sc, gene_sets = sets, annot = annot, domain = "1", min_units = 1,
             min_genes = 1, reps = 10, cores = 1, verbose = FALSE))
  testthat::expect_no_error(
    STenrich(sc, gene_sets = sets, min_units = 1, min_genes = 1, reps = 10,
             cores = 1, verbose = FALSE))
})

testthat::test_that("the on-disk backing never reports a present gene as absent", {
  ## The other half of the rule: silently accepting bad state, or silently
  ## answering nothing, is worse than an error because the caller believes it
  ## worked. All three defects this pins were that half rather than the error
  ## half -- show()/summary() printed blank dimensions, distribution_plots()
  ## messaged every gene as absent and then returned plots that failed only when
  ## printed, and STplot() stopped with "Data was not found in the specified
  ## slot" about data that was right there.
  ##
  ## A manifest answers rownames()/nrow() with NULL rather than erroring, so any
  ## read that leaves the accessors in R/utils_accessors.R fails this way by
  ## default. That is why this is a class guard and not three separate tests.
  st <- ge_fixture()
  bk <- STlist_to_disk(st, file.path(tempdir(), paste0("sge-ge-", as.integer(runif(1, 1, 1e6)))),
                       overwrite = TRUE)

  testthat::expect_no_warning(STplot(bk, genes = "Gene1"))
  testthat::expect_no_message(distribution_plots(bk, genes = "Gene1"))
  testthat::expect_false(any(grepl("( ROIs|spots|cells x  genes)",
                                   capture.output(methods::show(bk)), fixed = TRUE)))
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

testthat::test_that("STdiff rejects a numeric annotation with an actionable message", {
  ## `annot` is documented as "a column name ... containing the groups/clusters
  ## to test", but a perfectly ordinary numeric metadata column (total_counts,
  ## survival_months) used to surface only much later, from dplyr:
  ##   Can't join `x$orig_annot` with `y$orig_annot` due to incompatible types.
  ##   i `x$orig_annot` is a <double>
  ## which names neither STdiff nor the argument to change.
  st <- ge_fixture()
  acol <- "num_annot"
  for (s in names(st@spatial_meta)) {
    st@spatial_meta[[s]][[acol]] <- seq_len(nrow(st@spatial_meta[[s]]))
  }

  testthat::expect_error(
    STdiff(st, annot = acol, topgenes = 5, test_type = "t_test",
           cores = 1, verbose = FALSE),
    "spatialGE exception.*STdiff.*character or factor"
  )
})

testthat::test_that("STdiff rejects an annotation with a single group", {
  st <- ge_fixture()
  acol <- "one_group"
  for (s in names(st@spatial_meta)) {
    st@spatial_meta[[s]][[acol]] <- "only"
  }

  testthat::expect_error(
    STdiff(st, annot = acol, topgenes = 5, test_type = "t_test",
           cores = 1, verbose = FALSE),
    "spatialGE exception.*STdiff.*at least two distinct groups"
  )
})

testthat::test_that("STdiff still accepts a character annotation", {
  ## The guard above must not reject the normal case.
  st <- ge_fixture()
  acol <- "char_annot"
  for (s in names(st@spatial_meta)) {
    st@spatial_meta[[s]][[acol]] <- rep(c("a", "b"), length.out = nrow(st@spatial_meta[[s]]))
  }

  res <- suppressWarnings(
    STdiff(st, annot = acol, topgenes = 5, test_type = "t_test",
           cores = 1, verbose = FALSE)
  )

  testthat::expect_true(is.list(res))
})
