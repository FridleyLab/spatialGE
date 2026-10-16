## The expectation table: one record per difference we claim to already
## understand, each naming the commit that caused it.
##
## The table is used in BOTH directions, and the second one is the point:
##
##   forward  -- an observed difference that matches a record is EXPECTED_FIX
##   inverse  -- a record that matches NO observed difference is MISSING_EXPECTED
##
## The inverse check is what would have caught "the pseudobulk_heatmap CV fix is
## entirely unverified" a month ago. An expected fix that produces no difference
## means either the scenario does not reach it or the fix does not work, and both
## are failures rather than silence.
##
## Anything matching no record at all is a REGRESSION. Fail closed.

XV_EXPECT <- list(

  ## ---- delta 2: STenrich(annot=, domain=) always errored on main -----------
  list(id = "fix02_stenrich_domain", commit = "3cd478f", delta = 2,
       scenario = NA, path = "^S0b?/STenrich/",
       class = "EXPECTED_FIX",
       rationale = paste("Domain-restricted STenrich subsets spatial_meta rows, leaving",
                         "non-contiguous row names that tibble::column_to_rownames()",
                         "rejects. Every domain= call in 1.2.1/1.2.2 errored.")),

  ## ---- delta 3: Xenium MEX import ----------------------------------------
  list(id = "fix03_xenium_mex", commit = "a60651a", delta = 3,
       scenario = "S1", path = "^S1/STlist_xenium_mex/",
       class = "EXPECTED_FIX",
       rationale = paste("import_xenium.R defined a second function named import_visium,",
                         "shadowing the real one, and read_xenium_outs called a",
                         "non-existent import_xenium_mex.")),

  ## ---- delta 6: xlsx samples-file detection ------------------------------
  list(id = "fix06_xlsx_detect", commit = "777bb59", delta = 6,
       scenario = "S1", path = "^S1/STlist_list_dfs_samplesfile_xlsx/",
       class = "EXPECTED_FIX",
       rationale = "grepl('//.xlsx$') never matched, so an xlsx samples file was not detected."),

  ## ---- delta 8: undefined h5_file in the import warning path -------------
  list(id = "fix08_h5_file_warning", commit = "ab76c49", delta = 8,
       scenario = "S15", path = "^S15/s15_h5_multigroup/",
       class = "EXPECTED_FIX",
       rationale = paste("The warning path referenced an undefined h5_file, so raising the",
                         "warning itself errored. Only visible because warnings are",
                         "captured as payload here.")),

  ## ---- delta 5: rm_genes_expr did not filter @gene_meta ------------------
  list(id = "fix05_rm_genes_expr", commit = "d8a2e37", delta = 5,
       scenario = "S2", path = "^S2/s2_rm_genes_expr/.*gene_meta",
       class = "EXPECTED_FIX",
       rationale = paste("filter_data indexed @gene_meta with a malformed 3-subscript index",
                         "and matched rownames() of a tibble ('1'..'200'), so the stage was a",
                         "silent no-op. Deliberately absent from NEWS because a later stage",
                         "re-syncs @gene_meta -- which is the claim being tested here.")),

  ## ---- delta 9: inverted coefficient of variation ------------------------
  list(id = "fix09_pbulk_cv", commit = "e7917bc", delta = 9,
       scenario = "S4", path = "^S4/s4_heatmap",
       class = "EXPECTED_FIX",
       rationale = paste("CV was computed as mean/sd instead of sd/mean, so sorting",
                         "decreasing selected the LEAST variable genes. Requires >=4",
                         "samples, which the committed fixture does not have.")),
  list(id = "fix09_pbulk_cv_real", commit = "e7917bc", delta = 9,
       scenario = "S12", path = "^S12/s12_heatmap",
       class = "EXPECTED_FIX", optional = TRUE,
       rationale = "Same fix on real Visium data, which has >=4 samples natively."),

  ## ---- delta 10: STdiff_volcano's wrong numeric check -------------------
  list(id = "fix10_volcano_samples", commit = "b8de052", delta = 10,
       scenario = "S5", path = "^S5/s5_volcano_numeric/",
       class = "EXPECTED_FIX",
       rationale = "is.numeric(x) should have been is.numeric(samples), so numeric samples= was ignored."),

  ## ---- delta 11: STlist() without rnacounts -----------------------------
  list(id = "fix11_null_rnacounts", commit = "6954044", delta = 11,
       scenario = "S6", path = "^S6/s6_(null_rnacounts|null_both|no_args)/(error_message|ok)",
       class = "EXPECTED_FIX",
       rationale = "cores=NULL with rnacounts=NULL raised 'argument is of length zero'."),

  ## ---- delta 13: STgradient's pval_warn state leak ----------------------
  ## RESOLVED as unreachable-from-here, and verified elsewhere instead. This
  ## reported NOT_OBSERVED for three consecutive runs, which is the inverse check
  ## doing its job. The resolution is not to delete the record but to record why
  ## the integration route cannot reach it, and where the fix IS pinned.
  ##
  ## pval_warn is set only when cor.test(method='spearman') warns specifically
  ## about a zero standard deviation, i.e. when an input vector is constant. The
  ## only gene still constant after library-size normalisation is ZeroGene --
  ## ConstGene is not, because dividing a constant count by a varying library size
  ## gives sd 0.0389. And ZeroGene never reaches the correlation: the robust
  ## regression fails on it first, so its row carries 'rob_regr_no_convergence'.
  ## No fixture-based scenario can therefore produce 'zero_st_deviation' through
  ## STgradient(), whatever topgenes or robust setting is used.
  ##
  ## Verified directly instead, in tests/testthat/test-numeric-stgradient.R:
  ## "spearman_cor_safe() flags a zero-standard-deviation correlation" asserts the
  ## zero-SD branch, the ties branch (which must NOT be labelled zero-SD), and the
  ## clean branch.
  list(id = "fix13_pval_warn", commit = "f05c855", delta = 13,
       scenario = "S8", path = "^S8/s8_gradient",
       class = "EXPECTED_FIX", optional = TRUE,
       rationale = paste("pval_warn is initialised to NA up front and set only on cor.test's",
                         "zero-standard-deviation warning. Not reachable through STgradient()",
                         "from the fixture (see the comment above); pinned by a unit test on",
                         "spearman_cor_safe() instead.")),

  ## ---- order-induced value changes --------------------------------------
  ## These are NOT order-only and NOT regressions. STenrich's permutation null is
  ## set.seed(seed) then computeSubsampleSums(coords_mat, ...), and the C++ draws
  ## Rcpp::sample() over ROW POSITIONS, not names. Reordering the spots therefore
  ## makes the same seeded draws select different spots, so the null distribution
  ## and every p-value derived from it genuinely move. NEWS's claim that no
  ## expression value changed does not hold for these fields.
  list(id = "order_induced_stenrich_p", commit = "f278bd0", delta = 1,
       scenario = NA, path = "(p_value|adj_p_value|adj_p_val)$",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("Permutation null is drawn over row positions of coords_mat",
                         "(src/spatialge_helpers.cpp Rcpp::sample), so a spot reordering",
                         "changes which spots each seeded draw selects.")),
  list(id = "order_induced_cluster_labels", commit = "f278bd0", delta = 1,
       scenario = NA, path = "columns/stclust_",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("hclust() on the fixture's tie-saturated regular grid breaks ties by",
                         "input row order, so cluster labels permute even when the partition",
                         "is identical.")),
  list(id = "order_induced_vst", commit = "f278bd0", delta = 1,
       scenario = NA, path = "columns/vst\\.",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("Reductions over spots (colSums, var, the Rcpp VST) change summation",
                         "order under a column permutation, so the last bits of the",
                         "variance-stabilised statistics move.")),

  ## Two more instances of RNG-over-positions, both found by this comparison and
  ## both absent from NEWS. Same root cause as the STenrich permutation null: a
  ## seeded draw indexes ROW POSITIONS, so reordering the rows changes what the
  ## draw selects even though the seed is fixed.
  ## Tie-breaking amplifiers. Both were PROVEN order-induced by a single-version
  ## experiment (dev/xver/30_order_probe.R): permuting spot order within one build
  ## of alex-refactor reproduces the same magnitude of change, so neither is a
  ## regression -- the ordering fix merely changed the order these paths were
  ## always sensitive to. Neither is mentioned in NEWS.
  list(id = "order_induced_sthet_k", commit = "f278bd0", delta = 1,
       scenario = NA, path = "/s13_sthet_k/",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("SThet(k=) builds a k-nearest-neighbour graph. On a regular grid the",
                         "k-th distance is tied among several spots, and knearneigh breaks that",
                         "tie by position -- so reordering yields a different neighbour graph and",
                         "a materially different Moran's I (up to 0.018 on the fixture).",
                         "The default distance-based listw is unaffected (~1e-15).")),
  list(id = "order_induced_gradient_avg", commit = "f278bd0", delta = 1,
       scenario = NA, path = "/s13_gradient_avg/.*columns/avg_spearman",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("distsumm='avg' averages each spot's distances to the reference domain,",
                         "and on a regular grid many of those averages are exactly tied. A",
                         "summation-order change of ~1e-16 flips a tie, which flips a Spearman",
                         "rank, which moves the correlation by ~0.004. distsumm='min' does not",
                         "produce the same ties and is far less affected.")),

  list(id = "order_induced_plotcounts_subset", commit = "f278bd0", delta = 1,
       scenario = NA, path = "/[^/]*plot_counts[^/]*/value/.*/built/",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("plot_counts() does set.seed(subset_seed) then",
                         "df[sample(1:nrow(df), nrow(df)*distrib_subset), ] (R/plot_counts.R:88",
                         "and :99), so the seeded draw picks row POSITIONS. Reordering the spots",
                         "makes the same seed select a different 50% subset, changing the",
                         "plotted density. Cosmetic, but it IS an output change.")),
  list(id = "order_induced_sina_jitter", commit = "f278bd0", delta = 1,
       scenario = NA, path = "/built/\\[2\\]/(columns/(x|x_diff)|row_signatures)$",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("ggforce::geom_sina assigns each point a horizontal offset from its",
                         "position within the group, and draws that jitter inside",
                         "ggplot_build() rather than at plot construction (already noted in",
                         "dev/baseline/compare_baseline.R). Only x/x_diff move; the y values",
                         "carry identical multisets per group.")),
  list(id = "order_induced_row_signatures", commit = "f278bd0", delta = 1,
       scenario = NA, path = "/row_signatures$",
       class = "ORDER_INDUCED_VALUE", optional = TRUE,
       rationale = paste("Whole-row signatures of built plot data inherit whichever column",
                         "differences the two entries above describe.")),

  ## ---- delta 11 also surfaces in the S7 sweep ---------------------------
  list(id = "fix11_null_rnacounts_s7", commit = "6954044", delta = 11,
       scenario = "S7", path = "^S7/s7_STlist_(noargs|null)/error_message$",
       class = "EXPECTED_FIX",
       main_pattern = "argument is of length zero",
       alex_pattern = "must be provided",
       rationale = paste("Same fix as S6, reached through the no-args sweep: the cryptic",
                         "'argument is of length zero' becomes a named spatialGE error.")),

  ## ---- undocumented but real: which raw R error surfaces on bad input ----
  ## Not a fix and not a regression. resolve_samples() (R/utils_validate.R:92) reads
  ## the slot with names(methods::slot(x, slot)) where main used x@slot inline, and
  ## those two produce DIFFERENT R errors for a non-S4 input:
  ##   @             -> "no applicable method for `@` applied to an object of class ..."
  ##   methods::slot -> "cannot get a slot (\"...\") from an object of type ..."
  ## Both versions still error, so nothing is broken -- but the message a user sees
  ## changed across roughly ten exported functions, and NEWS does not mention it.
  ## Matched on the message text, not just the path, so a genuine behavior change in
  ## the same sweep cannot hide behind this record.
  ## SUPERSEDED by fix16_assert_stlist_msg above: alex no longer produces
  ## "cannot get a slot" on this path, so this record matches nothing. Marked
  ## optional rather than deleted, because the two records together are the
  ## record that this message went raw -> raw -> informative.
  list(id = "undoc_slot_error_text", commit = "6a35680", delta = NA,
       scenario = NA, path = "/error_message$",
       class = "UNDOCUMENTED_CHANGE", optional = TRUE,
       main_pattern = "no applicable method for",
       alex_pattern = "cannot get a slot",
       rationale = paste("Raw R slot-access error text changed because the guard now reads the",
                         "slot via methods::slot() rather than @. Cosmetic, but user-visible,",
                         "and it reveals that neither version gives a friendly error here:",
                         "resolve_samples() runs before any assert_stlist() class check.")),

  ## ---- delta 2 also surfaces wherever a domain= STenrich is run ----------
  list(id = "fix02_stenrich_domain_pinned", commit = "3cd478f", delta = 2,
       scenario = "S11", path = "^S11/s11_enrich_pinned/",
       class = "EXPECTED_FIX",
       rationale = paste("The same row-names collision, reached with an annotation pinned to the",
                         "fixture's own domain rather than to STclust output.")),

  ## =====================================================================
  ## Records below cover work committed after this table was first written
  ## (everything from a35b3cb onward). Before they existed, the run reported
  ## 9,853 REGRESSION paths in the main batch and 5,773 in the extras -- the
  ## harness failing closed on changes nobody had recorded, not on defects.
  ## =====================================================================

  ## ---- @misc$sthet$listws is no longer built on the default path ---------
  ## The single largest source of differences in both batches: 8,883 of 8,923
  ## main-batch paths and 5,345 of 5,773 extras paths, every one of them
  ## kind=only_main. Scoped to misc/sthet so it cannot absorb a moran_i/geary_c
  ## change, which lives under gene_meta and is a different question entirely.
  list(id = "change01_sthet_listw_absent", commit = "e631ddb", delta = NA,
       scenario = NA, path = "value/misc/sthet/",
       class = "EXPECTED_CHANGE",
       rationale = paste("SThet(k=NULL) computes Moran's I and Geary's C directly from the",
                         "coordinates in C++ rather than materializing a distance-weighted",
                         "complete graph, so @misc$sthet$listws is not populated on that path.",
                         "NEWS states this. The motive is capacity, not correctness: the old",
                         "listw needed 1.2 GB for one 4,096-spot sample and ~181 GB",
                         "extrapolated to 50,000 cells, i.e. it could not run on Xenium or",
                         "CosMx at all. SThet(k=) still builds a real spdep object, which is",
                         "why s13_sthet_k does not match this record.")),

  ## ---- transform_data() now populates the VST column --------------------
  ## Value-guarded on the alex side: the record only matches when the new column
  ## is what appeared, so an unrelated gene_meta column change cannot hide here.
  list(id = "change02_vst_column", commit = "a35b3cb", delta = NA,
       scenario = NA, path = "/gene_meta/[^/]+/colnames$",
       class = "EXPECTED_CHANGE",
       alex_pattern = "vst\\.variance\\.standardized",
       rationale = paste("transform_data() computes vst.variance.standardized itself instead of",
                         "leaving it to whichever analysis ran first, so @gene_meta gains one",
                         "column (3 -> 4). Does not change which genes any function SELECTS:",
                         "STclust, STgradient, gene_interpolation and STdiff all computed VST",
                         "themselves before ranking, on both versions.")),

  ## ---- pseudobulk_samples() died on a zero-variance gene ----------------
  ## Directional and unambiguous: main errors, alex succeeds. Verified from the
  ## producer payloads in both S0 (cores=1) and S0b (cores=3) --
  ## main ok=FALSE "cannot rescale a constant/zero column to unit variance",
  ## alex ok=TRUE. dim_plot and heatmap fail only because they depend on it.
  list(id = "fix14_pbulk_zero_variance", commit = "d00b0bf", delta = NA,
       scenario = NA, path = "^S0b?/pseudobulk_(samples|dim_plot|heatmap)/",
       class = "EXPECTED_FIX",
       rationale = paste("prcomp(scale.=TRUE) on a gene with zero variance across samples",
                         "errored, taking the whole pseudobulk family with it. The fixture",
                         "contains such a gene, so every pseudobulk step errored on main.")),

  ## ---- samples= matching nothing is an error, not a silent result -------
  ## main built a complete STlist keyed by a sample that does not exist
  ## ("not_a_sample", and "a" matched as a substring); alex raises. The whole
  ## object therefore differs, which is why these scenarios contribute counts
  ## across counts/tr_counts/gene_meta/spatial_meta rather than one path.
  list(id = "fix15_samples_match_nothing", commit = "b538c38", delta = NA,
       scenario = "S14",
       path = "^S14/s14_SThet_(missing_name|substring|out_of_range)/",
       class = "EXPECTED_FIX",
       rationale = paste("A samples= value matching nothing was a silent no-op that returned",
                         "results a caller would read as real -- including gene_meta and",
                         "spatial_meta entries named after the non-existent sample. Now an",
                         "error, so a typo no longer looks like a run that found nothing.")),

  ## ---- filter_data() reports what it removed ----------------------------
  ## Scoped by the message text rather than by scenario: filter_data is reached
  ## from S2, S4, S12 and S13, and a scenario-scoped record missed three of them.
  ## Note the counts appear blank in the payload because reduce.R's xv_scrub()
  ## strips 3-7 digit numbers at end-of-string so data-dependent magnitudes do
  ## not read as differences -- the counts are present in the real message.
  list(id = "change03_filter_verbose", commit = "a35b3cb", delta = NA,
       scenario = NA, path = "/messages$",
       class = "EXPECTED_CHANGE",
       alex_pattern = "filter_data: (removed|nothing was filtered)",
       rationale = paste("Every count threshold defaults to 0, so a call supplying none of",
                         "them filtered nothing while looking like a successful filter.",
                         "verbose=TRUE now says so, or reports the counts removed.")),

  ## ---- assert_stlist() replaces a raw slot-access error ------------------
  ## Supersedes undoc_slot_error_text below, which recorded an INTERMEDIATE state
  ## of this same message. At the time that record was written, alex answered a
  ## non-STlist input with "cannot get a slot (...) from an object of type ..." --
  ## a different raw R error, because the guard read the slot via methods::slot()
  ## instead of @. The error-catalog work replaced it with a real message, so the
  ## old record now matches nothing and reports MISSING_EXPECTED. Kept both, with
  ## the old one marked optional, because the progression is the evidence that
  ## this path went raw -> raw -> informative across two commits.
  list(id = "fix16_assert_stlist_msg", commit = "3a16a1c", delta = NA,
       scenario = "S7", path = "/error_message$",
       class = "EXPECTED_FIX",
       main_pattern = "no applicable method for",
       alex_pattern = "must be a STlist",
       rationale = paste("Passing NULL, no argument, or a non-STlist to any of ~15 exported",
                         "functions produced R's raw 'no applicable method for `@` applied to",
                         "an object of class \"NULL\"', which names nothing the caller can act",
                         "on. assert_stlist() now answers 'The input must be a STlist.'")),

  ## ---- which condition class a missing argument raises -------------------
  list(id = "undoc_missing_arg_class", commit = "6a35680", delta = NA,
       scenario = "S7", path = "/error_classes$",
       class = "UNDOCUMENTED_CHANGE",
       main_pattern = "evalError",
       alex_pattern = "getvarError",
       rationale = paste("For a missing required argument, the R condition class changed from",
                         "evalError to getvarError, because the value is now first touched by",
                         "a validation helper rather than inside the original expression.",
                         "Both still subclass missingArgError/error, so tryCatch on those is",
                         "unaffected; code matching the concrete class is not.")),

  ## ---- get_gene_meta() column ORDER changed -----------------------------
  ## Both versions return the same seven columns; the order differs, because
  ## transform_data() now inserts vst.variance.standardized before SThet appends
  ## moran_i/geary_c, where previously the analysis appended VST last. Recorded as
  ## undocumented: fidelity rule 8 treats column position as API, and NEWS
  ## mentions the new column but not that it moved the existing ones.
  list(id = "undoc_gene_meta_col_order", commit = "a35b3cb", delta = NA,
       scenario = NA, path = "^S0b?/get_gene_meta/value/colnames$",
       class = "UNDOCUMENTED_CHANGE",
       rationale = paste("main: sample, gene, gene_mean, gene_stdevs, moran_i, geary_c,",
                         "vst.variance.standardized. alex: the same seven with",
                         "vst.variance.standardized third from last instead of last. Content",
                         "is identical; anything selecting by position rather than name sees",
                         "a different column.")),

  ## ---- validation order: the message named the wrong argument ------------
  list(id = "fix17_stplot_interp_validation_order", commit = "3a16a1c", delta = NA,
       scenario = "S7", path = "^S7/s7_STplot_interpolation_",
       class = "EXPECTED_FIX",
       main_pattern = "enter one or more genes",
       alex_pattern = "must be a STlist",
       rationale = paste("STplot_interpolation() tested genes before validating x, so passing",
                         "NULL, no argument or a non-STlist was reported as 'Please, enter one",
                         "or more genes to plot.' -- naming an argument the caller had nothing",
                         "wrong with. The input is now validated first.")),

  ## ---- the VST loess now emits its numerical warnings from transform_data -
  ## Scenario-scoped to S9: only that scenario's data makes the loess warn, so a
  ## scenario=NA record would report MISSING_EXPECTED on every other batch.
  list(id = "change05_vst_loess_warnings", commit = "a35b3cb", delta = NA,
       scenario = "S9", path = "/warnings$",
       class = "EXPECTED_CHANGE",
       alex_pattern = paste("NaNs produced|near singularities|neighborhood radius",
                            "|pseudoinverse used|reciprocal condition number", sep = ""),
       rationale = paste("Moving vst.variance.standardized into transform_data() moved the",
                         "stats::loess fit there too, so its numerical warnings now surface",
                         "from the transform rather than from whichever analysis ran first.",
                         "The warnings are the loess reporting a marginal fit on this data;",
                         "they occurred on main as well, just later and only if an analysis",
                         "that needed VST was run at all.")),

  ## ---- pseudobulk zero-variance drop, as a warning on real data ---------
  ## Scenario-scoped to S12: on the fixture main errors outright (fix14 above),
  ## so the warning form is only reachable where the gene count permits it.
  list(id = "fix14_pbulk_zero_variance_warn", commit = "d00b0bf", delta = NA,
       scenario = "S12", path = "/warnings$",
       class = "EXPECTED_FIX",
       alex_pattern = "zero variance across samples",
       rationale = paste("Same fix as fix14, seen from the other side: where the gene count",
                         "permits, the degenerate genes are dropped with a warning instead of",
                         "erroring. On the real Visium data this is a warning rather than the",
                         "hard failure main produced.")),

  ## ---- filter_data() removing every gene leaked a dplyr error -----------
  ## Found by this harness: s13_filter_gene_minpct errored on BOTH versions, which
  ## the vacuous-coverage report surfaced once the intentional error-path steps
  ## were declared away. A 0-row dgCMatrix carries rownames = NULL rather than
  ## character(0), so sample_genes() returns NULL, tibble() DROPS the `gene`
  ## column it is handed, and the @gene_meta left_join(by='gene') fails with
  ## "Join columns in `x` must be present in the data" -- naming neither
  ## filter_data nor the threshold. Reachable from ordinary arguments
  ## (gene_minreads=1e9, gene_minspots=1e6) as well as from the pct scale
  ## confusion, and present in 1.2.2, so this is a long-standing crash rather
  ## than a refactor regression.
  list(id = "fix18_filter_all_genes_removed", commit = NA, delta = NA,
       scenario = "S13", path = "^S13/s13_filter_gene_minpct_all_removed/",
       class = "EXPECTED_FIX",
       main_pattern = "Join columns",
       alex_pattern = "removed every gene",
       rationale = paste("filter_data() now raises error0098, naming the sample and stating",
                         "that the pct thresholds are proportions in [0,1] rather than",
                         "percentages -- which is the mistake the old @param text invited by",
                         "calling them percentages. Pinned in",
                         "tests/testthat/test-filter_data.R and guarded by the leaked vector",
                         "in test-graceful-errors.R.")),

  ## The condition CLASS changes with the error source: dplyr raises an
  ## rlang_error, raise_err() uses stop(..., call.=FALSE) and so raises a
  ## simpleError. Given its own record rather than relaxing fix18's value guards,
  ## which deliberately refused to absorb a leaf whose values they did not
  ## describe -- that refusal is the guard working, not an inconvenience.
  list(id = "fix18_filter_all_genes_removed_class", commit = NA, delta = NA,
       scenario = "S13", path = "^S13/s13_filter_gene_minpct_all_removed/error_classes$",
       class = "EXPECTED_FIX",
       main_pattern = "rlang_error",
       alex_pattern = "simpleError",
       rationale = paste("Same fix as fix18, seen in the condition class: the dplyr failure was",
                         "an rlang_error, and the catalog error is a simpleError. Code",
                         "tryCatch-ing on rlang_error specifically would need updating; code",
                         "catching `error` is unaffected.")),

  ## ---- declared non-comparable ------------------------------------------
  list(id = "noncomparable_krige_call", commit = NA, delta = NA,
       scenario = NA, path = "gene_krige/.*/(call|fields)$",
       class = "NON_COMPARABLE", optional = TRUE,
       rationale = "gstat/krige objects carry calls and environments; only fitted parameters are a contract."),
  list(id = "noncomparable_elapsed", commit = NA, delta = NA,
       scenario = NA, path = "/elapsed$",
       class = "NON_COMPARABLE", optional = TRUE,
       rationale = "Wall clock. Fidelity rule 10: timing output is inherently non-reproducible.")
)

## Paths whose leaves are order-sensitive by algorithm rather than by layout, so a
## canonical sort cannot make them comparable and a difference is not a bug.
XV_ORDER_SENSITIVE <- c("(p_value|adj_p_value|adj_p_val)$",
                        "columns/stclust_",
                        "columns/vst\\.")

## A record matches a path only if its regex matches AND, when the record declares
## main_pattern/alex_pattern, the observed values match those too. Value-aware
## matching keeps a broad path regex from silently absorbing a genuine behavior
## change that happens to sit at the same kind of path.
xv_match_expect <- function(path, scenario, main_val = NULL, alex_val = NULL){
  for(rec in XV_EXPECT){
    if(!is.na(rec$scenario) && !identical(rec$scenario, scenario)) next
    if(!grepl(rec$path, path, perl = TRUE)) next
    if(!is.null(rec$main_pattern)){
      if(is.null(main_val) || !any(grepl(rec$main_pattern, as.character(main_val)))) next
    }
    if(!is.null(rec$alex_pattern)){
      if(is.null(alex_val) || !any(grepl(rec$alex_pattern, as.character(alex_val)))) next
    }
    return(rec)
  }
  NULL
}


## ---------------------------------------------------------------------------
## Steps that are SUPPOSED to error on both versions.
##
## A step erroring identically on both sides compares only an error string, and
## the vacuous-coverage report counts that as debt. But most of those steps exist
## precisely to assert that bad input errors -- the error IS the output being
## compared -- so an undifferentiated count buries the ones that matter. The
## first run of this table reported 86 such steps; 85 were intentional and 1 was
## a live crash in filter_data() (see below), which is exactly the signal the
## count was drowning.
##
## Declaring the intentional ones means an UNDECLARED error-on-both step is now
## a finding rather than a line in a list of 86. Add to this only when the step's
## purpose is to assert an error, and say why.
XV_ERROR_PATH <- c(
  ## S6, S7: the bad-input sweeps. Every step passes no argument, NULL, or a
  ## wrong-typed object to an exported function and asserts it errors.
  "^S6/s6_(null_rnacounts|null_both|null_spotcoords|no_args)$",
  "^S7/s7_[A-Za-z_]+_(noargs|null|wrong)$",

  ## S14: the samples= resolution sweep. For the plotting functions both versions
  ## error on a name that matches nothing, a substring, or an out-of-range index.
  ## (SThet is NOT here: main returned fabricated results, which is fix15.)
  "^S14/s14_(plot_counts|plot_image|STplot|distribution)_(substring|missing_name|out_of_range)$",

  ## S1: a samples file whose own columns point at per-sample counts/coords
  ## paths, with no rnacounts/spotcoords given. This format has never been
  ## supported -- process_sample_names_from_file() was dead code and was removed
  ## in Phase 4 -- and tests/testthat/test-fixtures.R pins that STlist() rejects
  ## it. Erroring on both versions is the correct outcome.
  "^S1/STlist_samplesfile_matrices$",

  ## S13: deliberately out-of-range gene_minpct, kept as a step to assert that a
  ## filter removing every gene errors informatively rather than leaking dplyr.
  "^S13/s13_filter_gene_minpct_all_removed$"
)

xv_is_error_path <- function(step_id){
  any(vapply(XV_ERROR_PATH, function(p) grepl(p, step_id, perl = TRUE), logical(1)))
}
