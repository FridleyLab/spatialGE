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
  list(id = "fix13_pval_warn", commit = "f05c855", delta = 13,
       scenario = "S8", path = "^S8/s8_gradient",
       class = "EXPECTED_FIX", optional = TRUE,
       rationale = paste("pval_warn is now initialised to NA up front rather than only inside",
                         "if(nrow(df_gene_range) > 1). Reachable only when a gene's",
                         "post-filter sample size is <= 1; if main does not error here the",
                         "delta is unreachable and stays unverified.")),

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
  list(id = "undoc_slot_error_text", commit = "6a35680", delta = NA,
       scenario = NA, path = "/error_message$",
       class = "UNDOCUMENTED_CHANGE",
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
