## The scenario matrix. Every scenario is version-neutral: it calls only exported
## functions, all of which have verified-identical formals on both refs, so the
## two producers run the same call sequence and any difference is the packages'.
##
## S0 is the anchor -- dev/baseline/pipeline.R's bl_pipeline() unmodified, which
## is also what tests/golden/ covers. S1..S14 exist because that anchor reaches
## only 2 of the 15 deliberate behavior changes: it uses one input format, one
## parameter set per function, three samples, and never inspects a warning.
##
## Each scenario returns a named list of xv_capture() records.

XV_SEED <- 20260115L

## ---------------------------------------------------------------------------
## S0 / S0b -- the anchor
## ---------------------------------------------------------------------------

xv_s0 <- function(cores){
  ## bl_pipeline() is sourced, not copied, so this is byte-for-byte the sequence
  ## the committed golden gate runs. bl_capture() has no warnings/messages
  ## fields, so normalise its records into xv_capture()'s shape.
  res <- bl_pipeline(cores = cores, verbose = FALSE)
  lapply(res, function(r){
    r$warnings <- character(0)
    r$messages <- character(0)
    r
  })
}

## ---------------------------------------------------------------------------
## S1 -- every input format
##
## Covers delta 3 (Xenium MEX import), delta 6 (xlsx samples-file detection) and
## delta 8 (the undefined h5_file in the import warning path). The anchor only
## ever builds an STlist from in-memory data frames, so all three import fixes are
## unreachable from it.
## ---------------------------------------------------------------------------

xv_s1 <- function(root){
  fx <- fixture_materialize_all(root)
  variants <- setdiff(names(fx), "base")
  out <- list()
  for(v in variants){
    spec <- fx[[v]]
    out[[paste0("STlist_", v)]] <- xv_capture({
      args <- list(cores = 1, verbose = FALSE)
      for(a in intersect(c("rnacounts", "spotcoords", "samples"), names(spec))){
        args[[a]] <- spec[[a]]
      }
      do.call(STlist, args)
    }, seed = XV_SEED)
  }
  out
}

## ---------------------------------------------------------------------------
## S2 -- rm_genes_expr reaching @gene_meta (delta 5)
##
## On main, filter_data()'s gene-expression stage indexes @gene_meta with a
## malformed 3-subscript index and matches against rownames() of a tibble, which
## are "1".."200" -- so the grepl matches nothing and the stage is a silent no-op.
## The public path is still correct end to end because a later stage re-syncs
## @gene_meta from @tr_counts, which is why this was deliberately left out of
## NEWS. Compared here anyway: "already correct downstream" is a claim to verify.
## ---------------------------------------------------------------------------

xv_s2 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  out <- list()
  out$s2_stlist <- xv_capture(
    STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE),
    seed = XV_SEED)
  if(!isTRUE(out$s2_stlist$ok)) return(out)
  st <- out$s2_stlist$value
  out$s2_tr <- xv_capture(transform_data(st, method = "log", cores = 1), seed = XV_SEED)
  if(!isTRUE(out$s2_tr$ok)) return(out)
  st2 <- out$s2_tr$value
  out$s2_rm_genes_expr <- xv_capture(filter_data(st2, rm_genes_expr = "^Grad"), seed = XV_SEED)
  out$s2_rm_genes      <- xv_capture(filter_data(st2, rm_genes = c("Gene1", "Gene2")), seed = XV_SEED)
  out
}

## ---------------------------------------------------------------------------
## S3 -- STclust's clobbered verbose argument (delta 4)
##
## The observable is process stdout, not a return value: sge_worker_msg keeps
## system(sprintf('echo "%s"', msg)) verbatim (fidelity rule 11), which writes to
## the process's own stdout where capture.output cannot see it. The producer
## redirects whole-process stdout to a file and the comparator diffs transcripts.
## cores=1 only -- output emitted from a fork is not attributable to a step.
## ---------------------------------------------------------------------------

xv_s3 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st, method = "log", cores = 1)
  out <- list()
  for(vb in list(TRUE, FALSE)){
    tag <- paste0("s3_stclust_verbose", vb)
    cat(sprintf("<<<XVMARK %s>>>\n", tag))
    out[[tag]] <- xv_capture(
      STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
              cores = 1, verbose = vb),
      seed = XV_SEED)
    cat(sprintf("<<<XVMARK /%s>>>\n", tag))
  }
  out
}

## ---------------------------------------------------------------------------
## S4 -- the pseudobulk chain on >= 4 samples (delta 9)
##
## pseudobulk_samples() refuses fewer than four samples, and the committed fixture
## has three, so on BOTH refs the anchor records pseudobulk_samples,
## pseudobulk_dim_plot and pseudobulk_heatmap as errors. That identical-error
## result reads as "no difference" while testing no code at all -- which is
## exactly how the inverted coefficient-of-variation fix stayed unverified.
##
## Shadowing fixture_sample_names() only APPENDS rpois/sample.int draws, so
## alpha/bravo/charlie stay bit-identical to the 3-sample fixture and the extra
## samples are new draws on the end.
## ---------------------------------------------------------------------------

xv_s4 <- function(){
  old <- fixture_sample_names
  on.exit(assign("fixture_sample_names", old, envir = globalenv()), add = TRUE)
  assign("fixture_sample_names",
         function() c("alpha", "bravo", "charlie", "delta", "echo"),
         envir = globalenv())

  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  out <- list()
  out$s4_stlist <- xv_capture(
    STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE),
    seed = XV_SEED)
  if(!isTRUE(out$s4_stlist$ok)) return(out)
  ## The fixture deliberately contains an all-zero gene and a constant gene.
  ## pseudobulk_samples() scales genes to unit variance, which fails outright on a
  ## zero-variance column ("cannot rescale a constant/zero column to unit
  ## variance"), so without dropping them the scenario never reaches the heatmap
  ## and the CV fix stays unverified. dev/baseline/pipeline.R drops the same two
  ## genes before SCT for the same class of reason.
  out$s4_drop_degenerate <- xv_capture(
    filter_data(out$s4_stlist$value, rm_genes = c("ZeroGene", "ConstGene")), seed = XV_SEED)
  if(!isTRUE(out$s4_drop_degenerate$ok)) return(out)
  out$s4_tr <- xv_capture(transform_data(out$s4_drop_degenerate$value, method = "log", cores = 1),
                          seed = XV_SEED)
  if(!isTRUE(out$s4_tr$ok)) return(out)
  st <- out$s4_tr$value

  ## calc_umap=FALSE is mandatory: uwot::umap is unseeded and multithreaded
  ## (fidelity rule 1), so a UMAP would be noise rather than evidence.
  out$s4_pbulk <- xv_capture(pseudobulk_samples(st, calc_umap = FALSE), seed = XV_SEED)
  if(!isTRUE(out$s4_pbulk$ok)) return(out)
  pb <- out$s4_pbulk$value

  out$s4_dim_plot <- xv_capture(pseudobulk_dim_plot(pb), seed = XV_SEED)
  out$s4_heatmap  <- xv_capture(pseudobulk_heatmap(pb, hm_display_genes = 5), seed = XV_SEED)
  out$s4_heatmap_30 <- xv_capture(pseudobulk_heatmap(pb, hm_display_genes = 30), seed = XV_SEED)

  ## The oracle's raw material: the matrix the gene ranking is computed from. If
  ## this differs between versions the CV oracle is void, so it must be compared
  ## on its own rather than inferred from the heatmap.
  out$s4_scaled_mtx <- xv_capture(pb@misc[["scaled_pbulk_mtx"]], seed = XV_SEED)
  out
}

## ---------------------------------------------------------------------------
## S5 -- STdiff_volcano's numeric samples= (delta 10)
##
## main tests is.numeric(x) where it means is.numeric(samples), so a numeric
## samples= is silently ignored and the plot is labelled from the wrong branch.
## ---------------------------------------------------------------------------

xv_s5 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st, method = "log", cores = 1)
  st <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                cores = 1, verbose = FALSE)
  out <- list()
  out$s5_stdiff <- xv_capture(
    STdiff(st, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 100,
           test_type = "t_test", cores = 1), seed = XV_SEED)
  if(!isTRUE(out$s5_stdiff$ok)) return(out)
  de <- out$s5_stdiff$value
  out$s5_volcano_numeric <- xv_capture(STdiff_volcano(de, samples = 1), seed = XV_SEED)
  out$s5_volcano_name    <- xv_capture(STdiff_volcano(de, samples = "alpha"), seed = XV_SEED)
  out$s5_volcano_default <- xv_capture(STdiff_volcano(de), seed = XV_SEED)
  out
}

## ---------------------------------------------------------------------------
## S6 -- STlist() called without rnacounts (delta 11)
## ---------------------------------------------------------------------------

xv_s6 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  list(
    s6_null_rnacounts  = xv_capture(STlist(rnacounts = NULL, spotcoords = fx$spotcoords), seed = XV_SEED),
    s6_null_both       = xv_capture(STlist(rnacounts = NULL, spotcoords = NULL), seed = XV_SEED),
    s6_null_spotcoords = xv_capture(STlist(rnacounts = fx$rnacounts, spotcoords = NULL), seed = XV_SEED),
    s6_no_args         = xv_capture(STlist(), seed = XV_SEED)
  )
}

## ---------------------------------------------------------------------------
## S7 -- error-message sweep across every export
##
## Zero coverage today. The refactor replaced hand-written guards at ~30 call
## sites with assert_stlist / assert_transformed / resolve_samples, deliberately
## preserving several pre-existing quirks (including an unanchored substring
## grep() for sample names). A wrong flag at one call site changes only the error
## a user sees, which no current test would notice.
## ---------------------------------------------------------------------------

xv_s7 <- function(){
  exports <- c("STlist", "filter_data", "transform_data", "SThet", "compare_SThet",
               "STclust", "STdiff", "STdiff_volcano", "STenrich", "STgradient",
               "gene_interpolation", "STplot_interpolation", "pseudobulk_samples",
               "pseudobulk_dim_plot", "pseudobulk_heatmap", "distribution_plots",
               "plot_counts", "plot_image", "STplot", "summarize_STlist",
               "get_gene_meta", "load_images", "spatial_metadata", "tissue_names")
  out <- list()
  for(fn in exports){
    f <- get(fn, envir = asNamespace("spatialGE"))
    out[[paste0("s7_", fn, "_noargs")]] <- xv_capture(f(), seed = XV_SEED)
    out[[paste0("s7_", fn, "_null")]]   <- xv_capture(f(NULL), seed = XV_SEED)
    out[[paste0("s7_", fn, "_wrong")]]  <- xv_capture(f(42L), seed = XV_SEED)
  }
  out
}

## ---------------------------------------------------------------------------
## S8 -- STgradient's pval_warn state leak (delta 13)
##
## In the default robust=TRUE branch, alex initialises pval_warn to NA up front
## while main only assigns it inside `if(nrow(df_gene_range) > 1)`, and main's
## rm(list=grep(...)) does not match the name. So with at most one non-reference
## spot, main should error with "object 'pval_warn' not found" and alex should
## return a row. If main does not error, the delta is unreachable and the report
## must say so rather than claiming coverage.
## ---------------------------------------------------------------------------

xv_s8 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st, method = "log", cores = 1)

  ## Inject a deterministic annotation directly. The STlist class definition is
  ## byte-identical on both refs, so this is version-neutral.
  for(s in names(st@spatial_meta)){
    d <- st@spatial_meta[[s]]
    lab <- rep("R", nrow(d))
    lab[1] <- "X"           # exactly one non-reference spot
    d[["xv_annot"]] <- lab
    st@spatial_meta[[s]] <- d
  }
  list(
    s8_gradient_1spot = xv_capture(
      STgradient(st, annot = "xv_annot", ref = "R", topgenes = 20,
                 cores = 1, verbose = FALSE), seed = XV_SEED),
    s8_gradient_ref_X = xv_capture(
      STgradient(st, annot = "xv_annot", ref = "X", topgenes = 20,
                 cores = 1, verbose = FALSE), seed = XV_SEED)
  )
}

## ---------------------------------------------------------------------------
## S11 -- pinned annotation
##
## STclust's cluster labels are order-sensitive: hclust() on the fixture's regular
## 12x10 integer grid has a distance matrix saturated with exact ties, and tie
## merge order follows input row order. The anchor feeds those labels straight
## into STenrich/STgradient/STdiff as annot=/domain=/ref=, so if they permute, the
## two refs are being asked different questions and their results are not
## comparable at all. Pinning the annotation to the fixture's own known domain
## isolates each function's fidelity from that.
## ---------------------------------------------------------------------------

xv_s11 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st, method = "log", cores = 1)

  ## Domain A/B is a property of the coordinates (xpos <= 6), so it is identical
  ## on both refs regardless of row order.
  for(s in names(st@spatial_meta)){
    d <- st@spatial_meta[[s]]
    d[["xv_domain"]] <- ifelse(d[["xpos"]] <= 6, "A", "B")
    st@spatial_meta[[s]] <- d
  }
  gs <- list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
             domain   = c("DomainGeneA", "DomainGeneB"))
  list(
    s11_enrich_pinned = xv_capture(
      STenrich(st, gene_sets = gs, annot = "xv_domain", domain = "A",
               reps = 20, min_units = 1, min_genes = 1, cores = 1, verbose = FALSE),
      seed = XV_SEED),
    s11_enrich_unpinned = xv_capture(
      STenrich(st, gene_sets = gs, reps = 20, min_units = 1, min_genes = 1,
               cores = 1, verbose = FALSE), seed = XV_SEED),
    s11_gradient_pinned = xv_capture(
      STgradient(st, annot = "xv_domain", ref = "A", topgenes = 100,
                 cores = 1, verbose = FALSE), seed = XV_SEED),
    s11_stdiff_pinned = xv_capture(
      STdiff(st, annot = "xv_domain", test_type = "t_test", topgenes = 100, cores = 1),
      seed = XV_SEED)
  )
}

## ---------------------------------------------------------------------------
## S13 -- non-default argument sweep
##
## The anchor uses exactly one parameter set per function, so whole branches are
## never executed: STclust(ks=<numeric>) means get_hier_clusters_ks() is never
## called at all, and the untaken branches in STdiff_prepare/_fit/_assemble,
## filter_data_stages and stlist_readers are precisely where a refactor bug would
## sit undetected.
## ---------------------------------------------------------------------------

xv_s13 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st0 <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st0, method = "log", cores = 1)

  out <- list()
  ## STclust: the numeric-ks branch, plus deepSplit variants of the dtc branch.
  out$s13_stclust_ks_numeric <- xv_capture(
    STclust(st, ws = 0.025, ks = 2:3, topgenes = 100, cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_stclust_deepsplit <- xv_capture(
    STclust(st, ws = 0.025, ks = "dtc", deepSplit = TRUE, topgenes = 100,
            cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_stclust_ws0 <- xv_capture(
    STclust(st, ws = 0, ks = "dtc", deepSplit = FALSE, topgenes = 100,
            cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_stclust_linkage <- xv_capture(
    STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, linkage = "average",
            dist_metric = "manhattan", topgenes = 100, cores = 1, verbose = FALSE), seed = XV_SEED)

  ## SThet: the Geary branch and an explicit k.
  het <- c("Gene1", "GradGene1", "DomainGeneA", "ZeroGene", "ConstGene")
  out$s13_sthet_geary <- xv_capture(
    SThet(st, genes = het, method = "geary", cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_sthet_k <- xv_capture(
    SThet(st, genes = het, method = "moran", k = 4, cores = 1, verbose = FALSE), seed = XV_SEED)

  ## filter_data: each numeric stage, one at a time, so a stage that stops
  ## filtering is attributable rather than masked by the others.
  out$s13_filter_spot_minreads <- xv_capture(filter_data(st, spot_minreads = 50), seed = XV_SEED)
  out$s13_filter_spot_maxreads <- xv_capture(filter_data(st, spot_maxreads = 900), seed = XV_SEED)
  out$s13_filter_gene_minspots <- xv_capture(filter_data(st, gene_minspots = 5), seed = XV_SEED)
  out$s13_filter_gene_minpct   <- xv_capture(filter_data(st, gene_minpct = 10), seed = XV_SEED)
  out$s13_filter_spot_mingenes <- xv_capture(filter_data(st, spot_mingenes = 20), seed = XV_SEED)

  ## STgradient: the non-default statistic and outlier/robust branches.
  stc <- STclust(st, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100,
                 cores = 1, verbose = FALSE)
  for(s in names(stc@spatial_meta)){
    d <- stc@spatial_meta[[s]]
    d[["xv_domain"]] <- ifelse(d[["xpos"]] <= 6, "A", "B")
    stc@spatial_meta[[s]] <- d
  }
  out$s13_gradient_avg <- xv_capture(
    STgradient(stc, annot = "xv_domain", ref = "A", distsumm = "avg", topgenes = 50,
               cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_gradient_notrobust <- xv_capture(
    STgradient(stc, annot = "xv_domain", ref = "A", robust = FALSE, topgenes = 50,
               cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s13_gradient_outrm_log <- xv_capture(
    STgradient(stc, annot = "xv_domain", ref = "A", out_rm = TRUE, log_dist = TRUE,
               topgenes = 50, cores = 1, verbose = FALSE), seed = XV_SEED)

  ## STplot / plotting with plot_meta rather than genes.
  out$s13_stplot_meta <- xv_capture(
    STplot(stc, plot_meta = "xv_domain"), seed = XV_SEED)
  out$s13_plot_counts_violin <- xv_capture(
    plot_counts(st, data_type = "tr", plot_type = "violin"), seed = XV_SEED)
  out$s13_distribution_box <- xv_capture(
    distribution_plots(st, genes = c("Gene1", "GradGene1"), data_type = "tr",
                       plot_type = "boxplot"), seed = XV_SEED)
  out
}

## ---------------------------------------------------------------------------
## S15 -- a multi-group .h5, which is the only way to reach delta 8
##
## import_visium_h5()/import_xenium_h5() warn when the file contains more than one
## HDF5 group. On main that warning() interpolates an undefined `h5_file`, so
## RAISING THE WARNING ITSELF ERRORS -- the import dies where it meant to carry on.
## The committed fixture writes a well-formed single-group file, so S1's h5 variants
## succeed on both refs with zero warnings and never touch this branch: without this
## scenario delta 8 is simply unreachable, not verified.
##
## Reached by adding a second group to a copy of the fixture's own .h5.
## ---------------------------------------------------------------------------

xv_s15 <- function(root){
  if(!requireNamespace("hdf5r", quietly = TRUE)){
    return(list(s15_unavailable = xv_capture(stop("hdf5r not installed"), seed = XV_SEED)))
  }
  base <- fixture_base_data()
  src <- file.path(root, "h5_single")
  fx <- fixture_write_visium_h5(base, src)

  ## Copy, then add an extra group so length(names(h5)) > 1.
  dst <- file.path(root, "h5_multigroup")
  unlink(dst, recursive = TRUE)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
  for(s in base$samples){
    fp <- file.path(dst, s, "filtered_feature_bc_matrix.h5")
    h5 <- hdf5r::H5File$new(fp, mode = "r+")
    ## Named so it sorts AFTER 'matrix': the reader takes names()[1], so the real
    ## group must still be chosen and only the warning branch is new.
    g <- h5$create_group("zz_extra")
    g[["note"]] <- "second group to trigger the multi-group warning"
    h5$close_all()
  }
  list(
    s15_h5_single = xv_capture(
      STlist(rnacounts = file.path(src, base$samples), samples = base$samples), seed = XV_SEED),
    s15_h5_multigroup = xv_capture(
      STlist(rnacounts = file.path(dst, base$samples), samples = base$samples), seed = XV_SEED)
  )
}

## ---------------------------------------------------------------------------
## S9 -- STdiff's spatial mixed-model path (delta 12, partial)
##
## test_type='mm' is STdiff's DEFAULT and has zero coverage: the anchor pins
## test_type='t_test' because a realistic mm fit costs ~21s and ~27GB per cluster.
## So the flagship function's default behaviour is currently unverified in both
## the golden gate and the test suite. This runs a deliberately tiny problem to
## check the numeric path agrees at all.
##
## It does NOT verify the mclapply-poisoning guard: that needs a spaMM fit that
## fails while sharing a preschedule chunk with fits that succeed, and
## ARCHITECTURE.md fidelity rule 2 states the fixture policy is the opposite
## ("design fixture data so every fit converges"). Reported as a gap.
## ---------------------------------------------------------------------------

xv_s9 <- function(){
  set.seed(XV_SEED)
  genes <- c(paste0("G", 1:12), "DomA", "DomB")
  grid <- expand.grid(xpos = 1:6, ypos = 1:5)
  n <- nrow(grid)
  mk <- function(){
    m <- matrix(0L, nrow = length(genes), ncol = n, dimnames = list(genes, NULL))
    for(i in seq_along(genes)){
      g <- genes[i]
      lam <- if(g == "DomA") ifelse(grid$xpos <= 3, 12, 1)
             else if(g == "DomB") ifelse(grid$xpos > 3, 12, 1)
             else 4
      m[i, ] <- stats::rpois(n, lam)
    }
    m
  }
  counts <- list(); coords <- list()
  for(s in c("s1", "s2")){
    bc <- paste0(s, "_c", sprintf("%02d", seq_len(n)))
    m <- mk(); colnames(m) <- bc
    counts[[s]] <- data.frame(gene = genes, m, check.names = FALSE, stringsAsFactors = FALSE)
    coords[[s]] <- data.frame(libname = bc, ypos = grid$ypos, xpos = grid$xpos,
                              stringsAsFactors = FALSE)
  }
  out <- list()
  out$s9_stlist <- xv_capture(STlist(rnacounts = counts, spotcoords = coords,
                                     cores = 1, verbose = FALSE), seed = XV_SEED)
  if(!isTRUE(out$s9_stlist$ok)) return(out)
  out$s9_tr <- xv_capture(transform_data(out$s9_stlist$value, method = "log", cores = 1),
                          seed = XV_SEED)
  if(!isTRUE(out$s9_tr$ok)) return(out)
  st <- out$s9_tr$value
  for(s in names(st@spatial_meta)){
    d <- st@spatial_meta[[s]]
    d[["xv_domain"]] <- ifelse(d[["xpos"]] <= 3, "A", "B")
    st@spatial_meta[[s]] <- d
  }
  out$s9_stdiff_mm <- xv_capture(
    STdiff(st, annot = "xv_domain", test_type = "mm", topgenes = 8, sp_topgenes = 1,
           cores = 1, verbose = 0L), seed = XV_SEED)
  out$s9_stdiff_wilcox <- xv_capture(
    STdiff(st, annot = "xv_domain", test_type = "wilcoxon", topgenes = 8, cores = 1,
           verbose = 0L), seed = XV_SEED)
  out
}

## ---------------------------------------------------------------------------
## S10 -- cores= and mc.cores (delta 7)
##
## main accepted cores= in STlist() and transform_data() and then ignored it, so
## both always ran at the auto-detected count. The fix is real but invisible in
## any return value -- it changes fork count and wall clock only. What CAN be
## verified is the invariant that must hold either way: for a given version,
## varying cores must not change the result. Any difference here other than under
## sct (fidelity rule 1: forked children discard the RNG) is a regression.
## ---------------------------------------------------------------------------

xv_s10 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  out <- list()
  for(cr in c(1, 3)){
    tag <- paste0("s10_cores", cr)
    out[[paste0(tag, "_stlist")]] <- xv_capture(
      STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = cr, verbose = FALSE),
      seed = XV_SEED)
    if(!isTRUE(out[[paste0(tag, "_stlist")]]$ok)) next
    st <- out[[paste0(tag, "_stlist")]]$value
    out[[paste0(tag, "_log")]] <- xv_capture(
      transform_data(st, method = "log", cores = cr), seed = XV_SEED)
  }
  ## Same question with the option rather than the argument: main leaves several
  ## mclapply sites with no mc.cores at all, falling through to
  ## getOption("mc.cores", 2L) -- fidelity rule 3 keeps that verbatim on purpose.
  old <- getOption("mc.cores")
  on.exit(options(mc.cores = old), add = TRUE)
  for(mc in c(1L, 2L)){
    options(mc.cores = mc)
    out[[paste0("s10_mcopt", mc, "_stlist")]] <- xv_capture(
      STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE),
      seed = XV_SEED)
  }
  options(mc.cores = 1L)
  out
}

## ---------------------------------------------------------------------------
## S12 -- real Visium data (Thrane et al. melanoma)
##
## The only test of the ordering fix against genuinely unsorted real coordinate
## files: the fixture's shuffling is synthetic and only 'charlie' is shuffled at
## all. This dataset also has >= 4 samples, so it independently exercises the
## pseudobulk path that the 3-sample fixture cannot reach.
##
## A download failure is recorded as a failure, never swallowed into a skip -- a
## tryCatch that reports "skipped" is how a coverage gap becomes a pass.
## ---------------------------------------------------------------------------

xv_s12 <- function(root){
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  lk <- "https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download="
  zipfp <- file.path(root, "melanoma_thrane.zip")
  if(!file.exists(zipfp)){
    utils::download.file(lk, destfile = zipfp, mode = "wb", quiet = TRUE)
  }
  utils::unzip(zipfile = zipfp, exdir = root)
  d <- file.path(root, "melanoma_thrane")
  counts <- sort(list.files(d, full.names = TRUE, pattern = "counts"))
  coords <- sort(list.files(d, full.names = TRUE, pattern = "mapping"))
  clin   <- list.files(d, full.names = TRUE, pattern = "clinical")
  if(!length(counts) || !length(coords)) stop("Thrane data did not unpack as expected")

  out <- list()
  out$s12_stlist <- xv_capture(
    STlist(rnacounts = counts, spotcoords = coords, samples = clin), seed = XV_SEED)
  if(!isTRUE(out$s12_stlist$ok)) return(out)
  st <- out$s12_stlist$value
  out$s12_n_samples <- xv_capture(length(st@counts), seed = XV_SEED)
  out$s12_filter <- xv_capture(filter_data(st, spot_minreads = 100, gene_minspots = 5),
                               seed = XV_SEED)
  st1 <- if(isTRUE(out$s12_filter$ok)) out$s12_filter$value else st
  out$s12_log <- xv_capture(transform_data(st1, method = "log", cores = 1), seed = XV_SEED)
  if(!isTRUE(out$s12_log$ok)) return(out)
  st2 <- out$s12_log$value
  out$s12_sthet <- xv_capture(
    SThet(st2, genes = c("HLA-A", "CD74", "PMEL"), cores = 1, verbose = FALSE), seed = XV_SEED)
  out$s12_summary <- xv_capture(summarize_STlist(st2), seed = XV_SEED)
  ## >= 4 samples, so unlike the synthetic fixture this reaches the CV fix.
  out$s12_pbulk <- xv_capture(pseudobulk_samples(st2, calc_umap = FALSE), seed = XV_SEED)
  if(isTRUE(out$s12_pbulk$ok)){
    pb <- out$s12_pbulk$value
    out$s12_heatmap <- xv_capture(pseudobulk_heatmap(pb, hm_display_genes = 10), seed = XV_SEED)
    out$s12_scaled_mtx <- xv_capture(pb@misc[["scaled_pbulk_mtx"]], seed = XV_SEED)
    out$s12_dim_plot <- xv_capture(pseudobulk_dim_plot(pb), seed = XV_SEED)
  }
  out
}

## ---------------------------------------------------------------------------
## S14 -- the samples= sweep
##
## resolve_samples() unified ~20 call sites that each resolved samples= slightly
## differently, and deliberately preserved the quirks -- including four plotting
## functions whose existence check is an unanchored grep() substring match rather
## than %in%. On a 3-sample fixture a wrong flag still returns all three samples,
## so without an explicit sweep this is a permanent false clean.
## ---------------------------------------------------------------------------

xv_s14 <- function(){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)
  st <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = 1, verbose = FALSE)
  st <- transform_data(st, method = "log", cores = 1)

  probes <- list(
    numeric_one   = 1,
    numeric_two   = c(1, 3),
    exact_name    = "alpha",
    two_names     = c("alpha", "charlie"),
    ## 'a' is a substring of alpha and charlie but an exact match for neither:
    ## the difference between grep() and %in% semantics.
    substring     = "a",
    missing_name  = "not_a_sample",
    out_of_range  = 99
  )
  fns <- list(
    plot_counts   = function(st, s) plot_counts(st, samples = s, data_type = "tr"),
    plot_image    = function(st, s) plot_image(st, samples = s),
    STplot        = function(st, s) STplot(st, genes = "Gene1", samples = s, data_type = "tr"),
    distribution  = function(st, s) distribution_plots(st, genes = "Gene1", samples = s, data_type = "tr"),
    SThet         = function(st, s) SThet(st, genes = c("Gene1","GradGene1"), samples = s, cores = 1, verbose = FALSE),
    summarize     = function(st, s) summarize_STlist(st)
  )
  out <- list()
  for(fn in names(fns)){
    for(p in names(probes)){
      out[[paste0("s14_", fn, "_", p)]] <- xv_capture(fns[[fn]](st, probes[[p]]), seed = XV_SEED)
    }
  }
  out
}
