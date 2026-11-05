## Defines the single, deterministic call sequence exercised by both
## dev/baseline/capture_baseline.R (writes tests/golden/*.rds) and
## tests/testthat/test-golden.R (recomputes fresh and diffs against them).
## Kept in its own file so capture and comparison always run identical code.
##
## bl_pipeline(cores) threads one STlist through every exported function in
## R/, in dependency order, and returns a named list of every intermediate
## and terminal result. Each step is wrapped in tryCatch: a function whose
## current, correct behavior is to error on this fixture (e.g.
## pseudobulk_samples needs >=4 samples) has that condition captured as the
## golden result instead of aborting the whole run.
##
## bl_pipeline() does not call devtools::load_all() or source the fixtures --
## callers are responsible for having the package and
## inst/testdata/make_fixtures.R already available.

bl_capture <- function(name, expr, results, cores){
  set.seed(20260115)
  val <- tryCatch(
    list(ok = TRUE, value = expr),
    error = function(e) list(ok = FALSE, value = conditionMessage(e))
  )
  results[[name]] <- val
  results
}

bl_pipeline <- function(cores, verbose = FALSE){
  base <- fixture_base_data()
  fx <- fixture_write_list_dfs(base)

  results <- list()

  results <- bl_capture("STlist", {
    STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = cores, verbose = verbose)
  }, results, cores)
  st0 <- results[["STlist"]]$value

  results <- bl_capture("filter_data", {
    filter_data(st0, gene_minreads = 0, spot_minreads = 0)
  }, results, cores)
  st1 <- if(results[["filter_data"]]$ok) results[["filter_data"]]$value else st0

  results <- bl_capture("transform_data_log", {
    transform_data(st1, method = "log", cores = cores)
  }, results, cores)
  st2 <- if(results[["transform_data_log"]]$ok) results[["transform_data_log"]]$value else st1

  ## sctransform's internal NB regression produces NA parameter estimates for
  ## a zero-variance gene (all-zero or constant), which trips an unguarded
  ## if() downstream -- a sctransform limitation, not a spatialGE bug. Those
  ## two genes stay in every other capture; only this SCT-specific branch
  ## drops them so the golden result reflects SCT's real, working output.
  results <- bl_capture("transform_data_sct", {
    st1_sct <- filter_data(st1, rm_genes = c("ZeroGene", "ConstGene"))
    transform_data(st1_sct, method = "sct", sct_min_cells = 1, cores = cores)
  }, results, cores)

  het_genes <- c("Gene1", "Gene2", "GradGene1", "GradGene2", "DomainGeneA", "DomainGeneB", "ZeroGene", "ConstGene")
  results <- bl_capture("SThet", {
    SThet(st2, genes = het_genes, cores = cores, verbose = verbose)
  }, results, cores)
  st3 <- if(results[["SThet"]]$ok) results[["SThet"]]$value else st2

  results <- bl_capture("compare_SThet", {
    compare_SThet(st3, genes = het_genes[1:2])
  }, results, cores)

  results <- bl_capture("STclust", {
    STclust(st3, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100, cores = cores, verbose = verbose)
  }, results, cores)
  st4 <- if(results[["STclust"]]$ok) results[["STclust"]]$value else st3

  stclust_col <- "stclust_spw0.025_dsplFalse"
  cluster_labels <- unique(unlist(lapply(st4@spatial_meta, function(d){
    if(stclust_col %in% colnames(d)) d[[stclust_col]] else NULL
  })))
  ref_cluster <- if(length(cluster_labels) > 0) cluster_labels[1] else NULL

  results <- bl_capture("STdiff", {
    STdiff(st4, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 100, test_type = "t_test", cores = cores)
  }, results, cores)

  results <- bl_capture("STdiff_volcano", {
    if(results[["STdiff"]]$ok) STdiff_volcano(results[["STdiff"]]$value) else stop("STdiff did not succeed")
  }, results, cores)

  results <- bl_capture("STenrich", {
    STenrich(st4,
             gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                               domain = c("DomainGeneA", "DomainGeneB"),
                               degenerate = c("ZeroGene", "ConstGene")),
             annot = stclust_col, domain = ref_cluster,
             reps = 20, min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
  }, results, cores)

  results <- bl_capture("STgradient", {
    STgradient(st4, annot = stclust_col, ref = ref_cluster, topgenes = 100, cores = cores, verbose = verbose)
  }, results, cores)

  results <- bl_capture("STenrich_gsva", {
    STenrich(st4,
             gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                               domain = c("DomainGeneA", "DomainGeneB")),
             score_type = "gsva", min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
  }, results, cores)

  krige_genes <- c("GradGene1", "DomainGeneA")
  results <- bl_capture("gene_interpolation", {
    gene_interpolation(st4, genes = krige_genes, cores = cores, verbose = verbose)
  }, results, cores)
  st5 <- if(results[["gene_interpolation"]]$ok) results[["gene_interpolation"]]$value else st4

  results <- bl_capture("STplot_interpolation", {
    STplot_interpolation(st5, genes = krige_genes)
  }, results, cores)

  ## The pseudobulk family runs on its own 6-sample fixture, because
  ## pseudobulk_samples() requires >= 4 samples and st5 has 3. Threading the
  ## 3-sample STlist through here captured an error string for all three
  ## functions, so their real output -- including the gene ranking that
  ## pseudobulk_heatmap() displays -- had no baseline coverage.
  ##
  ## calc_umap stays FALSE: uwot::umap is unseeded and multithreaded (fidelity
  ## rule 1), so it is not reproducible and cannot be baselined.
  pb_base <- fixture_pbulk_base_data()
  pb_fx <- fixture_write_pbulk(pb_base, file.path(tempdir(), paste0("bl-pbulk-cores", cores)))

  results <- bl_capture("pseudobulk_STlist", {
    STlist(rnacounts = pb_fx$rnacounts, spotcoords = pb_fx$spotcoords,
           samples = pb_fx$samples, cores = cores, verbose = verbose)
  }, results, cores)
  pb_st <- results[["pseudobulk_STlist"]]$value

  ## ZeroGene and ConstGene stay in deliberately: both have an identical
  ## pseudobulk total in every sample, which used to kill prcomp(). They are the
  ## fixture's standing regression case for that, so this step must not filter
  ## them out.
  results <- bl_capture("pseudobulk_samples", {
    pb_tr <- transform_data(pb_st, method = "log", cores = cores)
    pseudobulk_samples(pb_tr, calc_umap = FALSE)
  }, results, cores)

  pb_ok <- results[["pseudobulk_samples"]]$ok
  results <- bl_capture("pseudobulk_dim_plot", {
    if(!pb_ok) stop("pseudobulk_samples did not succeed")
    pseudobulk_dim_plot(results[["pseudobulk_samples"]]$value, plot_meta = "group")
  }, results, cores)

  results <- bl_capture("pseudobulk_heatmap", {
    if(!pb_ok) stop("pseudobulk_samples did not succeed")
    pseudobulk_heatmap(results[["pseudobulk_samples"]]$value, hm_display_genes = 5)
  }, results, cores)

  results <- bl_capture("distribution_plots", {
    distribution_plots(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  }, results, cores)

  results <- bl_capture("plot_counts", {
    plot_counts(st5, data_type = "tr")
  }, results, cores)

  results <- bl_capture("plot_image", {
    plot_image(st5)
  }, results, cores)

  results <- bl_capture("STplot", {
    STplot(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  }, results, cores)

  results <- bl_capture("summarize_STlist", {
    summarize_STlist(st5)
  }, results, cores)

  results <- bl_capture("get_gene_meta", {
    get_gene_meta(st5, sthet_only = TRUE)
  }, results, cores)

  img_dir <- file.path(tempdir(), paste0("bl-images-cores", cores))
  dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
  for(s in base$samples) fixture_write_image(file.path(img_dir, paste0(s, ".png")))
  results <- bl_capture("load_images", {
    load_images(st5, images = img_dir)
  }, results, cores)

  results <- bl_capture("spatial_metadata", {
    spatial_metadata(st5)
  }, results, cores)

  results <- bl_capture("tissue_names", {
    tissue_names(st5)
  }, results, cores)

  results
}
