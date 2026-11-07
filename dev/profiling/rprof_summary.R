## Phase 9 profiling harness: Rprof()/summaryRprof() text tables per exported
## function, diffable across runs so regressions show up as a table diff
## rather than requiring re-reading an HTML profile. Threads the same
## fixture and call sequence as dev/baseline/pipeline.R (mirrored
## independently, not via bl_pipeline() -- see profile_all.R's header for
## why).
##
## Rprof(memory.profiling = TRUE) is on deliberately: full-matrix
## expandSparse is primarily an allocation problem, not a time problem, and
## won't show up in a CPU-only (time) profile.
##
## Run from the package root:
##   Rscript dev/profiling/rprof_summary.R

devtools::load_all(".", quiet = TRUE)
source("inst/testdata/make_fixtures.R")

options(spatialGE.cores = 1)
cores <- 1
verbose <- FALSE

out_dir <- "dev/profiling/out/rprof"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Profiles one step with Rprof()+summaryRprof(), writes a diffable text
## table, and returns list(ok, value) in the same shape as
## dev/baseline/pipeline.R's bl_capture(), so a step whose current, correct
## behavior is to error on this fixture (e.g. pseudobulk_samples needs >=4
## samples) is captured rather than aborting the whole run.
rprof_step <- function(name, thunk){
  ok <- TRUE
  err_msg <- NULL

  prof_file <- tempfile(fileext = ".rprof")
  Rprof(prof_file, interval = 0.001, memory.profiling = TRUE)
  value <- tryCatch(thunk(), error = function(e){
    ok <<- FALSE
    err_msg <<- conditionMessage(e)
    NULL
  })
  Rprof(NULL)

  summary <- tryCatch(summaryRprof(prof_file, memory = "both"), error = function(e) NULL)
  unlink(prof_file)

  lines <- c(
    paste0("step: ", name),
    paste0("ok: ", ok),
    if(!ok) paste0("error: ", err_msg) else NULL,
    ""
  )
  if(is.null(summary)){
    lines <- c(lines, "(no samples collected -- step ran faster than the sampling interval)")
  } else{
    lines <- c(lines,
              sprintf("sample.interval: %s", summary$sample.interval),
              sprintf("sampling.time: %s", summary$sampling.time),
              "",
              "by.self:",
              utils::capture.output(print(summary$by.self)),
              "",
              "by.total:",
              utils::capture.output(print(summary$by.total)))
  }
  writeLines(lines, file.path(out_dir, paste0(name, ".txt")))

  list(ok = ok, value = if(ok) value else err_msg)
}

## Forces rendering of any plot object so its cost is captured in the
## profile. ComplexHeatmap objects are already drawn by the exported
## function itself, so they fall through here untouched.
force_plots <- function(p){
  if(inherits(p, "ggplot")){
    ggplot2::ggplot_build(p)
  } else if(is.list(p)){
    lapply(p, force_plots)
  }
  invisible(NULL)
}

base <- fixture_base_data()
fx <- fixture_write_list_dfs(base)

results <- list()

results[["STlist"]] <- rprof_step("STlist", function(){
  STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = cores, verbose = verbose)
})
st0 <- results[["STlist"]]$value

results[["filter_data"]] <- rprof_step("filter_data", function(){
  filter_data(st0, gene_minreads = 0, spot_minreads = 0)
})
st1 <- if(results[["filter_data"]]$ok) results[["filter_data"]]$value else st0

results[["transform_data_log"]] <- rprof_step("transform_data_log", function(){
  transform_data(st1, method = "log", cores = cores)
})
st2 <- if(results[["transform_data_log"]]$ok) results[["transform_data_log"]]$value else st1

results[["transform_data_sct"]] <- rprof_step("transform_data_sct", function(){
  st1_sct <- filter_data(st1, rm_genes = c("ZeroGene", "ConstGene"))
  transform_data(st1_sct, method = "sct", sct_min_cells = 1, cores = cores)
})

het_genes <- c("Gene1", "Gene2", "GradGene1", "GradGene2", "DomainGeneA", "DomainGeneB", "ZeroGene", "ConstGene")
results[["SThet"]] <- rprof_step("SThet", function(){
  SThet(st2, genes = het_genes, cores = cores, verbose = verbose)
})
st3 <- if(results[["SThet"]]$ok) results[["SThet"]]$value else st2

results[["compare_SThet"]] <- rprof_step("compare_SThet", function(){
  compare_SThet(st3, genes = het_genes[1:2])
})

results[["STclust"]] <- rprof_step("STclust", function(){
  STclust(st3, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100, cores = cores, verbose = verbose)
})
st4 <- if(results[["STclust"]]$ok) results[["STclust"]]$value else st3

stclust_col <- "stclust_spw0.025_dsplFalse"
cluster_labels <- unique(unlist(lapply(st4@spatial_meta, function(d){
  if(stclust_col %in% colnames(d)) d[[stclust_col]] else NULL
})))
ref_cluster <- if(length(cluster_labels) > 0) cluster_labels[1] else NULL

results[["STdiff"]] <- rprof_step("STdiff", function(){
  STdiff(st4, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 100, test_type = "t_test", cores = cores)
})

results[["STdiff_volcano"]] <- rprof_step("STdiff_volcano", function(){
  if(!results[["STdiff"]]$ok) stop("STdiff did not succeed")
  STdiff_volcano(results[["STdiff"]]$value)
})

results[["STenrich"]] <- rprof_step("STenrich", function(){
  STenrich(st4,
           gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                             domain = c("DomainGeneA", "DomainGeneB"),
                             degenerate = c("ZeroGene", "ConstGene")),
           annot = stclust_col, domain = ref_cluster,
           reps = 20, min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
})

results[["STgradient"]] <- rprof_step("STgradient", function(){
  STgradient(st4, annot = stclust_col, ref = ref_cluster, topgenes = 100, cores = cores, verbose = verbose)
})

results[["STenrich_gsva"]] <- rprof_step("STenrich_gsva", function(){
  STenrich(st4,
           gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                             domain = c("DomainGeneA", "DomainGeneB")),
           score_type = "gsva", min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
})

krige_genes <- c("GradGene1", "DomainGeneA")
results[["gene_interpolation"]] <- rprof_step("gene_interpolation", function(){
  gene_interpolation(st4, genes = krige_genes, cores = cores, verbose = verbose)
})
st5 <- if(results[["gene_interpolation"]]$ok) results[["gene_interpolation"]]$value else st4

results[["STplot_interpolation"]] <- rprof_step("STplot_interpolation", function(){
  p <- STplot_interpolation(st5, genes = krige_genes)
  force_plots(p)
  p
})

results[["pseudobulk_samples"]] <- rprof_step("pseudobulk_samples", function(){
  pseudobulk_samples(st5, calc_umap = FALSE)
})
pb_ok <- results[["pseudobulk_samples"]]$ok

results[["pseudobulk_dim_plot"]] <- rprof_step("pseudobulk_dim_plot", function(){
  if(!pb_ok) stop("pseudobulk_samples did not succeed")
  p <- pseudobulk_dim_plot(results[["pseudobulk_samples"]]$value)
  force_plots(p)
  p
})

results[["pseudobulk_heatmap"]] <- rprof_step("pseudobulk_heatmap", function(){
  if(!pb_ok) stop("pseudobulk_samples did not succeed")
  ## ComplexHeatmap::draw() inside pseudobulk_heatmap() already forces
  ## rendering, so no separate force_plots() call is needed here.
  pseudobulk_heatmap(results[["pseudobulk_samples"]]$value, hm_display_genes = 5)
})

results[["distribution_plots"]] <- rprof_step("distribution_plots", function(){
  p <- distribution_plots(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  force_plots(p)
  p
})

results[["plot_counts"]] <- rprof_step("plot_counts", function(){
  p <- plot_counts(st5, data_type = "tr")
  force_plots(p)
  p
})

results[["plot_image"]] <- rprof_step("plot_image", function(){
  p <- plot_image(st5)
  force_plots(p)
  p
})

results[["STplot"]] <- rprof_step("STplot", function(){
  p <- STplot(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  force_plots(p)
  p
})

results[["summarize_STlist"]] <- rprof_step("summarize_STlist", function(){
  summarize_STlist(st5)
})

results[["get_gene_meta"]] <- rprof_step("get_gene_meta", function(){
  get_gene_meta(st5, sthet_only = TRUE)
})

img_dir <- file.path(tempdir(), "rprof-summary-images")
dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
for(s in base$samples) fixture_write_image(file.path(img_dir, paste0(s, ".png")))
results[["load_images"]] <- rprof_step("load_images", function(){
  load_images(st5, images = img_dir)
})

results[["spatial_metadata"]] <- rprof_step("spatial_metadata", function(){
  spatial_metadata(st5)
})

results[["tissue_names"]] <- rprof_step("tissue_names", function(){
  tissue_names(st5)
})

failed <- names(results)[!vapply(results, function(r) r$ok, logical(1))]
if(length(failed) > 0){
  message("[rprof_summary] steps that errored (may be expected, e.g. pseudobulk_samples on this fixture): ",
          paste(failed, collapse = ", "))
}
message(sprintf("[rprof_summary] wrote %d table(s) to %s", length(results), out_dir))
