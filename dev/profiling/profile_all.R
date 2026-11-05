## Phase 9 profiling harness: runs profvis::profvis() around every exported
## function, threaded through the same fixture and call sequence as
## dev/baseline/pipeline.R. Does not call bl_pipeline() directly or modify
## pipeline.R -- that file's header comment states capture and comparison
## must always run identical code, and wrapping each of its steps in
## profvis() would inject profiling instrumentation into that shared path.
## This script instead mirrors pipeline.R's sequence and arguments
## independently.
##
## Run from the package root:
##   Rscript dev/profiling/profile_all.R

devtools::load_all(".", quiet = TRUE)
source("inst/testdata/make_fixtures.R")

## profvis can't see into forked children (parallel::mclapply), so at
## cores > 1 it would mostly measure fork() rather than the actual work.
options(spatialGE.cores = 1)
cores <- 1
verbose <- FALSE

out_dir <- "dev/profiling/out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Forces rendering of any plot object so its cost is captured in the
## profile rather than deferred to whenever a caller happens to print it.
## ComplexHeatmap objects (pseudobulk_heatmap's return value) are already
## drawn by the exported function itself via ComplexHeatmap::draw(), so
## they fall through here untouched.
force_plots <- function(p){
  if(inherits(p, "ggplot")){
    ggplot2::ggplot_build(p)
  } else if(is.list(p)){
    lapply(p, force_plots)
  }
  invisible(NULL)
}

## Profiles one step: runs `thunk` inside profvis::profvis(), saves the
## profile HTML, and returns list(ok, value) in the same shape as
## dev/baseline/pipeline.R's bl_capture(), so a step whose current, correct
## behavior is to error on this fixture (e.g. pseudobulk_samples needs >=4
## samples) is captured rather than aborting the whole run. Targets its own
## frame via an explicit environment (not `<<-`) so the result doesn't
## depend on which environment profvis happens to evaluate the block in.
##
## profvis itself throws ("No parsing data available") when `thunk` runs
## faster than the sampling interval -- on this small synthetic fixture,
## several steps do. That failure happens after `thunk` has already run, so
## `computed` distinguishes it from a real failure of `thunk`: if `thunk`
## already produced a value, a missing profile is a no-op (skip the HTML),
## not a reason to discard that value.
profile_step <- function(name, thunk){
  self_env <- environment()
  ok <- TRUE
  err_msg <- NULL
  value <- NULL
  computed <- FALSE

  prof <- tryCatch(
    profvis::profvis({
      assign("value", tryCatch(thunk(), error = function(e){
        assign("ok", FALSE, envir = self_env)
        assign("err_msg", conditionMessage(e), envir = self_env)
        NULL
      }), envir = self_env)
      assign("computed", TRUE, envir = self_env)
    }, interval = 0.001),
    error = function(e){
      if(!computed){
        assign("ok", FALSE, envir = self_env)
        assign("err_msg", paste("profiling failed:", conditionMessage(e)), envir = self_env)
      } else{
        message(sprintf("[profile_all] %s: no profile captured (ran faster than the sampling interval)", name))
      }
      NULL
    }
  )

  if(!is.null(prof)){
    htmlwidgets::saveWidget(prof, file.path(out_dir, paste0(name, ".html")), selfcontained = TRUE)
  }

  if(!ok){
    message(sprintf("[profile_all] %s errored: %s", name, err_msg))
  }

  list(ok = ok, value = if(ok) value else err_msg)
}

base <- fixture_base_data()
fx <- fixture_write_list_dfs(base)

results <- list()

results[["STlist"]] <- profile_step("STlist", function(){
  STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = cores, verbose = verbose)
})
st0 <- results[["STlist"]]$value

results[["filter_data"]] <- profile_step("filter_data", function(){
  filter_data(st0, gene_minreads = 0, spot_minreads = 0)
})
st1 <- if(results[["filter_data"]]$ok) results[["filter_data"]]$value else st0

results[["transform_data_log"]] <- profile_step("transform_data_log", function(){
  transform_data(st1, method = "log", cores = cores)
})
st2 <- if(results[["transform_data_log"]]$ok) results[["transform_data_log"]]$value else st1

results[["transform_data_sct"]] <- profile_step("transform_data_sct", function(){
  st1_sct <- filter_data(st1, rm_genes = c("ZeroGene", "ConstGene"))
  transform_data(st1_sct, method = "sct", sct_min_cells = 1, cores = cores)
})

het_genes <- c("Gene1", "Gene2", "GradGene1", "GradGene2", "DomainGeneA", "DomainGeneB", "ZeroGene", "ConstGene")
results[["SThet"]] <- profile_step("SThet", function(){
  SThet(st2, genes = het_genes, cores = cores, verbose = verbose)
})
st3 <- if(results[["SThet"]]$ok) results[["SThet"]]$value else st2

results[["compare_SThet"]] <- profile_step("compare_SThet", function(){
  compare_SThet(st3, genes = het_genes[1:2])
})

results[["STclust"]] <- profile_step("STclust", function(){
  STclust(st3, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100, cores = cores, verbose = verbose)
})
st4 <- if(results[["STclust"]]$ok) results[["STclust"]]$value else st3

stclust_col <- "stclust_spw0.025_dsplFalse"
cluster_labels <- unique(unlist(lapply(st4@spatial_meta, function(d){
  if(stclust_col %in% colnames(d)) d[[stclust_col]] else NULL
})))
ref_cluster <- if(length(cluster_labels) > 0) cluster_labels[1] else NULL

results[["STdiff"]] <- profile_step("STdiff", function(){
  STdiff(st4, w = 0.025, k = "dtc", deepSplit = FALSE, topgenes = 100, test_type = "t_test", cores = cores)
})

results[["STdiff_volcano"]] <- profile_step("STdiff_volcano", function(){
  if(!results[["STdiff"]]$ok) stop("STdiff did not succeed")
  STdiff_volcano(results[["STdiff"]]$value)
})

results[["STenrich"]] <- profile_step("STenrich", function(){
  STenrich(st4,
           gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                             domain = c("DomainGeneA", "DomainGeneB"),
                             degenerate = c("ZeroGene", "ConstGene")),
           annot = stclust_col, domain = ref_cluster,
           reps = 20, min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
})

results[["STgradient"]] <- profile_step("STgradient", function(){
  STgradient(st4, annot = stclust_col, ref = ref_cluster, topgenes = 100, cores = cores, verbose = verbose)
})

results[["STenrich_gsva"]] <- profile_step("STenrich_gsva", function(){
  STenrich(st4,
           gene_sets = list(gradient = c("GradGene1", "GradGene2", "GradGene3"),
                             domain = c("DomainGeneA", "DomainGeneB")),
           score_type = "gsva", min_units = 1, min_genes = 1, cores = cores, verbose = verbose)
})

krige_genes <- c("GradGene1", "DomainGeneA")
results[["gene_interpolation"]] <- profile_step("gene_interpolation", function(){
  gene_interpolation(st4, genes = krige_genes, cores = cores, verbose = verbose)
})
st5 <- if(results[["gene_interpolation"]]$ok) results[["gene_interpolation"]]$value else st4

results[["STplot_interpolation"]] <- profile_step("STplot_interpolation", function(){
  p <- STplot_interpolation(st5, genes = krige_genes)
  force_plots(p)
  p
})

results[["pseudobulk_samples"]] <- profile_step("pseudobulk_samples", function(){
  pseudobulk_samples(st5, calc_umap = FALSE)
})
pb_ok <- results[["pseudobulk_samples"]]$ok

results[["pseudobulk_dim_plot"]] <- profile_step("pseudobulk_dim_plot", function(){
  if(!pb_ok) stop("pseudobulk_samples did not succeed")
  p <- pseudobulk_dim_plot(results[["pseudobulk_samples"]]$value)
  force_plots(p)
  p
})

results[["pseudobulk_heatmap"]] <- profile_step("pseudobulk_heatmap", function(){
  if(!pb_ok) stop("pseudobulk_samples did not succeed")
  ## ComplexHeatmap::draw() inside pseudobulk_heatmap() already forces
  ## rendering, so no separate force_plots() call is needed here.
  pseudobulk_heatmap(results[["pseudobulk_samples"]]$value, hm_display_genes = 5)
})

results[["distribution_plots"]] <- profile_step("distribution_plots", function(){
  p <- distribution_plots(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  force_plots(p)
  p
})

results[["plot_counts"]] <- profile_step("plot_counts", function(){
  p <- plot_counts(st5, data_type = "tr")
  force_plots(p)
  p
})

results[["plot_image"]] <- profile_step("plot_image", function(){
  p <- plot_image(st5)
  force_plots(p)
  p
})

results[["STplot"]] <- profile_step("STplot", function(){
  p <- STplot(st5, genes = c("Gene1", "GradGene1"), data_type = "tr")
  force_plots(p)
  p
})

results[["summarize_STlist"]] <- profile_step("summarize_STlist", function(){
  summarize_STlist(st5)
})

results[["get_gene_meta"]] <- profile_step("get_gene_meta", function(){
  get_gene_meta(st5, sthet_only = TRUE)
})

img_dir <- file.path(tempdir(), "profile-all-images")
dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)
for(s in base$samples) fixture_write_image(file.path(img_dir, paste0(s, ".png")))
results[["load_images"]] <- profile_step("load_images", function(){
  load_images(st5, images = img_dir)
})

results[["spatial_metadata"]] <- profile_step("spatial_metadata", function(){
  spatial_metadata(st5)
})

results[["tissue_names"]] <- profile_step("tissue_names", function(){
  tissue_names(st5)
})

failed <- names(results)[!vapply(results, function(r) r$ok, logical(1))]
if(length(failed) > 0){
  message("[profile_all] steps that errored (may be expected, e.g. pseudobulk_samples on this fixture): ",
          paste(failed, collapse = ", "))
}
message(sprintf("[profile_all] wrote %d profile(s) to %s", length(results), out_dir))
