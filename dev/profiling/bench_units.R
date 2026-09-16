## Phase 9 profiling harness: bench::mark() timings for the individual named
## per-sample/per-gene/per-combo/per-cluster units of work extracted in
## Phase 7 (e.g. stenrich_perm_one_combo, interp_krige_one_gene,
## filter_stage_counts), isolated from the surrounding sge_lapply/
## parallelization layer so one unit's cost is measurable directly, without
## the multi-sample orchestration muddying the number.
##
## Mirrors the same fixture as profile_all.R/rprof_summary.R and reuses a
## short prefix of dev/baseline/pipeline.R's call sequence to reach a
## realistic intermediate STlist state for each target (kept independent of
## pipeline.R itself, for the same reason documented in profile_all.R's
## header). From that state, each target function is called directly with
## the one sample/gene/combo/cluster slice it would receive in production.
##
## Writes a diffable CSV (one row per unit) to dev/profiling/out/, consistent
## with rprof_summary.R's plain-text-over-HTML choice for regression tracking.
##
## Run from the package root:
##   Rscript dev/profiling/bench_units.R

devtools::load_all(".", quiet = TRUE)
source("inst/testdata/make_fixtures.R")

options(spatialGE.cores = 1)
cores <- 1
verbose <- FALSE

out_dir <- "dev/profiling/out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results <- list()

## Benchmarks one unit of work with bench::mark(), tolerating errors (a unit
## that cannot be exercised on this tiny fixture is skipped, not fatal).
bench_unit <- function(name, thunk, iterations = 20){
  res <- tryCatch(
    bench::mark(thunk(), iterations = iterations, check = FALSE, filter_gc = FALSE),
    error = function(e){
      message(sprintf("[bench_units] %s errored: %s", name, conditionMessage(e)))
      NULL
    }
  )
  if(is.null(res)) return(invisible(NULL))
  results[[name]] <<- data.frame(
    name = name,
    min_s = as.numeric(res$min),
    median_s = as.numeric(res$median),
    itr_per_sec = as.numeric(res$`itr/sec`),
    mem_alloc_bytes = as.numeric(res$mem_alloc),
    n_gc = as.numeric(res$n_gc),
    n_itr = as.numeric(res$n_itr)
  )
  invisible(NULL)
}

## ---- Build a realistic intermediate STlist state (short prefix of pipeline.R) ----

base <- fixture_base_data()
fx <- fixture_write_list_dfs(base)
sample1 <- base$samples[1]

st0 <- STlist(rnacounts = fx$rnacounts, spotcoords = fx$spotcoords, cores = cores, verbose = verbose)
st1 <- filter_data(st0, gene_minreads = 0, spot_minreads = 0)
st2 <- transform_data(st1, method = "log", cores = cores)

het_genes <- c("Gene1", "Gene2", "GradGene1", "GradGene2", "DomainGeneA", "DomainGeneB", "ZeroGene", "ConstGene")
st3 <- SThet(st2, genes = het_genes, cores = cores, verbose = verbose)
st3_vst <- calculate_vst(x = st3, samples = names(st3@counts), cores = cores)

st4 <- STclust(st3, ws = 0.025, ks = "dtc", deepSplit = FALSE, topgenes = 100, cores = cores, verbose = verbose)
stclust_col <- "stclust_spw0.025_dsplFalse"

krige_genes <- c("GradGene1", "DomainGeneA")

## ---- filter_data_stages.R ----

bench_unit("filter_stage_samples", function() filter_stage_samples(st0, base$samples[length(base$samples)]))
bench_unit("filter_stage_spots", function() filter_stage_spots(st1, sample1, colnames(st1@counts[[sample1]])[1:2]))
bench_unit("filter_stage_genes", function() filter_stage_genes(st1, sample1, c("ZeroGene", "ConstGene")))
bench_unit("filter_stage_genes_expr", function() filter_stage_genes_expr(st1, sample1, "^Zero"))
bench_unit("filter_stage_counts", function(){
  filter_stage_counts(st1, sample1, spot_minreads = 0, spot_maxreads = NULL, spot_mingenes = 0,
                      spot_maxgenes = NULL, spot_minpct = 0, spot_maxpct = NULL, gene_minreads = 0,
                      gene_maxreads = NULL, gene_minspots = 0, gene_maxspots = NULL, gene_minpct = 0,
                      gene_maxpct = NULL, spot_pct_expr = "^MT-")
})

## ---- sparse_utils.R ----

bench_unit("expandSparse", function() expandSparse(st0@counts[[sample1]]))
sparse_df <- tibble::rownames_to_column(expandSparse(st0@counts[[sample1]]), var = "gene")
bench_unit("makeSparse", function() makeSparse(sparse_df))

## ---- transform_data.R ----

bench_unit("log_transf_one_sample", function() log_transf_one_sample(st1@counts[[sample1]], sample1, 10000))
st1_sct <- filter_data(st1, rm_genes = c("ZeroGene", "ConstGene"))
bench_unit("sct_transf_one_sample", function() sct_transf_one_sample(st1_sct@counts[[sample1]], sample1, 3000, 1))
bench_unit("gene_meta_one_sample", function() gene_meta_one_sample(st2@tr_counts[[sample1]]))

## ---- SThet.R ----

sthet_genes <- intersect(het_genes, rownames(sample_tr_counts(st3, sample1)))
bench_unit("sthet_stat_one_sample", function(){
  sthet_stat_one_sample(sample_tr_counts(st3, sample1), sthet_genes,
                        st3@misc[["sthet"]][["listws"]][[sample1]], method = "moran")
})

## ---- STclust.R ----

bench_unit("stclust_topgenes_one_sample", function(){
  stclust_topgenes_one_sample(st3_vst@gene_meta[[sample1]], st3_vst@tr_counts[[sample1]], 100)
})
topgenes_mat <- stclust_topgenes_one_sample(st3_vst@gene_meta[[sample1]], st3_vst@tr_counts[[sample1]], 100)

bench_unit("stclust_one_sample", function(){
  stclust_one_sample(topgenes_mat, st3_vst@spatial_meta[[sample1]], "euclidean", 0.025, "hclust",
                     "dtc", FALSE, "ward.D2", 1L, sample1)
})

bench_unit("calculate_dist_matrices", function(){
  calculate_dist_matrices(expr_dat = topgenes_mat, coord_dat = st3_vst@spatial_meta[[sample1]], dist_metric = "euclidean")
})
scaled_dists <- calculate_dist_matrices(expr_dat = topgenes_mat, coord_dat = st3_vst@spatial_meta[[sample1]], dist_metric = "euclidean")

bench_unit("calculate_weighted_dist", function() calculate_weighted_dist(scaled_dists = scaled_dists, ws = 0.025))
weighted_dists <- calculate_weighted_dist(scaled_dists = scaled_dists, ws = 0.025)

bench_unit("get_hier_clusters_dtc", function(){
  get_hier_clusters_dtc(weighted_dists = weighted_dists, ws = 0.025, deepSplit = FALSE, linkage = "ward.D2")
})
bench_unit("get_hier_clusters_ks", function(){
  get_hier_clusters_ks(weighted_dists = weighted_dists, ws = 0.025, ks = c(2, 3), linkage = "ward.D2")
})

## ---- STdiff_fit.R ----

st4_vst <- stdiff_ensure_vst(st4)
combo_prep <- tryCatch(
  stdiff_build_combo_table(st4_vst, names(st4_vst@counts), stclust_col, 100L, NULL, FALSE, 0L),
  error = function(e){
    message("[bench_units] stdiff combo table setup failed, skipping stdiff_nonspatial_one_sample: ", conditionMessage(e))
    NULL
  }
)
if(!is.null(combo_prep)){
  combo_df <- combo_prep$combo_df
  meta_dict <- combo_prep$meta_dict
  combo_clust1 <- combo_df[combo_df$samplename == sample1, ]

  bench_unit("stdiff_nonspatial_one_sample", function(){
    stdiff_nonspatial_one_sample(sample_tr_counts(st4_vst, sample1), sample_spatial_meta(st4_vst, sample1),
                                 combo_clust1, sample1, stclust_col, meta_dict, FALSE, 0L, "t_test")
  }, iterations = 5)
}

## stdiff_spatial_one_cluster is only reached when test_type='mm' (STdiff.R's
## own gate: `if(test_type == 'mm' & sp_topgenes > 0)`), so building faithful
## inputs means re-running the non-spatial 'mm' step and the spatial-test
## setup inline (mirrors stdiff_run_spatial_tests()'s body). Kept fully
## isolated in its own tryCatch: spaMM convergence on this 3-sample fixture
## is not guaranteed, and a failure here should not block the other units.
stdiff_spatial_setup <- tryCatch({
  combo_prep_mm <- stdiff_build_combo_table(st4_vst, sample1, stclust_col, 100L, NULL, FALSE, 0L)
  combo_df_mm <- combo_prep_mm$combo_df
  meta_dict_mm <- combo_prep_mm$meta_dict
  combo_clust_mm <- combo_df_mm[combo_df_mm$samplename == sample1, ]

  non_sp_one <- stdiff_nonspatial_one_sample(sample_tr_counts(st4_vst, sample1), sample_spatial_meta(st4_vst, sample1),
                                             combo_clust_mm, sample1, stclust_col, meta_dict_mm, FALSE, 0L, "mm")
  result_de_one <- stdiff_assemble_nonspatial_results(list(sample1 = non_sp_one), meta_dict_mm, FALSE, "mm", "fdr")
  names(result_de_one) <- sample1

  models_keep <- result_de_one[[sample1]] %>%
    dplyr::filter(adj_p_val < 0.99) %>%
    dplyr::arrange(adj_p_val) %>%
    dplyr::group_by(cluster_1) %>%
    dplyr::slice_head(prop = 1) %>%
    dplyr::ungroup()

  ## Recode cluster_1 to the coded annotation BEFORE building genexcluster --
  ## matches stdiff_run_spatial_tests()'s own ordering, since non_sp_models's
  ## list keys (built in non_spatial_de()) are keyed on the coded meta1, not
  ## the original annotation label.
  for(spotrow in seq_len(nrow(models_keep))){
    models_keep[['cluster_1']][spotrow] <- meta_dict_mm[['coded_annot']][ meta_dict_mm[['orig_annot']] == models_keep[['cluster_1']][spotrow] ]
  }

  models_keep <- models_keep %>%
    dplyr::select(c('sample', 'gene', 'cluster_1')) %>%
    tibble::add_column(genexcluster = paste0(.[['sample']], '_', .[['gene']], '_', .[['cluster_1']], '_other'))

  expr_df_mm <- expandSparse(st4_vst@tr_counts[[sample1]])
  expr_df_mm <- t(expr_df_mm) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = 'libname') %>%
    dplyr::right_join(st4_vst@spatial_meta[[sample1]] %>%
                        tibble::add_column(group = 1, .after = 'libname') %>%
                        dplyr::select(c('libname', 'group', 'ypos', 'xpos'), orig_annot := !!stclust_col), ., by = 'libname') %>%
    dplyr::left_join(., meta_dict_mm, by = 'orig_annot') %>%
    dplyr::rename(meta = coded_annot) %>%
    tibble::column_to_rownames(var = 'libname') %>%
    dplyr::relocate(meta, .before = 1) %>%
    dplyr::select(-c('orig_annot'))

  list(models_keep = models_keep, cluster_name = unique(models_keep[['cluster_1']])[1],
      non_sp_models_sample = non_sp_one, expr_df = expr_df_mm, meta_dict = meta_dict_mm)
}, error = function(e){
  message("[bench_units] stdiff spatial-fit setup failed, skipping stdiff_spatial_one_cluster: ", conditionMessage(e))
  NULL
})

if(!is.null(stdiff_spatial_setup) && length(stdiff_spatial_setup$cluster_name) == 1 && !is.na(stdiff_spatial_setup$cluster_name)){
  bench_unit("stdiff_spatial_one_cluster", function(){
    stdiff_spatial_one_cluster(stdiff_spatial_setup$models_keep, stdiff_spatial_setup$cluster_name,
                               stdiff_spatial_setup$non_sp_models_sample, stdiff_spatial_setup$expr_df,
                               sample1, stdiff_spatial_setup$meta_dict, FALSE, 0L)
  }, iterations = 1)
}

## ---- stenrich_helpers.R ----

gradient_genes <- intersect(c("GradGene1", "GradGene2", "GradGene3"), rownames(st4@tr_counts[[sample1]]))
delayed_mat1 <- DelayedArray::DelayedArray(st4@tr_counts[[sample1]])

bench_unit("stenrich_score_one_combo", function(){
  stenrich_score_one_combo(delayed_mat1, gradient_genes, "gradient", 1L)
})
gs_score1 <- stenrich_score_one_combo(delayed_mat1, gradient_genes, "gradient", 1L)
coords_mat1 <- as.matrix(tibble::column_to_rownames(st4@spatial_meta[[sample1]][, c('libname', 'xpos', 'ypos')], var = 'libname'))

bench_unit("stenrich_perm_one_combo", function(){
  stenrich_perm_one_combo(unlist(gs_score1["gradient", ]), coords_mat1, sample1, "gradient",
                          length(gradient_genes), 3L, 1, 1L, 20L, 12345, 0L)
})

## ---- STgradient.R ----

## The stclust column is a per-sample factor with its own level set (samples
## can settle on different numbers of clusters), so the ref cluster must come
## from sample1's own levels directly -- unlist()-ing factors across samples
## silently reinterprets integer codes against a merged level set and hands
## stgradient_dists_one_sample a ref that isn't actually one of sample1's
## levels, tripping `Ops.factor: level sets of factors are different`.
sample1_clusters <- as.character(st4@spatial_meta[[sample1]][[stclust_col]])
ref_cluster <- if(length(sample1_clusters) > 0) sample1_clusters[1] else NULL

stgradient_setup <- if(!is.null(ref_cluster)) tryCatch({
  dists_summ <- stgradient_dists_one_sample(st4@spatial_meta[[sample1]], stclust_col, ref_cluster, NULL,
                                            "min", c(0.75, 1.25), 3, NULL)
  expr_res <- stgradient_expr_one_sample(st4@counts[[sample1]], st4@tr_counts[[sample1]], st4@spatial_meta[[sample1]],
                                         dists_summ, 100)
  list(dists_summ = dists_summ, expr_res = expr_res)
}, error = function(e){
  message("[bench_units] STgradient setup failed: ", conditionMessage(e))
  NULL
}) else NULL

if(!is.null(stgradient_setup)){
  bench_unit("stgradient_dists_one_sample", function(){
    stgradient_dists_one_sample(st4@spatial_meta[[sample1]], stclust_col, ref_cluster, NULL, "min", c(0.75, 1.25), 3, NULL)
  })
  bench_unit("stgradient_expr_one_sample", function(){
    stgradient_expr_one_sample(st4@counts[[sample1]], st4@tr_counts[[sample1]], st4@spatial_meta[[sample1]],
                               stgradient_setup$dists_summ, 100)
  })
  bench_unit("stgradient_one_sample", function(){
    stgradient_one_sample(st4@counts[[sample1]], st4@tr_counts[[sample1]], st4@spatial_meta[[sample1]], sample1,
                          stclust_col, ref_cluster, NULL, "min", c(0.75, 1.25), 3, NULL, 100, FALSE, TRUE, FALSE)
  })

  vargenes_expr <- stgradient_setup$expr_res[["vargenes_expr"]]
  gene_cols <- setdiff(colnames(vargenes_expr), c('ypos', 'xpos', 'dist2ref'))
  if(length(gene_cols) > 0){
    fit_gene <- gene_cols[1]
    df_gene <- vargenes_expr %>% dplyr::select(dist2ref, !!fit_gene)
    bench_unit("stgradient_fit_one_gene", function(){
      stgradient_fit_one_gene(df_gene, fit_gene, FALSE, TRUE, FALSE, NULL)
    })
  }
}

## ---- gene_interpolation_helpers.R ----

bench_unit("interp_krige_one_gene", function(){
  interp_krige_one_gene(st4@tr_counts[[sample1]], st4@spatial_meta[[sample1]], "GradGene1", sample1, 0L)
}, iterations = 5)
bench_unit("interp_one_sample", function(){
  interp_one_sample(st4@tr_counts[[sample1]], st4@spatial_meta[[sample1]], krige_genes, sample1, 0L)
}, iterations = 5)

## ---- Write results ----

out_df <- do.call(rbind, results)
out_path <- file.path(out_dir, "bench_units.csv")
if(is.null(out_df)){
  message("[bench_units] no units produced a result")
} else{
  utils::write.csv(out_df[order(out_df$name), ], out_path, row.names = FALSE)
  message(sprintf("[bench_units] wrote %d unit(s) to %s", nrow(out_df), out_path))
}
