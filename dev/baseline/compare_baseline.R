## Comparison mechanics for the golden baseline (see dev/baseline/pipeline.R
## and capture_baseline.R). Never use identical() on round-tripped objects --
## it compares attribute order and environment identity, which differ across
## otherwise-equivalent runs. waldo::compare(tolerance = NULL) gives an exact
## structural diff instead.
##
## SCOPE OF THIS GATE -- read before wiring it into CI.
##
## Comparison here is bit-exact (tolerance = NULL everywhere). The pipeline
## routes through prcomp (LAPACK), stats::loess, gstat kriging, spaMM::fitme
## and uwot::umap, all of which reach BLAS/LAPACK, so results differ in the
## low bits across BLAS implementations (macOS Accelerate vs Linux
## OpenBLAS/reference) and across package-version bumps in those
## dependencies. That makes this a SAME-MACHINE, SAME-TOOLCHAIN gate: it
## answers "did my edit change the numbers on this box", which is precisely
## the question a refactor needs answered, and it is why tests/golden/ is
## gitignored rather than committed.
##
## Do not enable it on a CI runner of a different platform than the one that
## ran capture_baseline.R -- it will fail on float noise and the only way to
## green it is loosening the tolerance, which destroys the drift detection
## that is the whole point. The platform-portable net is the committed
## tests/testthat/test-numeric-*.R tier, which asserts against literal values
## at explicit, stated tolerances.

## Recursively strips call/terms/srcref attributes (which capture the call
## site and change whenever code moves, even with no behavior change) and
## replaces environment/closure leaves with a fixed sentinel (closures over
## divergent environments are never bit-comparable and aren't the contract).
normalize_for_baseline <- function(x){
  if(is.environment(x) || is.function(x)){
    return("<<sentinel:env-or-closure>>")
  }
  attr(x, "call") <- NULL
  attr(x, "terms") <- NULL
  attr(x, "srcref") <- NULL
  if(is.list(x)){
    x[] <- lapply(x, normalize_for_baseline)
  }
  x
}

## dgCMatrix: drop0() both sides first (an already-sparse-matrix path can
## leave explicitly stored zeros that a freshly-built matrix doesn't), then
## compare slots directly plus class (dgC vs dge vs ngC is a real regression).
bl_compare_dgCMatrix <- function(a, b){
  a <- Matrix::drop0(a)
  b <- Matrix::drop0(b)
  diffs <- list(
    class = waldo::compare(class(a), class(b)),
    Dim = waldo::compare(a@Dim, b@Dim),
    Dimnames = waldo::compare(a@Dimnames, b@Dimnames),
    i = waldo::compare(a@i, b@i),
    p = waldo::compare(a@p, b@p),
    x = waldo::compare(a@x, b@x, tolerance = NULL)
  )
  Filter(function(d) length(d) > 0, diffs)
}

## spdep listw/nb: compare neighbours/weights/style; region.id and sym
## asserted separately since they're easy to silently drop. attr(,"call")
## changes whenever the call site moves, so it's excluded rather than diffed.
bl_compare_listw <- function(a, b){
  diffs <- list(
    neighbours = waldo::compare(unclass(a$neighbours), unclass(b$neighbours)),
    weights = waldo::compare(a$weights, b$weights, tolerance = NULL),
    style = waldo::compare(a$style, b$style),
    region.id = waldo::compare(attr(a$neighbours, "region.id"), attr(b$neighbours, "region.id")),
    sym = waldo::compare(attr(a$neighbours, "sym"), attr(b$neighbours, "sym"))
  )
  Filter(function(d) length(d) > 0, diffs)
}

bl_compare_sf <- function(a, b){
  diffs <- list(
    wkt = waldo::compare(sf::st_as_text(sf::st_geometry(a)), sf::st_as_text(sf::st_geometry(b))),
    crs = waldo::compare(sf::st_crs(a)$wkt, sf::st_crs(b)$wkt),
    data = waldo::compare(sf::st_drop_geometry(a), sf::st_drop_geometry(b), tolerance = NULL)
  )
  Filter(function(d) length(d) > 0, diffs)
}

bl_compare_sp <- function(a, b){
  diffs <- list(
    coords = waldo::compare(sp::coordinates(a), sp::coordinates(b), tolerance = NULL),
    data = waldo::compare(a@data, b@data, tolerance = NULL)
  )
  Filter(function(d) length(d) > 0, diffs)
}

## ggplot: never compare the object graph itself (environments, S7 internals
## make it noisy and non-substantive). Compare the actual rendered payload:
## the named list's keys, each layer's built data, labels, and scale classes.
## Some geoms (e.g. ggforce::geom_sina) draw fresh random jitter inside
## ggplot_build() itself, not when the plot object is constructed -- so two
## independently-built plots differ even when the upstream data is identical.
## Re-seeding immediately before each build makes the comparison
## deterministic without touching the pipeline's own RNG usage.
bl_compare_ggplot_signature <- function(p){
  if(inherits(p, "bl_ggplot_sig")) return(p)
  set.seed(20260115)
  built <- ggplot2::ggplot_build(p)
  structure(list(
    data = lapply(built$data, as.data.frame),
    labels = p$labels,
    scale_limits = lapply(p$scales$scales, function(s) s$get_limits()),
    scale_classes = lapply(p$scales$scales, function(s) class(s)),
    layer_geoms = vapply(p$layers, function(l) class(l$geom)[1], character(1))
  ), class = "bl_ggplot_sig")
}

bl_compare_ggplot <- function(a, b){
  waldo::compare(bl_compare_ggplot_signature(a), bl_compare_ggplot_signature(b), tolerance = NULL)
}

## ComplexHeatmap: pseudobulk_heatmap() ends in ComplexHeatmap::draw(), so it
## returns a HeatmapList and the displayed matrix is at @ht_list[[1]]@matrix,
## NOT @matrix. Without this branch the object falls through to raw waldo,
## where normalize_for_baseline() -- which recurses into lists but not S4
## slots -- leaves grobs and environments to be compared as-is. That
## comparison is pure noise, and it is part of why the inverted-CV gene
## ranking went unnoticed for two releases. Compare the matrix that is drawn,
## plus the dendrogram-determined row/column order, which is what the ranking
## actually controls.
bl_compare_heatmap_signature <- function(h){
  if(inherits(h, "bl_heatmap_sig")) return(h)
  ht <- h@ht_list[[1]]
  structure(list(
    matrix = as.matrix(ht@matrix),
    row_order = ComplexHeatmap::row_order(h),
    column_order = ComplexHeatmap::column_order(h),
    matrix_dimnames = dimnames(ht@matrix)
  ), class = "bl_heatmap_sig")
}

bl_compare_heatmap <- function(a, b){
  waldo::compare(bl_compare_heatmap_signature(a), bl_compare_heatmap_signature(b), tolerance = NULL)
}

## Reduces a captured result to exactly the form bl_compare() will compare,
## so capture_baseline.R stores that instead of the full object. ggplot and
## HeatmapList objects serialize enormously (a single STplot_interpolation
## result is 1.5 MB) while their comparable signature is ~2 KB; reducing at
## capture time takes the whole baseline from 13.9 MB to 3.1 MB and makes both
## capture and comparison faster.
##
## Idempotent, so it is safe to apply to either side: bl_compare() reduces
## both, which means a golden file captured before this existed still compares
## correctly against a freshly reduced result, and the cores1-vs-cores3 test
## can compare two already-reduced files.
##
## Tradeoff: a reduced golden file cannot answer a question the signature does
## not already capture. Widening bl_compare_ggplot_signature() therefore
## requires re-running capture_baseline.R rather than just re-comparing.
bl_reduce_value <- function(v){
  if(inherits(v, "ggplot")) return(bl_compare_ggplot_signature(v))
  if(is(v, "HeatmapList")) return(bl_compare_heatmap_signature(v))
  if(is.list(v) && !is.object(v) && length(v) > 0){
    if(all(vapply(v, function(e) inherits(e, "ggplot") || inherits(e, "bl_ggplot_sig"), logical(1)))){
      return(lapply(v, bl_compare_ggplot_signature))
    }
  }
  v
}

bl_reduce <- function(captured){
  if(!is.list(captured) || !("ok" %in% names(captured))) return(captured)
  if(!isTRUE(captured$ok)) return(captured)
  captured$value <- bl_reduce_value(captured$value)
  captured
}

## STlist: walk slots in a fixed order so a diff always names the slot, not
## just "objects differ". @counts / @tr_counts hold per-sample dgCMatrix
## lists; @spatial_meta / @sample_meta are data frames; @misc$sthet$listws
## (if present) holds spdep listw objects written by SThet().
bl_compare_stlist <- function(a, b){
  stopifnot(is(a, "STlist"), is(b, "STlist"))
  diffs <- list()

  cmp_matrix_list <- function(la, lb){
    nm <- union(names(la), names(lb))
    out <- list()
    for(n in nm){
      if(!(n %in% names(la)) || !(n %in% names(lb))){
        out[[n]] <- "present in only one side"
        next
      }
      d <- if(is(la[[n]], "dgCMatrix") || is(lb[[n]], "dgCMatrix")){
        bl_compare_dgCMatrix(as(la[[n]], "CsparseMatrix"), as(lb[[n]], "CsparseMatrix"))
      } else{
        Filter(length, list(waldo::compare(normalize_for_baseline(la[[n]]), normalize_for_baseline(lb[[n]]), tolerance = NULL)))
      }
      if(length(d) > 0) out[[n]] <- d
    }
    out
  }

  d_counts <- cmp_matrix_list(a@counts, b@counts)
  if(length(d_counts) > 0) diffs$counts <- d_counts

  d_tr <- cmp_matrix_list(a@tr_counts, b@tr_counts)
  if(length(d_tr) > 0) diffs$tr_counts <- d_tr

  d_spatial <- waldo::compare(normalize_for_baseline(a@spatial_meta), normalize_for_baseline(b@spatial_meta), tolerance = NULL)
  if(length(d_spatial) > 0) diffs$spatial_meta <- d_spatial

  d_gene <- waldo::compare(normalize_for_baseline(a@gene_meta), normalize_for_baseline(b@gene_meta), tolerance = NULL)
  if(length(d_gene) > 0) diffs$gene_meta <- d_gene

  d_sample <- waldo::compare(normalize_for_baseline(a@sample_meta), normalize_for_baseline(b@sample_meta), tolerance = NULL)
  if(length(d_sample) > 0) diffs$sample_meta <- d_sample

  a_listws <- a@misc$sthet$listws
  b_listws <- b@misc$sthet$listws
  if(!is.null(a_listws) || !is.null(b_listws)){
    nm <- union(names(a_listws), names(b_listws))
    d_listws <- list()
    for(n in nm){
      d <- bl_compare_listw(a_listws[[n]], b_listws[[n]])
      if(length(d) > 0) d_listws[[n]] <- d
    }
    if(length(d_listws) > 0) diffs$misc_sthet_listws <- d_listws
  }

  ## scaled_pbulk_mtx is the direct input to pseudobulk_heatmap()'s gene
  ## ranking and sp_images is the entire output of load_images(); neither was
  ## compared before, so a change in either was invisible to this gate. Both
  ## are small on the fixture (a genes x samples matrix, and 20x20x3 numeric
  ## arrays), so compare them outright rather than by signature.
  plain_misc <- c("scaled_pbulk_mtx", "pbulk_pca", "pbulk_pca_var", "pbulk_umap",
                  "sp_images", "platform", "transform", "image_scaling")
  for(key in plain_misc){
    ka <- a@misc[[key]]
    kb <- b@misc[[key]]
    if(is.null(ka) && is.null(kb)) next
    d_key <- waldo::compare(normalize_for_baseline(ka), normalize_for_baseline(kb), tolerance = NULL)
    if(length(d_key) > 0) diffs[[paste0("misc_", key)]] <- d_key
  }

  ## @misc$sthet holds the listw/nb objects (compared structurally above,
  ## since waldo on a listw diffs its call attribute) plus other entries.
  ## Compare those others rather than letting the listws stand in for the
  ## whole key.
  a_sthet_rest <- a@misc$sthet[setdiff(names(a@misc$sthet), "listws")]
  b_sthet_rest <- b@misc$sthet[setdiff(names(b@misc$sthet), "listws")]
  if(length(a_sthet_rest) > 0 || length(b_sthet_rest) > 0){
    d_rest <- waldo::compare(normalize_for_baseline(a_sthet_rest), normalize_for_baseline(b_sthet_rest), tolerance = NULL)
    if(length(d_rest) > 0) diffs$misc_sthet_other <- d_rest
  }

  ## Guard against a new @misc key being added without baseline coverage: the
  ## census in dev/ARCHITECTURE.md is the contract, and this gate silently
  ## ignoring a key is how scaled_pbulk_mtx went two releases uncompared.
  uncompared <- setdiff(union(names(a@misc), names(b@misc)),
                        c(plain_misc, "sthet", "gene_krige"))
  if(length(uncompared) > 0){
    diffs$misc_uncompared_keys <- sprintf(
      "@misc keys with no baseline coverage: %s -- add them to compare_baseline.R and to the census in dev/ARCHITECTURE.md",
      paste(uncompared, collapse = ", "))
  }

  d_krige <- waldo::compare(normalize_for_baseline(a@gene_krige), normalize_for_baseline(b@gene_krige), tolerance = NULL)
  if(length(d_krige) > 0) diffs$gene_krige <- d_krige

  diffs
}

## Top-level dispatcher used by tests/testthat/test-golden.R. Returns NULL on
## a match, or a non-empty diff object/string otherwise.
bl_compare <- function(golden, fresh){
  if(!identical(golden$ok, fresh$ok)){
    return(sprintf("ok mismatch: golden=%s fresh=%s (fresh msg: %s)",
                    golden$ok, fresh$ok, if(!fresh$ok) fresh$value else ""))
  }
  if(!golden$ok){
    ## Both sides are the *same kind* of error path (e.g. pseudobulk_samples'
    ## <4-samples guard) -- the condition class/message is the contract here,
    ## not incidental object structure.
    d <- waldo::compare(golden$value, fresh$value)
    return(if(length(d) > 0) d else NULL)
  }

  ## Reduce both sides so a golden file captured either before or after
  ## bl_reduce() existed compares identically, and so the cores1-vs-cores3
  ## test can hand two already-reduced files straight through.
  a <- bl_reduce_value(golden$value)
  b <- bl_reduce_value(fresh$value)
  d <- if(inherits(a, "bl_ggplot_sig") || inherits(b, "bl_ggplot_sig")){
    bl_compare_ggplot(a, b)
  } else if(inherits(a, "bl_heatmap_sig") || inherits(b, "bl_heatmap_sig")){
    bl_compare_heatmap(a, b)
  } else if(is.list(a) && length(a) > 0 && all(vapply(a, function(e) inherits(e, "bl_ggplot_sig"), logical(1)))){
    Filter(length, Map(bl_compare_ggplot, a, b))
  } else if(is(a, "STlist")){
    bl_compare_stlist(a, b)
  } else if(is(a, "dgCMatrix")){
    bl_compare_dgCMatrix(a, b)
  } else if(inherits(a, "sf")){
    bl_compare_sf(a, b)
  } else if(inherits(a, "Spatial")){
    bl_compare_sp(a, b)
  } else if(inherits(a, "listw")){
    bl_compare_listw(a, b)
  } else{
    waldo::compare(normalize_for_baseline(a), normalize_for_baseline(b), tolerance = NULL)
  }

  if(length(d) > 0) d else NULL
}
