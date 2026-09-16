## Comparison mechanics for the golden baseline (see dev/baseline/pipeline.R
## and capture_baseline.R). Never use identical() on round-tripped objects --
## it compares attribute order and environment identity, which differ across
## otherwise-equivalent runs. waldo::compare(tolerance = NULL) gives an exact
## structural diff instead.

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
  set.seed(20260115)
  built <- ggplot2::ggplot_build(p)
  list(
    data = lapply(built$data, as.data.frame),
    labels = p$labels,
    scale_limits = lapply(p$scales$scales, function(s) s$get_limits()),
    scale_classes = lapply(p$scales$scales, function(s) class(s)),
    layer_geoms = vapply(p$layers, function(l) class(l$geom)[1], character(1))
  )
}

bl_compare_ggplot <- function(a, b){
  waldo::compare(bl_compare_ggplot_signature(a), bl_compare_ggplot_signature(b), tolerance = NULL)
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

  a <- golden$value
  b <- fresh$value
  d <- if(is(a, "STlist")){
    bl_compare_stlist(a, b)
  } else if(is(a, "dgCMatrix")){
    bl_compare_dgCMatrix(a, b)
  } else if(inherits(a, "sf")){
    bl_compare_sf(a, b)
  } else if(inherits(a, "Spatial")){
    bl_compare_sp(a, b)
  } else if(inherits(a, "listw")){
    bl_compare_listw(a, b)
  } else if(inherits(a, "ggplot")){
    bl_compare_ggplot(a, b)
  } else if(is.list(a) && length(a) > 0 && all(vapply(a, function(e) inherits(e, "ggplot"), logical(1)))){
    Filter(length, Map(bl_compare_ggplot, a, b))
  } else{
    waldo::compare(normalize_for_baseline(a), normalize_for_baseline(b), tolerance = NULL)
  }

  if(length(d) > 0) d else NULL
}
