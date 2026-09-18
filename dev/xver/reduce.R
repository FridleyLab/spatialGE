## Version-neutral reduction of spatialGE results to plain, comparable payloads.
##
## This file must not name a single spatialGE symbol. It runs inside a producer
## process that has one specific version loaded, and its whole job is to hand the
## comparator something that carries no reference back to that version: no S4
## objects, no environments, no closures. The comparator then runs with spatialGE
## absent from its library path, so any leak fails loudly instead of silently
## deserializing against the wrong class table.
##
## Two modes:
##   xv_reduce(x, canon = FALSE)  as produced
##   xv_reduce(x, canon = TRUE)   order-insensitive
##
## Canonicalization happens here, on the live objects, rather than later on the
## reduced payload -- reordering a dgCMatrix's columns is a one-liner on the
## matrix and a re-derivation of i/p/x on the payload.
##
## Rule the whole design rests on: an unrecognized class is an ERROR, never a
## silent fall-through. A reducer that quietly returns NULL manufactures clean
## comparisons, because NULL == NULL.

## ---- identity keys --------------------------------------------------------

## Sorting a result table must never use value columns: sorting by all columns
## would reorder rows to make values line up, hiding exactly the value diffs this
## harness exists to find. These are the columns that identify a row across the
## codebase's result tables, in priority order.
XV_ID_COLS <- c("sample_name", "sample", "ST_sample", "libname", "gene",
                "gene_set", "cluster_1", "cluster_2", "cluster", "spatialmod")

xv_identity_cols <- function(df){
  intersect(XV_ID_COLS, names(df))
}

## Rounded sort key. A single 1e-16 difference in a double flips an ordering and
## turns a value diff into an apparent whole-table permutation, so the key is
## rounded while the compared values stay exact.
xv_sort_key <- function(v){
  if(is.numeric(v)) formatC(signif(v, 12), format = "e", digits = 12) else as.character(v)
}

xv_order_by <- function(df, cols){
  if(!length(cols) || !nrow(df)) return(seq_len(nrow(df)))
  keys <- lapply(cols, function(cn) xv_sort_key(df[[cn]]))
  do.call(order, c(keys, list(method = "radix")))
}

## ---- leaf helpers ---------------------------------------------------------

xv_plain_vec <- function(v){
  ## Factors compare as their labels; a level-order change that does not change
  ## any label is not a behavior change worth reporting.
  if(is.factor(v)) return(as.character(v))
  if(inherits(v, c("Date", "POSIXct"))) return(as.character(v))
  if(is.list(v)) return(vapply(v, function(e) paste(format(e), collapse = "|"), character(1)))
  if(is.atomic(v)) return(v)
  stop("xv_plain_vec: unsupported column type: ", paste(class(v), collapse = "/"))
}

## Paths differ between the two runs for reasons that have nothing to do with the
## package: separate scratch roots, separate libraries, separate PIDs.
xv_scrub <- function(x){
  if(!is.character(x) || !length(x)) return(x)
  subs <- list(
    c(Sys.getenv("XVER_ROOT"), "<XVER_ROOT>"),
    c(Sys.getenv("XVER_LIB_MAIN"), "<LIB>"),
    c(Sys.getenv("XVER_LIB_ALEX"), "<LIB>"),
    c(Sys.getenv("XVER_SRC_MAIN"), "<SRC>"),
    c(Sys.getenv("XVER_SRC_ALEX"), "<SRC>"),
    c(tempdir(), "<TMP>")
  )
  for(s in subs){
    if(nzchar(s[1])) x <- gsub(s[1], s[2], x, fixed = TRUE)
  }
  x <- gsub("\\b[0-9]{3,7}\\b(?=[^0-9]*$)", "", x, perl = TRUE)
  x
}

## ---- data frames ----------------------------------------------------------

xv_reduce_df <- function(d, canon){
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  idc <- xv_identity_cols(d)
  ambiguous <- NA
  if(canon && nrow(d)){
    if(length(idc)){
      key <- do.call(paste, c(lapply(idc, function(cn) xv_sort_key(d[[cn]])), list(sep = "\r")))
      ## A duplicated identity key means a value difference inside a duplicate
      ## cannot be attributed to a row. Report it rather than sorting anyway.
      ambiguous <- anyDuplicated(key) > 0
      ord <- xv_order_by(d, idc)
    } else {
      ## No declared key: fall back to a whole-row multiset comparison, keyed on
      ## rounded values. Flagged, because a value diff can hide here.
      ambiguous <- TRUE
      ord <- xv_order_by(d, names(d))
    }
    d <- d[ord, , drop = FALSE]
  }
  cols <- lapply(d, function(v) xv_scrub(xv_plain_vec(v)))
  names(cols) <- names(d)
  list(
    .type = "df",
    class = class(d),
    nrow = nrow(d),
    ## Column ORDER is API (fidelity rule 8), so it is compared as data and never
    ## canonicalized away.
    colnames = names(d),
    ## main's row-subsetting leaves non-contiguous row names ("1","5","9"), which
    ## is the actual mechanism of the STenrich domain= bug. They are attributes,
    ## not columns, so sorting rows would not normalize them: keep as their own
    ## field and let the classifier decide.
    row_names = if(canon) NULL else xv_scrub(as.character(attr(d, "row.names"))),
    identity_cols = idc,
    key_ambiguous = ambiguous,
    columns = cols
  )
}

## ---- matrices -------------------------------------------------------------

xv_is_sparse <- function(x) methods::is(x, "Matrix")

xv_reduce_matrix <- function(m, canon){
  sparse <- xv_is_sparse(m)
  cls <- class(m)
  if(sparse){
    m <- methods::as(m, "CsparseMatrix")
    ## An already-sparse input path can leave explicitly stored zeros that a
    ## freshly built matrix does not; drop0 makes the two comparable.
    m <- Matrix::drop0(m)
  }
  rn <- rownames(m); cn <- colnames(m)
  if(canon){
    if(!is.null(rn)) m <- m[order(rn, method = "radix"), , drop = FALSE]
    if(!is.null(cn)) m <- m[, order(cn, method = "radix"), drop = FALSE]
  }
  out <- list(
    .type = "matrix",
    class = cls,
    dim = dim(m),
    rownames = rownames(m),
    colnames = colnames(m)
  )
  if(sparse){
    out$i <- m@i
    out$p <- m@p
    out$x <- m@x
  } else {
    out$values <- as.vector(m)
  }
  out
}

## ---- spatial types --------------------------------------------------------

## `ids` supplies the spot names when the object cannot name them itself.
## spatialGE builds these via create_listw_from_dist() over a distance matrix in
## slot order, and the resulting nb carries NO region.id attribute -- so the
## neighbour lists are purely positional, and position is exactly what the
## ordering fix changes. Without names from outside there is nothing to canonicalize
## against, and every weight vector would differ for a reason that is not a bug.
xv_reduce_listw <- function(x, canon, ids = NULL){
  nb <- unclass(x$neighbours)
  rid <- attr(x$neighbours, "region.id")
  w <- x$weights
  if(is.null(rid) || length(rid) != length(nb)) rid <- ids
  ## Refuse to reorder against ids that cannot address the entries one-to-one:
  ## silently reordering by a bad key is how a canonicalizer invents a clean pass.
  canonicalizable <- !is.null(rid) && length(rid) == length(nb) &&
                     anyDuplicated(rid) == 0L
  if(canon && canonicalizable){
    ## A neighbour list is positional: nb[[i]] holds integer positions into
    ## region.id, and weights[[i]] is PARALLEL to nb[[i]]. So canonicalizing has
    ## to do three things together, and doing only the first two silently breaks
    ## the neighbour-to-weight correspondence:
    ##   1. reorder the outer list into region.id order
    ##   2. remap the inner positions into the new indexing and sort them
    ##   3. apply that same inner permutation to the weights
    ord <- order(rid, method = "radix")
    remap <- match(seq_along(rid), ord)   # old position -> new position
    nb_new <- vector("list", length(ord))
    w_new  <- if(is.null(w)) NULL else vector("list", length(ord))
    for(newi in seq_along(ord)){
      oldi <- ord[newi]
      v <- remap[nb[[oldi]]]
      inner <- order(v, method = "radix")
      nb_new[[newi]] <- v[inner]
      if(!is.null(w_new)) w_new[[newi]] <- w[[oldi]][inner]
    }
    nb <- nb_new; w <- w_new
    rid <- rid[ord]
  }
  ## Name the per-region entries so a difference reports which spot it is about
  ## rather than a bare position that means different things on the two sides.
  if(!is.null(rid) && length(rid) == length(nb)){
    names(nb) <- rid
    if(!is.null(w) && length(w) == length(rid)) names(w) <- rid
  }
  list(.type = "listw", class = class(x), style = x$style,
       region.id = rid, sym = attr(x$neighbours, "sym"),
       ## Recorded so the report can distinguish "compared by name" from
       ## "compared positionally because no usable key existed".
       canonicalizable = canonicalizable,
       n = length(nb),
       neighbours = nb, weights = w)
}

xv_reduce_nb <- function(x, canon){
  xv_reduce_listw(list(neighbours = x, weights = NULL, style = NA_character_), canon)
}

xv_reduce_sf <- function(x, canon){
  d <- sf::st_drop_geometry(x)
  wkt <- sf::st_as_text(sf::st_geometry(x))
  if(canon){
    ord <- order(wkt, method = "radix")
    wkt <- wkt[ord]; d <- d[ord, , drop = FALSE]
  }
  list(.type = "sf", class = class(x), wkt = wkt,
       crs = sf::st_crs(x)$wkt, data = xv_reduce_df(d, canon))
}

xv_reduce_sp <- function(x, canon){
  co <- sp::coordinates(x)
  d <- if(methods::.hasSlot(x, "data")) x@data else data.frame()
  if(canon && nrow(co)){
    key <- do.call(paste, c(lapply(seq_len(ncol(co)), function(j) xv_sort_key(co[, j])), list(sep = "\r")))
    ord <- order(key, method = "radix")
    co <- co[ord, , drop = FALSE]
    if(nrow(d)) d <- d[ord, , drop = FALSE]
  }
  list(.type = "sp", class = class(x),
       coords = as.vector(co), coord_dim = dim(co),
       data = xv_reduce_df(d, canon))
}

## ---- ggplot ---------------------------------------------------------------

## Reuses bl_compare_ggplot_signature() from dev/baseline/compare_baseline.R
## rather than reimplementing it (ARCHITECTURE.md: extend, don't duplicate), then
## adds facet/coord/panel information -- the existing signature covers built data,
## labels, scales and geoms, so a faceting or coord regression would read clean.
xv_reduce_ggplot <- function(p, canon){
  sig <- bl_compare_ggplot_signature(p)
  built <- lapply(sig$data, function(d){
    d <- as.data.frame(d)
    ## Built-data row names are the original row positions, so they encode input
    ## order and nothing else.
    rownames(d) <- NULL
    if(!canon || !nrow(d)) return(xv_reduce_df(d, canon = FALSE))

    ## Built plot data has no identity key, so it cannot be sorted into a
    ## canonical row order the way a result table can. Sorting it by its own
    ## columns actively backfires: as soon as ONE column genuinely differs (a
    ## jittered coordinate, say) the sort scrambles the rows, and every other
    ## column then compares against the wrong row and reports a large difference.
    ## That mistake turned a violin layer agreeing to 1.1e-16 into an apparent
    ## 96% discrepancy.
    ##
    ## So compare it as a multiset instead, in two complementary ways:
    ##   columns        -- each column sorted independently: does the layer draw
    ##                     the same set of values for each aesthetic?
    ##   row_signatures -- the sorted rounded whole-row strings: is the
    ##                     correspondence BETWEEN columns preserved?
    ## The second is what keeps the first from being fooled by two columns having
    ## their values swapped between rows.
    for(cn in intersect(c("group", "PANEL"), names(d))){
      ## Integer codes assigned from factor level order, i.e. from data order, so
      ## they renumber legitimately under a permutation.
      d[[cn]] <- match(d[[cn]], unique(d[[cn]]))
    }
    ## Signatures are rounded to 6 significant digits, deliberately COARSER than
    ## the 1e-12 tolerance the columns are compared at. At matching precision,
    ## float noise straddling the rounding cut would flag rows whose columns all
    ## agree -- reporting the rounding boundary rather than the data. Six digits
    ## still catches a genuine row-correspondence swap, which is all this is for.
    sigs <- do.call(paste, c(lapply(d, function(v){
      if(is.numeric(v)) formatC(signif(v, 6), format = "e", digits = 6) else as.character(v)
    }), list(sep = "\r")))
    cols <- lapply(d, function(v){
      v <- xv_scrub(xv_plain_vec(v))
      if(is.numeric(v)) sort(v, na.last = TRUE) else sort(as.character(v), na.last = TRUE)
    })
    list(.type = "df_multiset", class = class(d), nrow = nrow(d),
         colnames = names(d), columns = cols,
         row_signatures = sort(sigs))
  })
  facet <- tryCatch(class(p$facet), error = function(e) NA_character_)
  layout <- tryCatch({
    l <- ggplot2::ggplot_build(p)$layout$layout
    xv_reduce_df(as.data.frame(l), canon = FALSE)
  }, error = function(e) NULL)
  list(
    .type = "ggplot",
    class = class(p),
    built = built,
    labels = lapply(sig$labels, function(v) xv_scrub(as.character(v))),
    scale_limits = lapply(sig$scale_limits, function(v) xv_plain_vec(v)),
    scale_classes = sig$scale_classes,
    layer_geoms = sig$layer_geoms,
    facet_class = facet,
    coord_class = tryCatch(class(p$coordinates), error = function(e) NA_character_),
    panel_layout = layout
  )
}

## ---- ComplexHeatmap -------------------------------------------------------

## pseudobulk_heatmap() ends in ComplexHeatmap::draw(), so it returns a
## HeatmapList and the matrix lives at @ht_list[[i]]@matrix, not @matrix. Left to
## the generic path this object is pure noise: bl_compare() falls through to
## waldo, and normalize_for_baseline() recurses only into lists, not S4 slots, so
## grobs and environments get compared raw. The selected gene set is the whole
## observable that the inverted-CV fix changes.
xv_reduce_one_heatmap <- function(ht, nm, canon){
  ## An annotation or dendrogram component carries no matrix; recording its class
  ## is enough, and pretending otherwise would error on a legitimate object.
  if(!methods::.hasSlot(ht, "matrix")){
    return(list(.type = "heatmap_component", class = class(ht), name = nm))
  }
  m <- ht@matrix
  list(.type = "heatmap", class = class(ht), name = nm,
       ## Gene selection is the observable the inverted-CV fix changes, and it is
       ## a SET question, so record it both sorted and in display order.
       genes_display_order = rownames(m),
       genes_sorted = sort(rownames(m)),
       matrix = xv_reduce_matrix(m, canon),
       row_title = tryCatch(as.character(ht@row_title), error = function(e) NA_character_),
       column_title = tryCatch(as.character(ht@column_title), error = function(e) NA_character_))
}

xv_reduce_heatmaplist <- function(x, canon){
  hts <- x@ht_list
  nms <- names(hts)
  if(is.null(nms)) nms <- paste0("[", seq_along(hts), "]")
  out <- lapply(seq_along(hts), function(i) xv_reduce_one_heatmap(hts[[i]], nms[i], canon))
  names(out) <- nms
  list(.type = "heatmap_list", class = class(x),
       ht_names = nms,
       row_order = tryCatch(ComplexHeatmap::row_order(x), error = function(e) NULL),
       column_order = tryCatch(ComplexHeatmap::column_order(x), error = function(e) NULL),
       heatmaps = out)
}

## ---- STlist ---------------------------------------------------------------

XV_MISC_KEYS <- c("platform", "transform", "sp_images", "sthet", "gene_krige",
                  "pbulk_pca", "pbulk_pca_var", "pbulk_umap", "scaled_pbulk_mtx",
                  "image_scaling")

## @misc$sthet$listws holds one spdep listw per sample, built over the sample's
## spots in slot order and carrying no region.id. The only place the spot names
## are available is the STlist itself, so they are threaded in here rather than
## leaving the neighbour lists positional and therefore incomparable.
xv_reduce_misc <- function(x, canon){
  misc <- x@misc
  keys <- intersect(XV_MISC_KEYS, names(misc))
  out <- list()
  for(k in keys){
    if(identical(k, "sthet") && is.list(misc[[k]]) && !is.null(misc[[k]][["listws"]])){
      lws <- misc[[k]][["listws"]]
      red <- list()
      for(s in names(lws)){
        sm <- x@spatial_meta[[s]]
        ids <- if(!is.null(sm) && "libname" %in% names(sm)) as.character(sm[["libname"]]) else NULL
        red[[s]] <- xv_reduce_listw(lws[[s]], canon, ids = ids)
      }
      rest <- setdiff(names(misc[[k]]), "listws")
      out[[k]] <- list(.type = "sthet", listws = red,
                       other = if(length(rest)) xv_reduce_any(misc[[k]][rest], canon) else NULL)
    } else {
      out[[k]] <- xv_reduce_any(misc[[k]], canon)
    }
  }
  out
}

xv_reduce_stlist <- function(x, canon){
  samples <- names(x@counts)
  per_sample <- function(slot, fn){
    v <- methods::slot(x, slot)
    if(!length(v)) return(NULL)
    nms <- names(v)
    o <- lapply(nms, function(s) if(is.null(v[[s]])) NULL else fn(v[[s]]))
    names(o) <- nms
    o
  }
  list(
    .type = "STlist",
    ## Sample order is user-visible and is NOT expected to change (fidelity rule
    ## 7), so it is compared as data and never canonicalized.
    sample_order_counts = names(x@counts),
    sample_order_spatial_meta = names(x@spatial_meta),
    sample_order_tr_counts = names(x@tr_counts),
    counts = per_sample("counts", function(m) xv_reduce_matrix(m, canon)),
    tr_counts = per_sample("tr_counts", function(m) xv_reduce_matrix(m, canon)),
    spatial_meta = per_sample("spatial_meta", function(d) xv_reduce_df(d, canon)),
    gene_meta = per_sample("gene_meta", function(d) xv_reduce_df(d, canon)),
    sample_meta = if(length(x@sample_meta)) xv_reduce_df(x@sample_meta, canon) else NULL,
    gene_krige = xv_reduce_any(x@gene_krige, canon),
    ## bl_compare_stlist() walks only @misc$sthet$listws, leaving eight of the ten
    ## keys the package writes uncompared -- including scaled_pbulk_mtx, the direct
    ## input to the pseudobulk_heatmap CV fix, and sp_images, the entire output of
    ## load_images(). Walk all of them.
    misc = xv_reduce_misc(x, canon),
    misc_unexpected_keys = setdiff(names(x@misc), XV_MISC_KEYS)
  )
}

## ---- generic dispatch -----------------------------------------------------

xv_reduce_any <- function(x, canon){
  if(is.null(x)) return(NULL)

  ## Fail loudly on anything that cannot cross a process boundary meaningfully.
  if(is.function(x))    stop("xv_reduce: refusing to reduce a function")
  if(is.environment(x)) stop("xv_reduce: refusing to reduce an environment")

  if(methods::is(x, "STlist"))       return(xv_reduce_stlist(x, canon))
  if(methods::is(x, "HeatmapList"))  return(xv_reduce_heatmaplist(x, canon))
  ## An undrawn Heatmap, should any caller return one instead of draw()'s output.
  if(methods::is(x, "Heatmap"))      return(xv_reduce_one_heatmap(x, "ht", canon))
  if(inherits(x, "ggplot"))          return(xv_reduce_ggplot(x, canon))
  if(inherits(x, "sf"))              return(xv_reduce_sf(x, canon))
  if(methods::is(x, "Spatial"))      return(xv_reduce_sp(x, canon))
  if(inherits(x, "listw"))           return(xv_reduce_listw(x, canon))
  if(inherits(x, "nb"))              return(xv_reduce_nb(x, canon))
  if(xv_is_sparse(x) || is.matrix(x)) return(xv_reduce_matrix(x, canon))
  if(is.data.frame(x))               return(xv_reduce_df(x, canon))

  ## Image rasters (@misc$sp_images, written by load_images()). A 20x20x3 array of
  ## doubles compared elementwise turns any difference into a wall of pixel diffs,
  ## so reduce to shape plus a digest: enough to detect a real change, compact
  ## enough to read. This also means a failure of the fixture_write_image() seed
  ## shim shows up as one flagged field instead of thousands.
  if(is.array(x) && length(dim(x)) >= 3){
    return(list(.type = "raster", class = class(x), dim = dim(x),
                range = if(is.numeric(x)) range(x) else NA,
                mean = if(is.numeric(x)) mean(x) else NA,
                digest = paste(utils::capture.output(print(sum(as.numeric(x)))), collapse = "")))
  }

  ## gstat/krige objects carry formulas, calls and environments. Reduce to the
  ## parts that are actually a contract.
  if(inherits(x, c("gstat", "variogramModel", "variogram"))){
    return(list(.type = "gstat", class = class(x),
                fields = lapply(x[intersect(names(x), c("model","psill","range","kappa","nugget"))],
                                function(v) tryCatch(xv_plain_vec(v), error = function(e) NA))))
  }

  if(is.atomic(x) && is.null(dim(x))){
    return(list(.type = "atomic", class = class(x),
                names = names(x), value = xv_scrub(xv_plain_vec(x))))
  }

  if(is.list(x)){
    ## An empty list still has to reduce: several STlist slots are legitimately
    ## length 0 before the step that populates them (@gene_krige before
    ## gene_interpolation, @gene_meta before transform_data). Note that
    ## paste0("[", integer(0), "]") is "[]" of length 1, not character(0), so the
    ## empty case needs its own branch rather than falling through the recycling.
    o <- lapply(seq_along(x), function(i) xv_reduce_any(x[[i]], canon))
    names(o) <- if(!length(x)) character(0)
                else if(is.null(names(x))) paste0("[", seq_along(x), "]")
                else names(x)
    return(list(.type = "list", class = class(x), length = length(x), elements = o))
  }

  ## An S4 object we have no extractor for. Erroring is the point: a generic
  ## fall-through here would compare grobs and environments and call the noise a
  ## difference, or return NULL and call the absence a match.
  stop("xv_reduce: no extractor for class ", paste(class(x), collapse = "/"),
       " -- add one rather than letting it fall through")
}

## ---- captured step --------------------------------------------------------

## Reduces one xv_capture() record. Warnings and messages are first-class: a
## regression that adds, drops or downgrades a warning must be visible, and
## delta 8 (ab76c49) is exactly a warning path that errored.
xv_reduce_step <- function(rec, canon){
  out <- list(
    ok = isTRUE(rec$ok),
    ## Only conditionMessage, never conditionCall: extracting a guard into an
    ## assert_* helper legitimately moves the call frame without changing behavior.
    error_message = if(!isTRUE(rec$ok)) xv_scrub(rec$value) else NULL,
    error_classes = rec$error_classes,
    warnings = sort(xv_scrub(rec$warnings)),
    messages = sort(xv_scrub(rec$messages))
  )
  out$value <- if(isTRUE(rec$ok)) xv_reduce_any(rec$value, canon) else NULL
  out
}

## ---- flattening -----------------------------------------------------------

## Flattens a reduced payload to a path -> leaf dictionary. The comparator then
## diffs two dictionaries, which is what makes a difference report name an exact
## path like "STenrich/value/elements/alpha/columns/adj_p_value".
xv_flatten <- function(x, prefix = "", out = new.env(parent = emptyenv())){
  if(is.null(x)){
    assign(prefix, "<NULL>", envir = out)
  } else if(is.atomic(x)){
    assign(prefix, x, envir = out)
  } else if(is.list(x)){
    if(!length(x)) assign(prefix, "<empty-list>", envir = out)
    nms <- names(x)
    if(is.null(nms)) nms <- paste0("[", seq_along(x), "]")
    for(i in seq_along(x)){
      key <- if(nzchar(prefix)) paste(prefix, nms[i], sep = "/") else nms[i]
      xv_flatten(x[[i]], key, out)
    }
  } else {
    stop("xv_flatten: unexpected type at ", prefix, ": ", paste(class(x), collapse = "/"))
  }
  out
}

xv_flatten_list <- function(x){
  e <- xv_flatten(x)
  ks <- ls(e, all.names = TRUE)
  o <- lapply(ks, function(k) get(k, envir = e))
  names(o) <- ks
  o
}

## ---- the leak assertion ---------------------------------------------------

## The single most important correctness property of the harness. If a payload
## contains an environment, closure, S4 object or external pointer, the
## comparator would have to load some version of spatialGE to read it -- which
## defeats the isolation the whole design depends on.
xv_assert_plain <- function(x, path = "<root>"){
  if(is.null(x)) return(invisible(TRUE))
  if(is.function(x))    stop("payload leak: function at ", path)
  if(is.environment(x)) stop("payload leak: environment at ", path)
  if(typeof(x) == "externalptr") stop("payload leak: externalptr at ", path)
  if(isS4(x))           stop("payload leak: S4 ", paste(class(x), collapse="/"), " at ", path)
  if(inherits(x, c("formula", "call", "name", "expression")))
    stop("payload leak: language object at ", path)
  if(is.list(x)){
    nms <- names(x); if(is.null(nms)) nms <- paste0("[", seq_along(x), "]")
    for(i in seq_along(x)) xv_assert_plain(x[[i]], paste(path, nms[i], sep = "/"))
  }
  invisible(TRUE)
}
