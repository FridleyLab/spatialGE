#' @importFrom methods .hasSlot is
#' @importFrom stats aggregate cor setNames median sd
#' @importFrom FNN get.knn
#' @importFrom circlize colorRamp2
NULL

utils::globalVariables(c(
  "Step", "clusterA", "clusterB", "dominant_pathway", "similarity",
  "STniche_DTC", "cell_type", "n", "pct", "mean_Zi", "total",
  "GeneSet", "niche", "MeanZ", "label", "xlab", "lab",
  "Boundary", "is_hl", "libkey", "Zi", "sampleA", "sampleB",
  "Pair", "Pathway", "Stratum", "Label", "Weight", "Sample"
))

#' @title Plot STniche results
#' @description Fifteen plot types: clusters, Zi, ranking, composition, noise,
#'   violin, zi_heatmap, celltype_heatmap, celltype_composition,
#'   celltype_stackedbar, cluster_annotation, cluster_matching,
#'   celltype_spatial, niche_signature, multi_sample_alluvial.
#'
#' @details
#' \strong{Why a stacked bar of cell types per niche?}
#' A niche is only biologically interpretable if its cellular make-up is known.
#' The \code{celltype_heatmap} shows the same information as a grid, which is
#' efficient for many cell types but makes it hard to judge how much of a niche
#' is occupied by each type. \code{type = "celltype_stackedbar"} draws one bar
#' per niche whose segments sum to 100\% (or to the raw cell count), so that
#' (i) the dominant cell type(s) of every niche, and therefore its purity, can
#' be read directly; (ii) niches can be compared side by side within a sample;
#' and (iii) the same niche can be compared across samples (e.g. IO-naive vs
#' IO-exposed) when several samples are supplied, which exposes sample-specific
#' composition shifts that pooled summaries hide. Because the bars are
#' compositional, a high share of one type is not evidence of enrichment
#' relative to the tissue background; interpret them together with
#' \code{celltype_heatmap} and with the niche's absolute size (printed as
#' \code{n} under each bar).
#'
#' @param x An STlist object.
#' @param samples Character vector of sample names (or indices) to plot.
#'   Default = all samples. For \code{type = "multi_sample_alluvial"} a named
#'   list of sample groups is also accepted.
#' @param type Character. One of the supported plot types (see details).
#' @param zi Matrix of Local Moran Z values (rows = gene sets,
#'   columns = "sample&&libname"). Required for several plot types.
#' @param gene_set Character. Gene set / pathway name required for
#'   \code{type} in \code{c("Zi", "ranking", "violin")}.
#' @param palette Colour palette. Can be a named character vector, a function
#'   that returns colours, or a length-3 vector/named vector for Zi gradients.
#' @param ptsize Numeric. Point size for spatial plots.
#' @param stroke Numeric. Point stroke width (0 = no stroke).
#' @param reverse_y Logical. Reverse the y-axis.
#' @param overlay_samples Logical. Overlay multiple samples in one plot when
#'   possible.
#' @param title Character. Plot title (overrides the default title).
#' @param legend_name Character. Legend title.
#' @param pearson_thresh Numeric. Pearson correlation threshold
#'   (currently unused; reserved for future matching logic).
#' @param jaccard_thresh Numeric. Jaccard index threshold
#'   (currently unused; reserved for future matching logic).
#' @param top_n_spots Integer. Number of top spots
#'   (currently unused; reserved for future ranking logic).
#' @param pie_top_k Integer. Number of top gene sets shown in pie charts.
#' @param out_dir Character. Output directory
#'   (currently unused; reserved for optional file export).
#' @param highlight_types Cell types to highlight (only for
#'   \code{type = "celltype_spatial"}). Either an unnamed character vector, or
#'   a \emph{named} vector whose names are display labels and whose values are
#'   regular expressions, e.g. \code{c(CD8 = "CD8", M2 = "M2|CD163")}. Patterns
#'   are matched with \code{grepl(..., ignore.case = TRUE)} against
#'   \code{cell_type}.
#' @param highlight_cols Named character of colours per highlighted type
#'   (names must equal the display labels).
#' @param pool_samples Logical; for \code{type = "celltype_stackedbar"}, pool
#'   all selected samples into one set of bars (default \code{FALSE}: one panel
#'   per sample).
#' @param bar_scale \code{"percent"} (default; each bar sums to 100) or
#'   \code{"count"} (absolute number of cells per niche) for
#'   \code{type = "celltype_stackedbar"}.
#' @param top_n_celltypes Optional integer; keep the \code{top_n_celltypes}
#'   most abundant cell types and collapse the rest into "Other"
#'   (\code{type = "celltype_stackedbar"}).
#' @param min_cells Minimum number of cells a niche needs to be drawn
#'   (\code{type = "celltype_stackedbar"}).
#' @param label_min_pct Segments at least this large (in percent) are labelled
#'   with their value (\code{type = "celltype_stackedbar"}).
#'
#' @return A \code{ggplot} object, a named list of them (one per sample), a
#'   \code{ComplexHeatmap} object/list, or a data frame, depending on
#'   \code{type}. \code{"celltype_stackedbar"} and \code{"niche_signature"}
#'   return a single \code{ggplot}.
#'
#' @export
plot_stniche <- function(x = NULL,
                         samples = NULL,
                         type = c("clusters", "Zi", "ranking",
                                  "composition", "noise", "violin",
                                  "zi_heatmap",
                                  "celltype_heatmap", "celltype_composition",
                                  "cluster_annotation", "cluster_matching",
                                  "celltype_spatial",
                                  "celltype_stackedbar",
                                  "niche_signature",
                                  "multi_sample_alluvial"),
                         zi = NULL,
                         gene_set = NULL,
                         palette = NULL,
                         ptsize = 0.9,
                         stroke = 0,
                         reverse_y = FALSE,
                         overlay_samples = FALSE,
                         title = NULL,
                         legend_name = NULL,
                         pearson_thresh = 0.60,
                         jaccard_thresh = 0.15,
                         top_n_spots    = 200L,
                         pie_top_k      = 3L,
                         out_dir        = NULL,
                         highlight_types = NULL,
                         highlight_cols  = NULL,
                         pool_samples    = FALSE,
                         bar_scale       = c("percent", "count"),
                         top_n_celltypes = NULL,
                         min_cells       = 1L,
                         label_min_pct   = 5) {
  if (is.null(x) || !methods::is(x, "STlist"))
    stop("x must be an STlist.", call. = FALSE)
  type <- match.arg(type)
  
  if (is.null(samples)) {
    samples <- names(x@spatial_meta)
  } else if (is.numeric(samples)) {
    samples <- names(x@spatial_meta)[samples]
  } else if (is.list(samples)) {
    # named list of groups is allowed for multi_sample_alluvial
    if (type != "multi_sample_alluvial") {
      samples <- unique(unlist(samples))
      samples <- intersect(samples, names(x@spatial_meta))
    }
  } else {
    samples <- intersect(samples, names(x@spatial_meta))
  }
  
  if (length(samples) == 0 && type != "multi_sample_alluvial")
    stop("No matching samples in the STlist.", call. = FALSE)
  
  if (type %in% c("Zi", "ranking", "violin", "zi_heatmap",
                  "cluster_annotation", "cluster_matching", "niche_signature",
                  "multi_sample_alluvial")) {
    if (is.null(zi))
      stop("`zi` is required for type = '", type, "'.", call. = FALSE)
  }
  if (type %in% c("Zi", "ranking", "violin")) {
    if (is.null(gene_set) || !nzchar(gene_set))
      stop("`gene_set` is required for type = '", type, "'.", call. = FALSE)
    if (!gene_set %in% rownames(zi))
      stop("gene_set '", gene_set, "' not found in rownames(zi).", call. = FALSE)
  }
  if (type %in% c("celltype_heatmap", "celltype_composition",
                  "celltype_stackedbar", "celltype_spatial")) {
    has_ct <- vapply(samples, function(s) {
      "cell_type" %in% colnames(x@spatial_meta[[s]]) &&
        any(!is.na(x@spatial_meta[[s]]$cell_type))
    }, logical(1))
    if (!any(has_ct))
      stop("No sample has a non-empty 'cell_type' column; ",
           "type = '", type, "' requires cell_type.",
           call. = FALSE)
    samples <- samples[has_ct]
  }
  if (type == "celltype_spatial" &&
      (is.null(highlight_types) || !length(highlight_types)))
    stop("`highlight_types` is required for type = 'celltype_spatial'.",
         call. = FALSE)
  
  # ------------------------------------------------------------------
  # helpers
  # ------------------------------------------------------------------
  sort_clusters <- function(v) {
    v <- unique(as.character(v))
    n <- suppressWarnings(as.integer(sub("^[A-Za-z_]+", "", v)))
    if (all(!is.na(n))) v[order(n)] else sort(v)
  }
  pal_for <- function(lv) {
    lv <- unique(as.character(lv))
    if (anyDuplicated(lv)) lv <- make.unique(lv)
    n <- length(lv)
    if (!is.null(palette)) {
      if (is.function(palette)) return(stats::setNames(palette(n), lv))
      if (is.character(palette) && !is.null(names(palette))) {
        miss <- setdiff(lv, names(palette))
        extra <- if (length(miss) > 0)
          stats::setNames(grDevices::hcl.colors(length(miss), "Dynamic"), miss)
        else character(0)
        return(c(palette[intersect(names(palette), lv)], extra)[lv])
      }
      if (length(palette) == n) return(stats::setNames(palette, lv))
    }
    stats::setNames(grDevices::hcl.colors(n, "Dynamic"), lv)
  }
  pal_zi_gradient <- function() {
    defaults <- c(low = "#0571b0", mid = "grey92", high = "#ca0020")
    if (is.null(palette)) return(defaults)
    if (is.character(palette) &&
        all(c("low", "mid", "high") %in% names(palette)))
      return(palette[c("low", "mid", "high")])
    if (length(palette) == 3 && is.null(names(palette))) {
      names(palette) <- c("low", "mid", "high")
      return(palette)
    }
    defaults
  }
  spatial_void_theme <- function(base_size = 11) {
    ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        plot.title        = ggplot2::element_text(
          face = "bold", size = base_size + 1, hjust = 0,
          margin = ggplot2::margin(b = 6)),
        legend.title      = ggplot2::element_text(
          face = "bold", size = base_size),
        legend.text       = ggplot2::element_text(size = base_size - 1),
        legend.key.size   = grid::unit(0.9, "lines"),
        legend.background = ggplot2::element_rect(
          colour = NA, fill = "transparent"),
        plot.margin       = ggplot2::margin(6, 8, 6, 6),
        strip.text        = ggplot2::element_text(
          face = "bold", size = base_size - 1),
        strip.background  = ggplot2::element_rect(
          colour = NA, fill = "grey95"))
  }
  add_boundary_flag <- function(d) {
    if (nrow(d) < 3) {
      d$Boundary <- "Interior"
      return(d)
    }
    nn <- FNN::get.knn(as.matrix(d[, c("xpos", "ypos")]), k = 1)
    nn_idx <- as.integer(nn$nn.index[, 1])
    d$Boundary <- ifelse(d$STniche_DTC != d$STniche_DTC[nn_idx],
                         "Boundary", "Interior")
    d
  }
  
  # ------------------------------------------------------------------
  # multi_sample_alluvial
  # ------------------------------------------------------------------
  if (type == "multi_sample_alluvial") {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
      stop("Package 'ggalluvial' is required for type = 'multi_sample_alluvial'.",
           call. = FALSE)
    }
    
    if (is.list(samples)) {
      sample_groups <- samples
    } else {
      sample_groups <- as.list(samples)
      names(sample_groups) <- samples
    }
    
    group_names <- names(sample_groups)
    profiles <- list()
    
    for (grp in group_names) {
      fovs <- sample_groups[[grp]]
      fov_profiles <- list()
      
      for (f in fovs) {
        if (!f %in% names(x@spatial_meta)) next
        d <- x@spatial_meta[[f]]
        keys <- paste0(f, "&&", d$libname)
        keys <- intersect(keys, colnames(zi))
        if (length(keys) < 10) next
        
        clus <- as.character(d$STniche_DTC)
        top_clus <- names(sort(table(clus), decreasing = TRUE))[1:8]
        valid_idx <- which(clus %in% top_clus)
        if (length(valid_idx) < 10) next
        
        prof <- sapply(sort(unique(clus[valid_idx])), function(cl) {
          idx <- valid_idx[clus[valid_idx] == cl]
          rowMeans(zi[, keys[idx], drop = FALSE], na.rm = TRUE)
        })
        fov_profiles[[f]] <- prof
      }
      
      if (length(fov_profiles) == 0) next
      
      all_clusters <- unique(unlist(lapply(fov_profiles, colnames)))
      consensus <- list()
      for (cl in all_clusters) {
        mat <- do.call(cbind, lapply(fov_profiles, function(p) {
          if (cl %in% colnames(p)) p[, cl] else rep(NA, nrow(zi))
        }))
        if (ncol(mat) > 0) consensus[[cl]] <- rowMeans(mat, na.rm = TRUE)
      }
      profiles[[grp]] <- do.call(cbind, consensus)
    }
    
    if (length(profiles) < 2)
      stop("Need at least 2 valid groups with data.", call. = FALSE)
    
    cosine_sim <- function(a, b) {
      if (any(is.na(a)) || any(is.na(b))) return(0)
      sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
    }
    
    match_list <- list()
    valid_groups <- names(profiles)
    for (i in seq_len(length(valid_groups) - 1)) {
      gA <- valid_groups[i]
      gB <- valid_groups[i + 1]
      prof_A <- profiles[[gA]]
      prof_B <- profiles[[gB]]
      
      for (cA in colnames(prof_A)) {
        for (cB in colnames(prof_B)) {
          sim <- cosine_sim(prof_A[, cA], prof_B[, cB])
          mean_profile <- (prof_A[, cA] + prof_B[, cB]) / 2
          max_score <- max(mean_profile, na.rm = TRUE)
          rich_pws <- names(which(mean_profile >= max_score * 0.8))
          if (length(rich_pws) == 0)
            rich_pws <- names(which.max(mean_profile))
          
          for (pw in rich_pws) {
            match_list[[length(match_list) + 1]] <- data.frame(
              Step = i, sampleA = gA, sampleB = gB,
              clusterA = cA, clusterB = cB,
              similarity = sim, dominant_pathway = pw,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    match_df <- do.call(rbind, match_list)
    
    top_matches <- match_df %>%
      dplyr::group_by(.data$Step, .data$clusterA, .data$dominant_pathway) %>%
      dplyr::slice_max(.data$similarity, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    plot_list <- list()
    for (i in unique(top_matches$Step)) {
      m <- top_matches[top_matches$Step == i, ]
      df_A <- data.frame(
        Sample = m$sampleA,
        Stratum = paste0(m$sampleA, "_", m$clusterA),
        Label = m$clusterA,
        Weight = m$similarity,
        Pair = paste0(i, "_", seq_len(nrow(m))),
        Pathway = m$dominant_pathway,
        stringsAsFactors = FALSE
      )
      df_B <- data.frame(
        Sample = m$sampleB,
        Stratum = paste0(m$sampleB, "_", m$clusterB),
        Label = m$clusterB,
        Weight = m$similarity,
        Pair = paste0(i, "_", seq_len(nrow(m))),
        Pathway = m$dominant_pathway,
        stringsAsFactors = FALSE
      )
      plot_list[[as.character(i)]] <- rbind(df_A, df_B)
    }
    plot_df <- do.call(rbind, plot_list)
    plot_df <- plot_df[!duplicated(plot_df[, c("Sample", "Stratum", "Pair")]), ]
    plot_df$Sample <- factor(plot_df$Sample, levels = valid_groups)
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data$Sample, stratum = .data$Stratum,
                   alluvium = .data$Pair, y = .data$Weight,
                   label = .data$Label)) +
      ggalluvial::geom_flow(
        stat = "alluvium", lode.guidance = "frontback",
        ggplot2::aes(fill = .data$Pathway),
        alpha = 0.7, curve_type = "cubic", width = 0.15) +
      ggalluvial::geom_stratum(
        alpha = 1, width = 0.15,
        colour = "black", linewidth = 0.5, fill = "grey90") +
      ggalluvial::stat_stratum(geom = "text", size = 3, fontface = "bold") +
      ggplot2::scale_fill_viridis_d(option = "turbo", name = "Shared Pathway") +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::labs(
        title = if (is.null(title)) "Multi-Sample Cluster Correspondence" else title,
        subtitle = "STniche domains tracked sequentially",
        y = "Similarity (Cosine)", x = NULL) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          face = "bold", size = 12, colour = "black",
          angle = 45, hjust = 1),
        legend.position = "right",
        plot.margin = ggplot2::margin(10, 10, 10, 10))
    return(p)
  }
  
  # ------------------------------------------------------------------
  # niche_signature
  # ------------------------------------------------------------------
  if (type == "niche_signature") {
    meta <- do.call(rbind, lapply(samples, function(s) {
      d <- x@spatial_meta[[s]]
      if (is.null(d) || !"STniche_DTC" %in% colnames(d)) return(NULL)
      data.frame(
        key   = paste0(s, "&&", d$libname),
        niche = as.character(d$STniche_DTC),
        stringsAsFactors = FALSE)
    }))
    if (is.null(meta))
      stop("No sample has an STniche_DTC column.", call. = FALSE)
    meta <- meta[meta$key %in% colnames(zi) & !is.na(meta$niche), , drop = FALSE]
    if (nrow(meta) == 0)
      stop("No spots/cells could be matched between the STlist and `zi`.",
           call. = FALSE)
    
    z   <- zi[, meta$key, drop = FALSE]
    idx <- split(seq_len(nrow(meta)), meta$niche)
    m   <- vapply(idx, function(i)
      rowMeans(z[, i, drop = FALSE], na.rm = TRUE),
      numeric(nrow(z)))
    if (!is.matrix(m))
      m <- matrix(m, nrow = nrow(z),
                  dimnames = list(rownames(z), names(idx)))
    rownames(m) <- rownames(z)
    
    sig <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
    names(sig) <- c("GeneSet", "niche", "MeanZ")
    n_tab   <- table(meta$niche)
    sig$n   <- as.integer(n_tab[sig$niche])
    lev     <- sort_clusters(sig$niche)
    sig$label <- factor(
      paste0(sig$niche, "\n(n=", sig$n, ")"),
      levels = paste0(lev, "\n(n=", as.integer(n_tab[lev]), ")"))
    sig$GeneSet <- factor(sig$GeneSet, levels = rev(rownames(z)))
    
    grad <- pal_zi_gradient()
    p <- ggplot2::ggplot(
      sig, ggplot2::aes(x = .data$label, y = .data$GeneSet,
                        fill = .data$MeanZ)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f", .data$MeanZ)), size = 2.6) +
      ggplot2::scale_fill_gradient2(
        low = grad[["low"]], mid = "white", high = grad[["high"]],
        midpoint = 0,
        name = if (is.null(legend_name)) "Mean Zi" else legend_name) +
      ggplot2::labs(
        x = NULL, y = NULL,
        title = if (is.null(title))
          "Gene-set signature of each niche" else title) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid  = ggplot2::element_blank(),
        plot.title  = ggplot2::element_text(face = "bold", size = 11),
        axis.text.x = ggplot2::element_text(size = 8))
    return(p)
  }
  
  # ------------------------------------------------------------------
  # Build per_sample list
  # ------------------------------------------------------------------
  per_sample <- list()
  for (s in samples) {
    if (!"STniche_DTC" %in% colnames(x@spatial_meta[[s]])) {
      warning("Sample ", s, " has no STniche_DTC column; skipping.")
      next
    }
    d <- x@spatial_meta[[s]]
    d$STniche_DTC <- as.character(d$STniche_DTC)
    d$sample      <- s
    if (type %in% c("Zi", "ranking", "violin", "zi_heatmap")) {
      d$libkey <- paste0(s, "&&", d$libname)
      common   <- intersect(d$libkey, colnames(zi))
      if (length(common) < 5) {
        warning("Sample ", s, ": fewer than 5 Zi values matched; skipping.")
        next
      }
      if (type %in% c("Zi", "ranking", "violin")) {
        zi_df <- data.frame(
          libkey = common,
          Zi     = as.numeric(zi[gene_set, common]),
          stringsAsFactors = FALSE)
        d <- merge(d, zi_df, by = "libkey", all.x = FALSE, sort = FALSE)
        if (nrow(d) == 0) next
      } else {
        d <- d[d$libkey %in% common, , drop = FALSE]
        if (nrow(d) < 5) next
      }
    }
    per_sample[[s]] <- d
  }
  if (length(per_sample) == 0)
    return(if (overlay_samples) NULL else list())
  
  # ------------------------------------------------------------------
  # cluster_annotation
  # ------------------------------------------------------------------
  if (type == "cluster_annotation") {
    out <- list()
    for (s in names(per_sample)) {
      d <- per_sample[[s]]
      d$libkey <- paste0(s, "&&", d$libname)
      d <- d[d$libkey %in% colnames(zi), , drop = FALSE]
      if (nrow(d) < 5) next
      lv <- sort_clusters(d$STniche_DTC)
      
      zi_long <- do.call(rbind, lapply(rownames(zi), function(gs) {
        data.frame(
          STniche_DTC = d$STniche_DTC,
          gene_set    = gs,
          Zi          = as.numeric(zi[gs, d$libkey]),
          stringsAsFactors = FALSE)
      }))
      zi_long <- zi_long[is.finite(zi_long$Zi), , drop = FALSE]
      if (nrow(zi_long) == 0) next
      
      zi_pos <- zi_long
      zi_pos$Zi <- pmax(zi_pos$Zi, 0)
      zi_sums <- stats::aggregate(Zi ~ STniche_DTC, data = zi_pos, FUN = sum)
      zi_pos$pct <- 100 * zi_pos$Zi /
        zi_sums$Zi[match(zi_pos$STniche_DTC, zi_sums$STniche_DTC)]
      zi_pos$pct[!is.finite(zi_pos$pct)] <- 0
      
      zi_top <- do.call(rbind, lapply(lv, function(cl) {
        sub <- zi_pos[zi_pos$STniche_DTC == cl, , drop = FALSE]
        if (nrow(sub) == 0) return(NULL)
        sub  <- sub[order(-sub$pct), , drop = FALSE]
        keep <- sub[seq_len(min(pie_top_k, nrow(sub))),
                    c("STniche_DTC", "gene_set", "pct")]
        other <- 100 - sum(keep$pct)
        if (other > 0.5)
          keep <- rbind(keep, data.frame(
            STniche_DTC = cl, gene_set = "other", pct = other,
            stringsAsFactors = FALSE))
        keep
      }))
      zi_top$STniche_DTC <- factor(zi_top$STniche_DTC, levels = lv)
      zi_top$gene_set    <- factor(
        zi_top$gene_set,
        levels = c(setdiff(unique(zi_top$gene_set), "other"), "other"))
      
      ncol_grid <- min(4, length(lv))
      p_zi <- ggplot2::ggplot(
        zi_top, ggplot2::aes(x = "", y = .data$pct, fill = .data$gene_set)) +
        ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.3) +
        ggplot2::coord_polar(theta = "y", start = 0) +
        ggplot2::facet_wrap(~ STniche_DTC, ncol = ncol_grid) +
        ggplot2::scale_fill_manual(
          values = stats::setNames(
            grDevices::hcl.colors(length(levels(zi_top$gene_set)), "Dynamic"),
            levels(zi_top$gene_set)),
          name = if (is.null(legend_name)) "Gene set" else legend_name) +
        ggplot2::labs(
          title = if (is.null(title))
            sprintf("Per-cluster gene-set Zi contribution | %s", s) else title,
          subtitle = sprintf(
            "Each pie = %% of the cluster's positive Zi (top %d gene sets + other)",
            pie_top_k),
          x = NULL, y = NULL) +
        ggplot2::theme_void(base_size = 11) +
        ggplot2::theme(
          strip.text      = ggplot2::element_text(face = "bold", size = 11),
          legend.position = "right",
          plot.title      = ggplot2::element_text(face = "bold", size = 13),
          plot.subtitle   = ggplot2::element_text(size = 10),
          plot.margin     = ggplot2::margin(6, 8, 6, 6))
      out[[paste0(s, "_gene_sets")]] <- p_zi
      
      if ("cell_type" %in% colnames(d) && any(!is.na(d$cell_type))) {
        dd <- d[!is.na(d$cell_type), , drop = FALSE]
        dd$cell_type <- as.character(dd$cell_type)
        ct_counts <- dd %>%
          dplyr::count(.data$STniche_DTC, .data$cell_type, name = "n") %>%
          dplyr::group_by(.data$STniche_DTC) %>%
          dplyr::mutate(pct = 100 * .data$n / sum(.data$n)) %>%
          dplyr::ungroup()
        ct_counts$STniche_DTC <- factor(ct_counts$STniche_DTC, levels = lv)
        
        p_ct <- ggplot2::ggplot(
          ct_counts,
          ggplot2::aes(x = "", y = .data$pct, fill = .data$cell_type)) +
          ggplot2::geom_col(width = 1, colour = "white", linewidth = 0.3) +
          ggplot2::coord_polar(theta = "y", start = 0) +
          ggplot2::facet_wrap(~ STniche_DTC, ncol = ncol_grid) +
          ggplot2::scale_fill_manual(
            values = stats::setNames(
              grDevices::hcl.colors(length(unique(ct_counts$cell_type)), "Dynamic"),
              unique(ct_counts$cell_type)),
            name = if (is.null(legend_name)) "Cell type" else legend_name) +
          ggplot2::labs(
            title = if (is.null(title))
              sprintf("Per-cluster cell-type composition | %s", s) else title,
            subtitle = "Each pie = % of cells per cell type within the cluster",
            x = NULL, y = NULL) +
          ggplot2::theme_void(base_size = 11) +
          ggplot2::theme(
            strip.text      = ggplot2::element_text(face = "bold", size = 11),
            legend.position = "right",
            plot.title      = ggplot2::element_text(face = "bold", size = 13),
            plot.subtitle   = ggplot2::element_text(size = 10),
            plot.margin     = ggplot2::margin(6, 8, 6, 6))
        out[[paste0(s, "_cell_types")]] <- p_ct
      }
    }
    return(out)
  }
  
  # ------------------------------------------------------------------
  # celltype_spatial
  # ------------------------------------------------------------------
  if (type == "celltype_spatial") {
    # normalise highlight_types into patterns + labels
    if (is.null(names(highlight_types))) {
      hl_lab <- highlight_types
      hl_pat <- highlight_types
    } else {
      hl_lab <- names(highlight_types)
      hl_pat <- unname(highlight_types)
    }
    if (is.null(highlight_cols)) {
      hl_cols <- stats::setNames(
        grDevices::hcl.colors(length(hl_lab), "Dynamic"), hl_lab)
    } else {
      hl_cols <- highlight_cols
    }
    
    out <- list()
    for (s in names(per_sample)) {
      d <- per_sample[[s]]
      d <- d[!is.na(d$cell_type), , drop = FALSE]
      if (nrow(d) < 5) next
      d$cell_type <- as.character(d$cell_type)
      d$is_hl     <- rep("Other", nrow(d))
      for (k in seq_along(hl_pat)) {
        d$is_hl[grepl(hl_pat[k], d$cell_type, ignore.case = TRUE)] <- hl_lab[k]
      }
      d$is_hl <- factor(d$is_hl, levels = c("Other", hl_lab))
      cols <- c(stats::setNames("#C8C8C8", "Other"), hl_cols)
      
      p <- ggplot2::ggplot(d) +
        ggplot2::geom_point(
          data = d[d$is_hl == "Other", ],
          ggplot2::aes(x = .data$xpos, y = .data$ypos),
          size = ptsize * 0.6, colour = "#D0D0D0", alpha = 0.7) +
        ggplot2::geom_point(
          data = d[d$is_hl != "Other", ],
          ggplot2::aes(x = .data$xpos, y = .data$ypos, colour = .data$is_hl),
          size = ptsize, alpha = 0.9) +
        ggplot2::scale_colour_manual(
          values = cols,
          name = if (is.null(legend_name)) "Cell type" else legend_name,
          drop = FALSE) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::labs(title = if (is.null(title))
          sprintf("Cell types - %s", s) else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      out[[s]] <- p
    }
    return(out)
  }
  
  # ------------------------------------------------------------------
  # shared geometric helpers
  # ------------------------------------------------------------------
  dot_layer <- function() {
    if (stroke > 0)
      ggplot2::geom_point(
        shape = 21, size = ptsize, colour = "white", stroke = stroke,
        mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                               fill = .data$STniche_DTC))
    else
      ggplot2::geom_point(
        size = ptsize,
        mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                               colour = .data$STniche_DTC))
  }
  dot_layer_zi <- function() {
    if (stroke > 0)
      ggplot2::geom_point(
        shape = 21, size = ptsize, colour = "white", stroke = stroke,
        mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                               fill = .data$Zi))
    else
      ggplot2::geom_point(
        size = ptsize,
        mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                               colour = .data$Zi))
  }
  
  # ------------------------------------------------------------------
  # overlay_samples shortcuts
  # ------------------------------------------------------------------
  if (overlay_samples && type %in% c("clusters", "Zi", "noise")) {
    if (type == "noise") {
      per_sample_noise <- lapply(names(per_sample), function(s)
        add_boundary_flag(per_sample[[s]]))
      names(per_sample_noise) <- names(per_sample)
      map_dat <- do.call(rbind, per_sample_noise)
    } else {
      map_dat <- do.call(rbind, per_sample)
    }
    map_dat$sample <- factor(map_dat$sample, levels = names(per_sample))
    
    if (type == "clusters") {
      lv <- sort_clusters(map_dat$STniche_DTC)
      map_dat$STniche_DTC <- factor(map_dat$STniche_DTC, levels = lv)
      cols <- pal_for(lv)
      p <- ggplot2::ggplot(map_dat) + dot_layer() +
        (if (stroke > 0)
          ggplot2::scale_fill_manual(
            values = cols,
            name = if (is.null(legend_name)) "STniche cluster" else legend_name)
         else
           ggplot2::scale_colour_manual(
             values = cols,
             name = if (is.null(legend_name)) "STniche cluster" else legend_name)) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::facet_wrap(~ sample,
                            labeller = ggplot2::label_wrap_gen(width = 25)) +
        ggplot2::labs(title = if (is.null(title)) "STniche clusters" else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      return(p)
    }
    if (type == "Zi") {
      grad <- pal_zi_gradient()
      p <- ggplot2::ggplot(map_dat) + dot_layer_zi() +
        (if (stroke > 0)
          ggplot2::scale_fill_gradient2(
            low = grad[["low"]], mid = grad[["mid"]],
            high = grad[["high"]], midpoint = 0,
            name = if (is.null(legend_name)) "Local Moran Z" else legend_name)
         else
           ggplot2::scale_colour_gradient2(
             low = grad[["low"]], mid = grad[["mid"]],
             high = grad[["high"]], midpoint = 0,
             name = if (is.null(legend_name)) "Local Moran Z" else legend_name)) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::facet_wrap(~ sample,
                            labeller = ggplot2::label_wrap_gen(width = 25)) +
        ggplot2::labs(title = if (is.null(title))
          sprintf("%s - Local Moran Z", gene_set) else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      return(p)
    }
    # noise
    p <- ggplot2::ggplot(map_dat) +
      ggplot2::geom_point(
        size = ptsize,
        mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                               colour = .data$Boundary)) +
      ggplot2::scale_colour_manual(
        values = c(Boundary = "#D62728", Interior = "grey80"),
        name = if (is.null(legend_name)) "Spot type" else legend_name) +
      ggplot2::coord_equal() + spatial_void_theme() +
      ggplot2::facet_wrap(~ sample,
                          labeller = ggplot2::label_wrap_gen(width = 25)) +
      ggplot2::labs(title = if (is.null(title))
        "STniche cluster boundaries" else title)
    if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
    return(p)
  }
  
  # ------------------------------------------------------------------
  # zi_heatmap
  # ------------------------------------------------------------------
  if (type == "zi_heatmap") {
    ht_list <- list()
    for (s in names(per_sample)) {
      d <- per_sample[[s]]
      d <- d[!is.na(d$STniche_DTC), , drop = FALSE]
      if (nrow(d) < 5) next
      lv    <- sort_clusters(d$STniche_DTC)
      cols  <- pal_for(lv)
      ord   <- order(match(d$STniche_DTC, lv))
      d_ord <- d[ord, , drop = FALSE]
      Zi_u  <- zi[, d_ord$libkey, drop = FALSE]
      keep  <- apply(Zi_u, 1, function(v) any(is.finite(v)))
      Zi_u  <- Zi_u[keep, , drop = FALSE]
      if (nrow(Zi_u) == 0) next
      
      ha <- ComplexHeatmap::HeatmapAnnotation(
        STniche_DTC = factor(d_ord$STniche_DTC, levels = lv),
        col = list(STniche_DTC = cols),
        show_annotation_name = TRUE)
      ht <- ComplexHeatmap::Heatmap(
        Zi_u,
        name = if (is.null(legend_name)) "Local Moran Z" else legend_name,
        top_annotation = ha,
        show_column_names = FALSE,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_title = if (is.null(title))
          sprintf("Local Moran Z | %s", s) else title)
      ht_list[[s]] <- ht
    }
    if (length(ht_list) == 0) return(NULL)
    if (length(ht_list) == 1) return(ht_list[[1]])
    return(ht_list)
  }
  
  # ------------------------------------------------------------------
  # celltype_heatmap
  # ------------------------------------------------------------------
  if (type == "celltype_heatmap") {
    ht_list <- list()
    for (s in names(per_sample)) {
      d_ct <- per_sample[[s]]
      d_ct <- d_ct[!is.na(d_ct$cell_type), ]
      if (nrow(d_ct) == 0) next
      tab   <- table(d_ct$STniche_DTC, d_ct$cell_type)
      lv_ct <- sort_clusters(rownames(tab))
      tab   <- tab[lv_ct, , drop = FALSE]
      mat   <- sweep(tab, 1, rowSums(tab), "/")
      mat[is.na(mat)] <- 0
      
      ht <- ComplexHeatmap::Heatmap(
        mat * 100,
        name = "% of cluster",
        column_title = if (is.null(title))
          sprintf("Cell types within STniche clusters | %s", s) else title,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        col = circlize::colorRamp2(c(0, 50, 100),
                                   c("white", "#FDB863", "#B2182B")),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if (mat[i, j] >= 0.05)
            grid::grid.text(sprintf("%.0f", mat[i, j] * 100), x, y,
                            gp = grid::gpar(fontsize = 8))
        })
      ht_list[[s]] <- ht
    }
    if (length(ht_list) == 0) return(NULL)
    if (length(ht_list) == 1) return(ht_list[[1]])
    return(ht_list)
  }
  
  # ------------------------------------------------------------------
  # celltype_stackedbar
  # ------------------------------------------------------------------
  if (type == "celltype_stackedbar") {
    bar_scale <- match.arg(bar_scale)
    all_d <- do.call(rbind, lapply(names(per_sample), function(s) {
      d <- per_sample[[s]]
      d <- d[!is.na(d$cell_type) & !is.na(d$STniche_DTC),
             c("sample", "STniche_DTC", "cell_type"), drop = FALSE]
      d$cell_type <- as.character(d$cell_type)
      d
    }))
    if (is.null(all_d) || nrow(all_d) == 0) return(NULL)
    if (isTRUE(pool_samples)) all_d$sample <- "All samples (pooled)"
    
    if (!is.null(top_n_celltypes)) {
      tot  <- sort(table(all_d$cell_type), decreasing = TRUE)
      keep <- names(tot)[seq_len(min(as.integer(top_n_celltypes), length(tot)))]
      all_d$cell_type[!all_d$cell_type %in% keep] <- "Other"
    }
    
    counts <- all_d %>%
      dplyr::count(.data$sample, .data$STniche_DTC, .data$cell_type, name = "n") %>%
      dplyr::group_by(.data$sample, .data$STniche_DTC) %>%
      dplyr::mutate(total = sum(.data$n), pct = 100 * .data$n / .data$total) %>%
      dplyr::ungroup()
    counts <- counts[counts$total >= min_cells, , drop = FALSE]
    if (nrow(counts) == 0) {
      warning("No niche has at least `min_cells` = ", min_cells, " cells.")
      return(NULL)
    }
    
    ct_ord <- names(sort(tapply(counts$n, counts$cell_type, sum), decreasing = TRUE))
    ct_lv  <- c(setdiff(ct_ord, "Other"), intersect("Other", ct_ord))
    ct_cols <- pal_for(ct_lv)
    if ("Other" %in% ct_lv &&
        (is.null(palette) || is.null(names(palette)) || !"Other" %in% names(palette)))
      ct_cols[["Other"]] <- "grey70"
    counts$cell_type <- factor(counts$cell_type, levels = ct_lv)
    
    s_lv <- if (isTRUE(pool_samples)) unique(counts$sample)
    else intersect(names(per_sample), unique(counts$sample))
    counts$sample <- factor(counts$sample, levels = s_lv)
    cl_lv <- sort_clusters(counts$STniche_DTC)
    key   <- unique(counts[, c("sample", "STniche_DTC", "total")])
    key   <- key[order(match(key$sample, s_lv), match(key$STniche_DTC, cl_lv)), ]
    key$xlab <- paste0(key$STniche_DTC, "\n(n=", key$total, ")")
    counts$xlab <- key$xlab[match(paste(counts$sample, counts$STniche_DTC),
                                  paste(key$sample, key$STniche_DTC))]
    counts$xlab <- factor(counts$xlab, levels = unique(key$xlab))
    
    yvar  <- if (bar_scale == "percent") "pct" else "n"
    ylab  <- if (bar_scale == "percent") "% of cells in niche" else "Number of cells"
    lab_d <- counts[counts$pct >= label_min_pct, , drop = FALSE]
    lab_d$lab <- if (bar_scale == "percent") sprintf("%.0f", lab_d$pct)
    else as.character(lab_d$n)
    
    p <- ggplot2::ggplot(
      counts,
      ggplot2::aes(x = .data$xlab, y = .data[[yvar]], fill = .data$cell_type)) +
      ggplot2::geom_col(width = 0.8, colour = "grey20", linewidth = 0.2,
                        position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::geom_text(
        data = lab_d,
        ggplot2::aes(label = .data$lab),
        position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE),
        size = 2.6, colour = "black") +
      ggplot2::scale_fill_manual(
        values = ct_cols, drop = FALSE,
        name = if (is.null(legend_name)) "Cell type" else legend_name) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(
        x = "STniche cluster", y = ylab,
        title = if (is.null(title))
          "Cell-type composition of STniche niches" else title,
        subtitle = if (bar_scale == "percent")
          "Each bar = 100% of the cells in a niche; n = cells in the niche"
        else
          "Each bar = number of cells in a niche, split by cell type") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        axis.text.x        = ggplot2::element_text(size = 8),
        plot.title         = ggplot2::element_text(face = "bold", size = 13),
        strip.text         = ggplot2::element_text(face = "bold"))
    if (length(s_lv) > 1)
      p <- p + ggplot2::facet_wrap(
        ~ sample, scales = "free_x", ncol = min(length(s_lv), 3),
        labeller = ggplot2::label_wrap_gen(width = 25))
    return(p)
  }
  
  # ------------------------------------------------------------------
  # celltype_composition
  # ------------------------------------------------------------------
  if (type == "celltype_composition") {
    out <- list()
    for (s in names(per_sample)) {
      d_ct <- per_sample[[s]]
      d_ct$cell_type <- as.character(d_ct$cell_type)
      d_ct <- d_ct[!is.na(d_ct$cell_type), ]
      if (nrow(d_ct) == 0) next
      lv_cl <- sort_clusters(d_ct$STniche_DTC)
      d_ct$STniche_DTC <- factor(d_ct$STniche_DTC, levels = lv_cl)
      
      counts <- d_ct %>%
        dplyr::count(.data$STniche_DTC, .data$cell_type, name = "n") %>%
        dplyr::group_by(.data$STniche_DTC) %>%
        dplyr::mutate(pct = 100 * .data$n / sum(.data$n)) %>%
        dplyr::ungroup()
      
      ct_levels <- unique(as.character(counts$cell_type))
      ct_levels <- ct_levels[!is.na(ct_levels)]
      if (anyDuplicated(ct_levels)) ct_levels <- make.unique(ct_levels)
      ct_cols <- stats::setNames(
        grDevices::hcl.colors(length(ct_levels), "Dynamic"), ct_levels)
      counts$cell_type <- factor(counts$cell_type, levels = ct_levels)
      
      p <- ggplot2::ggplot(
        counts,
        ggplot2::aes(x = .data$STniche_DTC, y = .data$pct,
                     fill = .data$cell_type)) +
        ggplot2::geom_col(width = 0.75, colour = "grey20", linewidth = 0.2) +
        ggplot2::scale_fill_manual(
          values = ct_cols,
          name = if (is.null(legend_name)) "Cell type" else legend_name) +
        ggplot2::labs(
          x = "STniche cluster", y = "% of cluster",
          title = if (is.null(title))
            sprintf("Cell-type composition per cluster | %s", s) else title) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor   = ggplot2::element_blank(),
          plot.title         = ggplot2::element_text(face = "bold", size = 13))
      out[[s]] <- p
    }
    return(out)
  }
  
  # ------------------------------------------------------------------
  # remaining types (clusters, Zi, ranking, composition, noise, violin)
  # ------------------------------------------------------------------
  out <- list()
  for (s in names(per_sample)) {
    d <- per_sample[[s]]
    
    if (type == "clusters") {
      lv <- sort_clusters(d$STniche_DTC)
      d$STniche_DTC <- factor(d$STniche_DTC, levels = lv)
      cols <- pal_for(lv)
      p <- ggplot2::ggplot(d) + dot_layer() +
        (if (stroke > 0)
          ggplot2::scale_fill_manual(
            values = cols,
            name = if (is.null(legend_name)) "STniche cluster" else legend_name)
         else
           ggplot2::scale_colour_manual(
             values = cols,
             name = if (is.null(legend_name)) "STniche cluster" else legend_name)) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::labs(title = if (is.null(title))
          sprintf("STniche clusters | %s", s) else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      out[[s]] <- p
      next
    }
    
    if (type == "Zi") {
      grad <- pal_zi_gradient()
      p <- ggplot2::ggplot(d) + dot_layer_zi() +
        (if (stroke > 0)
          ggplot2::scale_fill_gradient2(
            low = grad[["low"]], mid = grad[["mid"]],
            high = grad[["high"]], midpoint = 0,
            name = if (is.null(legend_name)) "Local Moran Z" else legend_name)
         else
           ggplot2::scale_colour_gradient2(
             low = grad[["low"]], mid = grad[["mid"]],
             high = grad[["high"]], midpoint = 0,
             name = if (is.null(legend_name)) "Local Moran Z" else legend_name)) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::labs(title = if (is.null(title))
          sprintf("%s | %s", gene_set, s) else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      out[[s]] <- p
      next
    }
    
    if (type == "noise") {
      d <- add_boundary_flag(d)
      p <- ggplot2::ggplot(d) +
        ggplot2::geom_point(
          size = ptsize,
          mapping = ggplot2::aes(x = .data$xpos, y = .data$ypos,
                                 colour = .data$Boundary)) +
        ggplot2::scale_colour_manual(
          values = c(Boundary = "#D62728", Interior = "grey80"),
          name = if (is.null(legend_name)) "Spot type" else legend_name) +
        ggplot2::coord_equal() + spatial_void_theme() +
        ggplot2::labs(title = if (is.null(title))
          sprintf("STniche cluster boundaries | %s", s) else title)
      if (isTRUE(reverse_y)) p <- p + ggplot2::scale_y_reverse()
      out[[s]] <- p
      next
    }
    
    if (type == "composition") {
      comp <- as.data.frame(table(d$STniche_DTC), stringsAsFactors = FALSE)
      names(comp) <- c("STniche_DTC", "n")
      comp$pct <- comp$n / sum(comp$n)
      lv <- sort_clusters(comp$STniche_DTC)
      comp$STniche_DTC <- factor(comp$STniche_DTC, levels = lv)
      cols <- pal_for(lv)
      p <- ggplot2::ggplot(
        comp,
        ggplot2::aes(x = .data$STniche_DTC, y = .data$pct,
                     fill = .data$STniche_DTC)) +
        ggplot2::geom_col(width = 0.7, colour = "grey20", linewidth = 0.25) +
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.1f%%", 100 * .data$pct)),
          vjust = -0.5, size = 3) +
        ggplot2::scale_fill_manual(values = cols, guide = "none") +
        ggplot2::scale_y_continuous(
          labels = scales::percent_format(accuracy = 1),
          expand = ggplot2::expansion(mult = c(0, 0.10))) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor   = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            colour = "grey70", fill = NA, linewidth = 0.3),
          plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
          plot.margin = ggplot2::margin(6, 8, 6, 6)) +
        ggplot2::labs(
          x = "STniche cluster", y = "Fraction of spots",
          title = if (is.null(title))
            sprintf("Cluster composition | %s", s) else title)
      out[[s]] <- p
      next
    }
    
    if (type == "violin") {
      lv <- sort_clusters(d$STniche_DTC)
      d$STniche_DTC <- factor(d$STniche_DTC, levels = lv)
      cols <- pal_for(lv)
      p <- ggplot2::ggplot(
        d,
        ggplot2::aes(x = .data$STniche_DTC, y = .data$Zi,
                     fill = .data$STniche_DTC)) +
        ggplot2::geom_violin(scale = "width", trim = TRUE,
                             colour = "grey30", linewidth = 0.25, alpha = 0.85) +
        ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.2,
                              colour = "grey20", fill = "white", alpha = 0.9) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey50", linewidth = 0.3) +
        ggplot2::scale_fill_manual(values = cols, guide = "none") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor   = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            colour = "grey70", fill = NA, linewidth = 0.3),
          plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
          plot.margin = ggplot2::margin(6, 8, 6, 6)) +
        ggplot2::labs(
          x = "STniche cluster", y = "Local Moran Z",
          title = if (is.null(title))
            sprintf("%s - Local Moran Z per cluster | %s", gene_set, s)
          else title)
      out[[s]] <- p
      next
    }
    
    # ranking (default fall-through)
    rk <- stats::aggregate(
      Zi ~ STniche_DTC, data = d,
      FUN = function(v) c(mean = mean(v, na.rm = TRUE),
                          sd   = stats::sd(v, na.rm = TRUE)))
    rk <- data.frame(
      STniche_DTC = rk$STniche_DTC,
      n           = as.vector(table(d$STniche_DTC)[rk$STniche_DTC]),
      mean_Zi     = rk$Zi[, "mean"],
      sd_Zi       = rk$Zi[, "sd"])
    rk <- rk[order(-rk$mean_Zi), , drop = FALSE]
    rk$STniche_DTC <- factor(rk$STniche_DTC, levels = rk$STniche_DTC)
    cols <- pal_for(levels(rk$STniche_DTC))
    
    p <- ggplot2::ggplot(
      rk,
      ggplot2::aes(x = .data$STniche_DTC, y = .data$mean_Zi,
                   fill = .data$STniche_DTC)) +
      ggplot2::geom_col(width = 0.72, colour = "grey20", linewidth = 0.25) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_Zi - .data$sd_Zi,
                     ymax = .data$mean_Zi + .data$sd_Zi),
        width = 0.28, colour = "grey30", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          colour = "grey50", linewidth = 0.3) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$mean_Zi)),
        vjust = ifelse(rk$mean_Zi >= 0, -0.7, 1.4), size = 3) +
      ggplot2::scale_fill_manual(values = cols, guide = "none") +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0.10, 0.10))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          colour = "grey70", fill = NA, linewidth = 0.3),
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
        plot.margin = ggplot2::margin(6, 8, 6, 6)) +
      ggplot2::labs(
        x = "STniche cluster", y = "Mean Local Moran Z",
        title = if (is.null(title))
          sprintf("%s-high cluster ranking | %s", gene_set, s) else title)
    out[[s]] <- p
  }
  out
}