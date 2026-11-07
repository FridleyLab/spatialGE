#' @title STniche (mclust version): Detect functional spatial domains with model-based clustering
#' @description Unsupervised clustering of spots/cells in one or more ST samples using
#' gene-set GSVA scores -> local Moran's I (Z) as features, then Gaussian mixture modeling
#' via mclust with BIC-based model/cluster selection.
#'
#' @details
#' Workflow:
#' 1) For each sample, compute GSVA scores for provided gene sets.
#' 2) For each gene set within each sample, compute Local Moran's I Z-score over a kNN graph (k = `knns`).
#' 3) Concatenate all samples into a single feature matrix (spots x gene-sets, entries are Z.Ii).
#' 4) Run mclust::Mclust with candidate cluster counts `ks` (e.g., 2:10). Best model & G selected by BIC.
#' 5) Write cluster labels back into `x@spatial_meta[[sample]]` under columns `stniche_dsplFalse` and `STniche_DTC` for plotting compatibility.
#'
#' @param x an STlist with normalized expression data in `x@tr_counts[[sample]]` (genes x spots)
#' @param gene_sets a *named* list of gene sets (each an atomic character vector of gene symbols)
#' @param samples vector of sample names or indices to run (default: all)
#' @param ks numeric vector of candidate cluster numbers to evaluate with BIC (default: 2:10)
#' @param knns integer k for kNN used in Local Moran's I (default: 6)
#' @param min_genes minimum overlapping genes required to compute a GSVA score (default: 1)
#' @param out_i logical; if TRUE return list(STlist, Z.Ii matrix), else STlist only
#' @param verbose integer verbosity level (0=quiet, 1=basic, 2=chatty)
#' @param cores integer number of workers for parallel scoring (Unix only; default chooses sensibly)
#'
#' @return If `out_i=TRUE`: `list(updated_STlist, Zi_matrix)` where `Zi_matrix` has rows = gene sets,
#' columns = "sample&&libname" spots. Otherwise, the updated STlist.
#'
#' @importFrom mclust Mclust mclustBIC
#' @importFrom methods is
#' @importFrom stats median
#' @importFrom mclust Mclust mclustBIC
#' @export
STniche <- function(x = NULL,
                    gene_sets = NULL,
                    samples = NULL,
                    ks = 2:10,
                    knns = 6,
                    min_genes = 1,
                    out_i = FALSE,
                    verbose = 1L,
                    cores = NULL) {
  
  # ---- Small internal helpers ----
  count_cores <- function(n) {
    cc <- 1L
    if (.Platform$OS.type == "unix") {
      total <- parallel::detectCores(logical = TRUE)
      cc <- max(1L, min(max(1L, floor(total / 2)), n))
    }
    # Respect R CMD check core limit
    if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
      cc <- min(cc, 2L)
    }
    return(as.integer(cc))
  }
  
  stop_bad <- function(msg) stop(msg, call. = FALSE)
  
  # ---- Input checks ----
  if (is.null(x) || !methods::is(x, "STlist")) stop_bad("x must be an STlist")
  if (is.null(gene_sets) || !length(gene_sets)) stop_bad("gene_sets must be a named list of gene vectors")
  if (is.null(names(gene_sets))) stop_bad("gene_sets must be a *named* list")
  if (length(x@tr_counts) < 1) stop_bad("No transformed counts found in x@tr_counts. Run transform/normalize first.")
  
  # ---- Resolve samples ----
  if (is.null(samples)) {
    samples <- names(x@spatial_meta)
  } else if (is.numeric(samples)) {
    samples <- names(x@spatial_meta)[samples]
  } else {
    samples <- intersect(samples, names(x@spatial_meta))
  }
  if (length(samples) == 0) stop_bad("None of the requested samples are present in the STlist.")
  
  # ---- Cores ----
  if (is.null(cores)) {
    cores <- count_cores(length(samples))
  } else {
    cores <- as.integer(cores)
    if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
      cores <- min(cores, 2L)
    }
  }
  if (verbose) cat(sprintf("STniche (mclust) starting with %d sample(s); using %d core(s)\n", length(samples), cores))
  
  # ---- Extract expression and coordinates ----
  trcounts_df <- x@tr_counts[samples]
  coord_dat <- lapply(samples, function(s) {
    df <- x@spatial_meta[[s]][, c("libname", "xpos", "ypos")]
    df$libname <- paste0(s, "&&", df$libname)
    df
  }) |> dplyr::bind_rows()
  
  # ---- Build (sample, gene-set) combinations ----
  combo <- data.frame(expand.grid(samples, names(gene_sets), stringsAsFactors = FALSE))
  colnames(combo) <- c("sample", "gene_set")
  combo <- combo[order(combo$sample), ]
  pw_genes <- lapply(seq_len(nrow(combo)), function(i) {
    s <- combo$sample[i]
    gs <- combo$gene_set[i]
    intersect(rownames(trcounts_df[[s]]), gene_sets[[gs]])
  })
  
  # ---- Compute GSVA scores ----
  if (verbose) cat("Computing GSVA scores...\n")
  gsea_scores <- calculate_gs_gsva_score(
    x = trcounts_df,
    combo = combo,
    gene_sets = gene_sets,
    pw_genes = pw_genes,
    min_genes = min_genes,
    cores = cores,
    verbose = (verbose > 1)
  )
  
  # Join into one long table
  gsea_scores <- lapply(seq_along(gsea_scores), function(ix) {
    sname <- names(gsea_scores)[ix]
    df <- data.frame(t(gsea_scores[[ix]]), check.names = FALSE)
    df$libname <- rownames(df)
    df$libname <- paste0(sname, "&&", df$libname)
    df
  }) |> dplyr::bind_rows()
  
  gsea_scores <- dplyr::left_join(coord_dat, gsea_scores, by = "libname")
  
  # ---- Local Moran's I (Z.Ii) ----
  if (verbose) cat(sprintf("Computing Local Moran's I with k = %d...\n", knns))
  set.seed(1L)
  
  output_list <- parallel::mclapply(seq_len(nrow(combo)), function(i) {
    s  <- combo$sample[i]
    gs <- combo$gene_set[i]
    
    df <- gsea_scores[, c("libname", "xpos", "ypos", gs)]
    df <- df[grepl(paste0("^", s, "&&"), df$libname), ]
    
    # ---- Guard against all-NA or zero-variance gene sets ----
    vals <- df[[gs]]
    if (all(is.na(vals)) || length(unique(vals[!is.na(vals)])) < 2L) {
      # Return a dummy result with NA Z.Ii so downstream code stays consistent
      out <- data.frame(
        libname = df$libname,
        Z.Ii = NA_real_,
        stringsAsFactors = FALSE
      )
      names(out)[2] <- gs
      return(out)
    }
    
    # Build kNN graph & weights
    spdf <- sf::st_as_sf(df, coords = c("xpos", "ypos"))
    knn  <- spdep::knearneigh(sf::st_coordinates(spdf), k = knns)
    nb   <- spdep::knn2nb(knn)
    lw   <- spdep::nb2listw(nb, style = "W")
    
    loc <- data.frame(spdep::localmoran(vals, lw), check.names = FALSE)
    # Sign the Local Moran Z by the deviation of the spot's value from the
    # gene set's cross-spot mean. Vanilla Z.Ii is positive for BOTH
    # "high-high clusters" and "low-low clusters" (it only measures
    # autocorrelation, not direction), which makes it unable to distinguish
    # two domains that have symmetric high/low gene-set expression.
    # Multiplying by sign(x_i - mean(x)) resolves the high-high vs low-low
    # ambiguity: high-high spots get positive signed Zi, low-low spots get
    # negative signed Zi, and the two domains separate cleanly in feature
    # space. This is the standard domain-detection variant of Local Moran's
    # I (see Getis-Ord Gi* for the equivalent for hotspots).
    vals_mean <- mean(vals, na.rm = TRUE)
    signed_Zi <- sign(vals - vals_mean) * loc[["Z.Ii"]]
    out <- tibble::add_column(loc, libname = df$libname,
                              signed_Zi = signed_Zi) |>
      dplyr::select("signed_Zi", "libname") |>
      dplyr::rename(!!gs := "signed_Zi")
    out
  }, mc.cores = cores, mc.allow.recursive = TRUE)
  
  # ---- Assemble Z.Ii matrix ----
  by_sample <- lapply(unique(combo$sample), function(s) {
    idx <- which(combo$sample == s)
    df  <- Reduce(function(a, b) merge(a, b, by = "libname", all = TRUE), output_list[idx])
    
    if (!is.null(rownames(df))) rownames(df) <- NULL
    df <- tibble::as_tibble(df)
    df <- tibble::column_to_rownames(df, "libname")
    df
  })
  
  Zi <- t(dplyr::bind_rows(by_sample))
  keep_cols <- intersect(coord_dat$libname, colnames(Zi))
  Zi <- Zi[, keep_cols, drop = FALSE]
  
  # ---- Model-based clustering ----
  if (verbose) cat("Clustering with mclust (BIC selection over ks)...\n")
  if (!is.numeric(ks) || any(ks < 1)) stop_bad("ks must be a numeric vector of candidate cluster counts")
  
  feat <- t(Zi)
  
  # Remove spots with all-NA features
  keep <- rowSums(is.na(feat)) < ncol(feat)
  feat <- feat[keep, , drop = FALSE]
  if (nrow(feat) < 5) stop_bad("Too few valid spots after filtering to run mclust")
  
  # Drop gene-set columns that carry too few non-NA values to define a
  # covariance for mclust. A column with only 1 or 2 non-NA entries is
  # effectively constant and would also drive mclust to a degenerate fit.
  keep_cols <- colSums(!is.na(feat)) >= 3
  if (sum(keep_cols) == 0) {
    stop_bad(paste0(
      "All gene sets produced all-NA Local Moran Z. ",
      "Check that gene sets have >= 2 genes and that the data has spatial signal."
    ))
  }
  dropped <- colnames(feat)[!keep_cols]
  if (length(dropped) > 0L) {
    warning(
      "Dropping gene set(s) with all-NA Local Moran Z before mclust: ",
      paste(dropped, collapse = ", "),
      call. = FALSE
    )
  }
  feat <- feat[, keep_cols, drop = FALSE]
  if (ncol(feat) < 1L) stop_bad("All gene sets have all-NA Local Moran Z; cannot run mclust")
  
  # Impute remaining NAs to 0. NA in the Zi matrix means "this gene set had
  # no valid Local Moran's I at this spot", which happens when GSVA returned
  # NA for that spot (all gene-set genes at zero, or insufficient rank
  # information). The correct value for "gene set locally absent here" is
  # Zi = 0, since a spot surrounded by a uniform absence of the gene set
  # has zero local autocorrelation for that set. Median imputation would
  # flatten all such spots to the GLOBAL median -- which for a domain-
  # restricted gene set equals the dominant-domain value, destroying the
  # very signal the clustering needs to find.
  if (anyNA(feat)) {
    feat[is.na(feat)] <- 0
  }
  
  # Final safety check
  if (anyNA(feat)) {
    stop_bad("Internal error: NA values remain in feature matrix after imputation.")
  }
  
  bic  <- mclust::mclustBIC(feat, G = ks, verbose = FALSE)
  mfit <- mclust::Mclust(feat, x = bic, G = ks, verbose = FALSE)
  
  cl_int <- as.integer(mfit$classification)
  clusters <- data.frame(
    libname = rownames(feat),
    stniche_dsplFalse = paste0("C", cl_int),
    stringsAsFactors = FALSE
  )
  clusters$STniche_DTC <- clusters$stniche_dsplFalse
  
  # ---- Merge back into STlist ----
  cl_df <- tidyr::separate_wider_delim(clusters, "libname", delim = "&&",
                                       names = c("samplename", "libname"))
  if (verbose) cat(sprintf("Selected %d cluster(s) by BIC; updating STlist...\n", mfit$G))
  
  for (s in samples) {
    drop_cols <- intersect(c("stniche_dsplFalse", "STniche_DTC"), colnames(x@spatial_meta[[s]]))
    if (length(drop_cols)) {
      x@spatial_meta[[s]][, drop_cols] <- NULL
    }
    add <- cl_df[cl_df$samplename == s, c("libname", "stniche_dsplFalse", "STniche_DTC"), drop = FALSE]
    x@spatial_meta[[s]] <- dplyr::full_join(x@spatial_meta[[s]], add, by = "libname")
  }
  
  if (verbose) cat("Done.\n")
  if (out_i) return(list(x, Zi)) else return(x)
}


# -----------------------------------------------------------------------------
# Embedded helper: calculate_gs_gsva_score
# -----------------------------------------------------------------------------
calculate_gs_gsva_score <- function(x = NULL, combo = NULL, gene_sets = NULL,
                                    pw_genes = NULL, min_genes = 1,
                                    cores = 1L, verbose = FALSE) {
  
  # Respect check-time core limit inside the helper as well
  if (identical(Sys.getenv("_R_CHECK_LIMIT_CORES_"), "TRUE")) {
    cores <- min(as.integer(cores), 2L)
  }
  
  samp_levels <- unique(combo$sample)
  res <- lapply(seq_along(samp_levels), function(ii) {
    s <- samp_levels[ii]
    pw_idx <- which(combo$sample == s)
    pw_map <- pw_genes[pw_idx]
    names(pw_map) <- combo$gene_set[pw_idx]
    
    keep_sets <- names(pw_map)[vapply(pw_map, length, integer(1)) >= min_genes]
    gs_for_s  <- gene_sets[keep_sets]
    
    if (length(gs_for_s) > 0) {
      # minSize = 2: a gene set of size 1 produces a degenerate GSVA score
      # (constant or near-constant) which then yields all-NA Local Moran's Z,
      # which in turn breaks mclust. Filtering at size 2 before GSVA
      # guarantees every surviving gene set carries real spatial signal.
      gp <- GSVA::gsvaParam(as.matrix(x[[s]]), geneSets = gs_for_s,
                            minSize = 2)
      sc <- GSVA::gsva(gp,
                       BPPARAM = BiocParallel::MulticoreParam(workers = cores),
                       verbose = verbose)
      sc <- as.data.frame(sc)
    } else {
      sc <- data.frame(matrix(nrow = length(pw_map), ncol = ncol(x[[s]])))
      rownames(sc) <- names(pw_map)
      colnames(sc) <- colnames(x[[s]])
    }
    
    # Ensure all gene_sets are present (fill missing with NA)
    if (length(gs_for_s) != length(gene_sets)) {
      sc_full <- lapply(names(gene_sets), function(gs) {
        if (gs %in% rownames(sc)) {
          sc[gs, , drop = FALSE]
        } else {
          df <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = ncol(sc)))
          rownames(df) <- gs
          colnames(df) <- colnames(sc)
          df
        }
      })
      sc <- do.call(rbind, sc_full)
    }
    sc
  })
  names(res) <- samp_levels
  return(res)
}