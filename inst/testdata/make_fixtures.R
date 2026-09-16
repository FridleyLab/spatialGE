## Synthetic, seeded, offline fixtures for spatialGE's test suite.
##
## This file only defines generator functions -- it has no side effects at
## source time. tests/testthat/helper-fixtures.R sources this and materializes
## the fixtures into tempdir() once per test session; nothing here is ever
## committed as a binary. See plans/this-is-an-r-velvety-music.md Phase 2 for
## the rationale behind the specific degenerate cases below.

FIXTURE_SEED <- 20260115

## 3 samples (not 1 or 2): `mclapply` short-circuits to `lapply` when
## length(X) < 2, so fewer than 3 samples never exercises real chunking.
## Names avoid leading digits, use only the STlist()-allowed character class,
## and no name is a substring of another.
fixture_sample_names <- function() c("alpha", "bravo", "charlie")

## ~200 genes, comfortably under sctransform's sct_n_regr_genes (3000) so
## `sctransform::vst` skips its internal gene-sampling step and stays
## reproducible across mc.cores. Includes deliberately degenerate/adversarial
## names: an all-zero gene, a constant gene, one with a shell metacharacter,
## and one with a syntactically-invalid R name (leading digit) -- together
## these prove `check.names = FALSE` is preserved end to end.
fixture_gene_names <- function() {
  c(
    paste0("Gene", seq_len(184)),
    paste0("GradGene", seq_len(10)), # spatial-gradient markers
    "DomainGeneA", "DomainGeneB",    # 2-domain markers
    "ZeroGene",
    "ConstGene",
    "Gene$special",
    "1BadGene"
  )
}

## 12 x 10 spatial grid (120 spots), split into two domains along x.
fixture_spot_grid <- function() {
  grid <- expand.grid(xpos = seq_len(12), ypos = seq_len(10))
  grid$domain <- ifelse(grid$xpos <= 6, "A", "B")
  grid
}

fixture_counts_matrix <- function(genes, grid) {
  n_genes <- length(genes)
  n_spots <- nrow(grid)
  mat <- matrix(0L, nrow = n_genes, ncol = n_spots, dimnames = list(genes, NULL))
  for (gi in seq_len(n_genes)) {
    gene <- genes[gi]
    if (gene == "ZeroGene") {
      next # stays all-zero
    } else if (gene == "ConstGene") {
      mat[gi, ] <- 7L # constant -> sd == 0, hits the cor.test degenerate branch
    } else if (grepl("^GradGene", gene)) {
      lambda <- 1 + 0.5 * grid$xpos
      mat[gi, ] <- rpois(n_spots, lambda)
    } else if (gene == "DomainGeneA") {
      lambda <- ifelse(grid$domain == "A", 14, 1)
      mat[gi, ] <- rpois(n_spots, lambda)
    } else if (gene == "DomainGeneB") {
      lambda <- ifelse(grid$domain == "B", 14, 1)
      mat[gi, ] <- rpois(n_spots, lambda)
    } else {
      mat[gi, ] <- rpois(n_spots, 4)
    }
  }
  mat
}

## The shared synthetic dataset every format-specific writer below serializes
## from. Baking degenerate content in here (rather than per-writer) guarantees
## every input layout exercises the same edge cases.
fixture_base_data <- function() {
  set.seed(FIXTURE_SEED)
  samples <- fixture_sample_names()
  genes <- fixture_gene_names()
  grid <- fixture_spot_grid()
  n_spots <- nrow(grid)

  counts_list <- list()
  coords_list <- list()

  for (s in samples) {
    barcodes <- paste0(s, "_spot", sprintf("%03d", seq_len(n_spots)))
    mat <- fixture_counts_matrix(genes, grid)
    if (s == "bravo") mat[, 1] <- 0L # zero-library-size spot
    colnames(mat) <- barcodes

    counts_df <- data.frame(gene = genes, mat, check.names = FALSE, stringsAsFactors = FALSE)
    coords_df <- data.frame(libname = barcodes, ypos = grid$ypos, xpos = grid$xpos, stringsAsFactors = FALSE)
    if (s == "charlie") coords_df <- coords_df[sample.int(n_spots), ] # shuffled vs. counts column order

    counts_list[[s]] <- counts_df
    coords_list[[s]] <- coords_df
  }

  list(samples = samples, genes = genes, grid = grid, counts = counts_list, coords = coords_list)
}

## ---- format writers -------------------------------------------------------
## Each returns a list with (a subset of) rnacounts/spotcoords/samples ready
## to hand straight to STlist().

fixture_write_list_dfs <- function(base) {
  list(rnacounts = base$counts, spotcoords = base$coords, samples = NULL)
}

fixture_write_delim <- function(base, dir, delim = c(",", "\t")) {
  delim <- match.arg(delim)
  ext <- if (delim == ",") "csv" else "tsv"
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  rna_paths <- character(0)
  coord_paths <- character(0)
  for (s in base$samples) {
    rna_fp <- file.path(dir, paste0(s, "_counts.", ext))
    coord_fp <- file.path(dir, paste0(s, "_coords.", ext))
    utils::write.table(base$counts[[s]], rna_fp, sep = delim, row.names = FALSE, quote = FALSE)
    utils::write.table(base$coords[[s]], coord_fp, sep = delim, row.names = FALSE, quote = FALSE)
    rna_paths[s] <- rna_fp
    coord_paths[s] <- coord_fp
  }
  list(rnacounts = rna_paths, spotcoords = coord_paths, samples = base$samples)
}

## samples file alone (csv or xlsx) alongside list-of-data.frames counts/coords
## -- exercises detect_input_samples()'s csv-vs-xlsx branches specifically.
fixture_write_samplesfile <- function(base, dir, format = c("csv", "xlsx")) {
  format <- match.arg(format)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ## detect_input_samples() sniffs line 2 for a literal comma/tab, so a single-
  ## column file is never classified as delimited -- add a harmless 2nd column.
  samples_df <- data.frame(sample = base$samples, group = "all", stringsAsFactors = FALSE)
  fp <- file.path(dir, paste0("samples.", format))
  if (format == "csv") {
    utils::write.csv(samples_df, fp, row.names = FALSE, quote = FALSE)
  } else {
    writexl::write_xlsx(samples_df, fp)
  }
  list(rnacounts = base$counts, spotcoords = base$coords, samples = fp)
}

## samplesfile_matrices: a samples file whose own columns point at per-sample
## counts/coords file paths (rnacounts/spotcoords are NOT passed to STlist()).
fixture_write_samplesfile_matrices <- function(base, dir, delim = ",") {
  matdir <- file.path(dir, "matrices")
  written <- fixture_write_delim(base, matdir, delim = delim)
  samples_df <- data.frame(
    sample = base$samples,
    counts_path = unname(written$rnacounts[base$samples]),
    coords_path = unname(written$spotcoords[base$samples]),
    stringsAsFactors = FALSE
  )
  samples_fp <- file.path(dir, "samples.csv")
  utils::write.csv(samples_df, samples_fp, row.names = FALSE, quote = FALSE)
  list(rnacounts = NULL, spotcoords = NULL, samples = samples_fp)
}

## MEX trio (features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz) shared by the
## Visium- and Xenium-MEX writers below.
fixture_write_mex_trio <- function(mat, genes, barcodes, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  features_fp <- file.path(dir, "features.tsv.gz")
  barcodes_fp <- file.path(dir, "barcodes.tsv.gz")
  matrix_fp <- file.path(dir, "matrix.mtx.gz")

  features_df <- data.frame(
    ensembl = paste0("ENSG", sprintf("%08d", seq_along(genes))),
    gene = genes,
    dtype = "Gene Expression"
  )
  con <- gzfile(features_fp, "wt")
  utils::write.table(features_df, con, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  close(con)

  con <- gzfile(barcodes_fp, "wt")
  writeLines(barcodes, con)
  close(con)

  ## Matrix::writeMM emits exactly the 3-line MatrixMarket header that
  ## import_visium()'s `fread(..., skip = 3)` expects, so the gz can just
  ## wrap its plain-text output line for line.
  sp_mat <- Matrix::Matrix(mat, sparse = TRUE)
  mtx_tmp <- tempfile(fileext = ".mtx")
  Matrix::writeMM(sp_mat, mtx_tmp)
  mtx_lines <- readLines(mtx_tmp)
  unlink(mtx_tmp)
  con <- gzfile(matrix_fp, "wt")
  writeLines(mtx_lines, con)
  close(con)

  invisible(list(features = features_fp, barcodes = barcodes_fp, matrix = matrix_fp))
}

## Headerless 6-column tissue_positions_list.csv:
## barcode, in_tissue, array_row, array_col, imagerow(px), imagecol(px)
fixture_write_visium_positions <- function(barcodes, grid, fp) {
  df <- data.frame(
    barcode = barcodes,
    in_tissue = 1L,
    array_row = grid$ypos,
    array_col = grid$xpos,
    imagerow = grid$ypos * 100L,
    imagecol = grid$xpos * 100L
  )
  utils::write.table(df, fp, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

fixture_write_scalefactors <- function(fp) {
  jsonlite::write_json(
    list(
      tissue_hires_scalef = 1, tissue_lowres_scalef = 1,
      fiducial_diameter_fullres = 10, spot_diameter_fullres = 10
    ),
    fp, auto_unbox = TRUE
  )
}

fixture_write_image <- function(fp, side = 20) {
  arr <- array(runif(side * side * 3), dim = c(side, side, 3))
  png::writePNG(arr, fp)
}

## Visium, MEX layout: one directory per sample, each containing a dummy
## filtered_feature_bc_matrix.tar.gz (only its filename is ever checked) plus
## the real unpacked filtered_feature_bc_matrix/ trio and a spatial/ folder.
fixture_write_visium_mex <- function(base, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in base$samples) {
    sdir <- file.path(dir, s)
    dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
    file.create(file.path(sdir, "filtered_feature_bc_matrix.tar.gz"))

    barcodes <- colnames(base$counts[[s]])[-1]
    mat <- as.matrix(base$counts[[s]][, -1, drop = FALSE])
    rownames(mat) <- base$genes
    fixture_write_mex_trio(mat, base$genes, barcodes, file.path(sdir, "filtered_feature_bc_matrix"))

    spdir <- file.path(sdir, "spatial")
    dir.create(spdir, showWarnings = FALSE, recursive = TRUE)
    fixture_write_visium_positions(barcodes, base$grid, file.path(spdir, "tissue_positions_list.csv"))
    fixture_write_scalefactors(file.path(spdir, "scalefactors_json.json"))
    fixture_write_image(file.path(spdir, "tissue_lowres_image.png"))
  }
  list(rnacounts = file.path(dir, base$samples), samples = base$samples)
}

## Xenium, MEX layout. The only fixture that drives read_xenium_outs()'s
## xenium_mex branch (import_xenium_mex()).
fixture_write_xenium_mex <- function(base, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in base$samples) {
    sdir <- file.path(dir, s)
    dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
    file.create(file.path(sdir, "cell_feature_matrix.tar.gz"))

    barcodes <- colnames(base$counts[[s]])[-1]
    mat <- as.matrix(base$counts[[s]][, -1, drop = FALSE])
    rownames(mat) <- base$genes
    fixture_write_mex_trio(mat, base$genes, barcodes, file.path(sdir, "cell_feature_matrix"))

    coords_df <- data.frame(
      cell_id = barcodes,
      x_centroid = base$grid$xpos * 10,
      y_centroid = base$grid$ypos * 10
    )
    arrow::write_parquet(coords_df, file.path(sdir, "cells.parquet"))
  }
  list(rnacounts = file.path(dir, base$samples), samples = base$samples)
}

## 10x-style HDF5 matrix (/matrix/{data,indices,indptr,shape,barcodes,
## features/name}). Not committed anywhere -- hdf5 isn't byte-reproducible --
## only ever generated fresh into tempdir() by the test helper.
fixture_write_10x_h5 <- function(base, sample, fp) {
  mat <- as.matrix(base$counts[[sample]][, -1, drop = FALSE])
  rownames(mat) <- base$genes
  sp_mat <- Matrix::Matrix(mat, sparse = TRUE)
  barcodes <- colnames(mat)

  h5 <- hdf5r::H5File$new(fp, mode = "w")
  grp <- h5$create_group("matrix")
  grp[["data"]] <- sp_mat@x
  grp[["indices"]] <- sp_mat@i
  grp[["indptr"]] <- sp_mat@p
  grp[["shape"]] <- dim(sp_mat)
  grp[["barcodes"]] <- barcodes
  feat_grp <- grp$create_group("features")
  feat_grp[["name"]] <- base$genes
  h5$close_all()
  invisible(fp)
}

fixture_write_visium_h5 <- function(base, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in base$samples) {
    sdir <- file.path(dir, s)
    dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
    fixture_write_10x_h5(base, s, file.path(sdir, "filtered_feature_bc_matrix.h5"))

    spdir <- file.path(sdir, "spatial")
    dir.create(spdir, showWarnings = FALSE, recursive = TRUE)
    barcodes <- colnames(base$counts[[s]])[-1]
    fixture_write_visium_positions(barcodes, base$grid, file.path(spdir, "tissue_positions_list.csv"))
    fixture_write_scalefactors(file.path(spdir, "scalefactors_json.json"))
    fixture_write_image(file.path(spdir, "tissue_lowres_image.png"))
  }
  list(rnacounts = file.path(dir, base$samples), samples = base$samples)
}

fixture_write_xenium_h5 <- function(base, dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in base$samples) {
    sdir <- file.path(dir, s)
    dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
    fixture_write_10x_h5(base, s, file.path(sdir, "cell_feature_matrix.h5"))

    barcodes <- colnames(base$counts[[s]])[-1]
    coords_df <- data.frame(
      cell_id = barcodes,
      x_centroid = base$grid$xpos * 10,
      y_centroid = base$grid$ypos * 10
    )
    arrow::write_parquet(coords_df, file.path(sdir, "cells.parquet"))
  }
  list(rnacounts = file.path(dir, base$samples), samples = base$samples)
}

## CosMx exprMat/metadata csv. Each of the 3 base samples becomes one FOV of
## a single slide, matching import_smi()'s `<slide>_fov_<n>` naming and
## process_meta()'s FOV-expansion, which matches metadata sample names by
## slide-name prefix.
fixture_write_cosmx <- function(base, dir, slide = "slideA") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  expr_rows <- list()
  meta_rows <- list()
  for (fov_i in seq_along(base$samples)) {
    s <- base$samples[fov_i]
    mat <- t(as.matrix(base$counts[[s]][, -1, drop = FALSE])) # cells x genes
    colnames(mat) <- base$genes
    n_cells <- nrow(mat)

    expr_rows[[s]] <- data.frame(
      fov = fov_i, cell_ID = seq_len(n_cells), cell = seq_len(n_cells),
      mat, check.names = FALSE
    )
    meta_rows[[s]] <- data.frame(
      fov = fov_i,
      cell_ID = seq_len(n_cells),
      CenterX_local_px = base$grid$xpos * 10,
      CenterY_local_px = base$grid$ypos * 10
    )
  }
  expr_all <- do.call(rbind, expr_rows)
  meta_all <- do.call(rbind, meta_rows)

  expr_fp <- file.path(dir, paste0(slide, "_exprMat_file.csv"))
  meta_fp <- file.path(dir, paste0(slide, "_metadata_file.csv"))
  utils::write.csv(expr_all, expr_fp, row.names = FALSE, quote = FALSE)
  utils::write.csv(meta_all, meta_fp, row.names = FALSE, quote = FALSE)

  samples_fp <- file.path(dir, "samples.csv")
  utils::write.csv(data.frame(sample = slide, group = "all", stringsAsFactors = FALSE), samples_fp, row.names = FALSE, quote = FALSE)

  list(rnacounts = expr_fp, spotcoords = meta_fp, samples = samples_fp)
}

## Materializes every synthetic-tier fixture variant under `root` (normally a
## fresh tempdir()) and returns their STlist()-ready argument lists. The h5
## variants are dev/network-tier per the plan, but are cheap and offline here,
## so they're included too -- callers that want the CRAN-safe subset only
## should skip `visium_h5`/`xenium_h5`/`xenium_mex` explicitly.
fixture_materialize_all <- function(root) {
  base <- fixture_base_data()
  list(
    base = base,
    list_dfs = fixture_write_list_dfs(base),
    list_dfs_samplesfile_csv = fixture_write_samplesfile(base, file.path(root, "samplesfile_csv"), "csv"),
    list_dfs_samplesfile_xlsx = fixture_write_samplesfile(base, file.path(root, "samplesfile_xlsx"), "xlsx"),
    delim_csv = fixture_write_delim(base, file.path(root, "delim_csv"), ","),
    delim_tsv = fixture_write_delim(base, file.path(root, "delim_tsv"), "\t"),
    samplesfile_matrices = fixture_write_samplesfile_matrices(base, file.path(root, "samplesfile_matrices")),
    visium_mex = fixture_write_visium_mex(base, file.path(root, "visium_mex")),
    xenium_mex = fixture_write_xenium_mex(base, file.path(root, "xenium_mex")),
    cosmx = fixture_write_cosmx(base, file.path(root, "cosmx")),
    visium_h5 = fixture_write_visium_h5(base, file.path(root, "visium_h5")),
    xenium_h5 = fixture_write_xenium_h5(base, file.path(root, "xenium_h5"))
  )
}
