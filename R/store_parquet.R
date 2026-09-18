# Optional on-disk backing for an STlist's count matrices, as parquet in COO
# (long) form.
#
# OPT-IN BY DESIGN. Nothing here runs unless a user asks for it. An STlist built
# the normal way keeps dgCMatrix objects in @counts/@tr_counts exactly as before;
# STlist_to_disk(x, root) swaps those elements for manifests pointing at parquet
# files, and STlist_to_memory(x) reverses it. Both representations satisfy the same
# slot contract, because every read goes through the accessors in
# R/utils_accessors.R.
#
# WHY COO RATHER THAN A MATRIX FORMAT
#
# One parquet file per sample per slot, schema (gene_idx int32, cell_idx int32,
# value double), sorted by (gene_idx, cell_idx).
#
# Sorting by gene gives parquet row-group min/max statistics on gene_idx, so
# filter(gene_idx %in% ...) prunes whole row groups without reading them -- which
# is the dominant access pattern in this package (STdiff, STgradient, STenrich and
# the gene plots all want ~100 genes by all spots). Physically hive-partitioning by
# gene would give finer pruning but create tens of thousands of tiny files, which
# arrow handles badly.
#
# Including cell_idx in the sort key is not cosmetic. It makes a per-gene
# group_by/sum visit values in the same order Matrix::rowMeans walks a row, which
# is what allows an arrow-side aggregate to agree with the in-memory one
# bit-for-bit. Sorting by gene alone leaves intra-gene order unspecified.
#
# THE MANIFEST IS THE ONLY AUTHORITY for dimensions, dimnames and order. COO omits
# structural zeros, so an all-zero gene or cell has no rows at all -- and
# filter_data()'s defaults are all >= 0, so such genes and cells are retained and
# real. Reconstructing dimnames from the data would silently shrink and reorder the
# matrix, which fidelity rule 13(c) makes a genuine regression.
#
# WHAT IS NOT DONE HERE. The aggregate call sites still materialize and compute in
# memory; pushing colSums/rowMeans-style work down into arrow queries is follow-on
# work. So this buys capacity -- an STlist larger than RAM can be held and read
# from in gene-subset slices -- not yet a speedup on the aggregate paths.


##
# @title sge_manifest
# @description Describes one sample's on-disk matrix. Held in the slot element in
# place of a dgCMatrix, so sample_counts() is the only seam that needs to know.
# @return an object of class 'sge_manifest'.
#
new_sge_manifest = function(sample, store_id, path, generation, genes, cells,
                            nnz, value_type, explicit_zeros, fingerprint){
  structure(list(
    sample = sample,
    store_id = store_id,
    path = path,
    generation = generation,
    genes = genes,
    cells = cells,
    n_genes = length(genes),
    n_cells = length(cells),
    nnz = nnz,
    gene_dup = anyDuplicated(genes) > 0,
    cell_dup = anyDuplicated(cells) > 0,
    value_type = value_type,
    sort_key = c('gene_idx', 'cell_idx'),
    explicit_zeros = explicit_zeros,
    fingerprint = fingerprint
  ), class = 'sge_manifest')
}

is_sge_manifest = function(x) inherits(x, 'sge_manifest')

#' @export
print.sge_manifest = function(x, ...){
  cat(sprintf('<sge_manifest> %s: %d genes x %d cells, %.0f nonzero, gen %d\n',
              x$sample, x$n_genes, x$n_cells, x$nnz, x$generation))
  cat(sprintf('  %s\n', x$path))
  invisible(x)
}

##
# Dimension/dimname accessors used by R/utils_accessors.R. base::nrow() and
# friends return NULL on a manifest, so the accessors dispatch here instead.
#
manifest_dim = function(m) c(m$n_genes, m$n_cells)
manifest_genes = function(m) m$genes
manifest_cells = function(m) m$cells


##
# @title sge_store_id
# @description A filesystem-safe directory segment for a sample.
#
# Sample names are NOT usable as paths. Fidelity rule 11 already documents that
# real names contain ", $, backticks and backslashes; paths add /, :, leading dots,
# macOS NFC/NFD normalization, the 255-byte component limit, and case-insensitive
# collision on APFS (Sample1 and sample1 are distinct samples but one file). The
# index keeps ordering stable and the hash keeps it unique; the human-readable name
# lives only in the manifest.
#
sge_store_id = function(sample, index){
  h = substr(digest_string(sample), 1, 12)
  sprintf('%04d-%s', index, h)
}

## Small dependency-free string hash, so this does not pull in digest().
digest_string = function(s){
  bytes = utf8ToInt(enc2utf8(s))
  if(length(bytes) == 0) bytes = 0L
  h1 = 2166136261
  for(b in bytes){
    h1 = bitwAnd(as.integer((bitwXor(as.integer(h1 %% 2^31), b) * 16777619) %% 2^31), .Machine$integer.max)
  }
  h2 = sum(bytes * seq_along(bytes)) %% 2^31
  sprintf('%08x%08x', as.integer(h1), as.integer(h2))
}


##
# @title sge_store_write_matrix
# @description Writes one matrix to parquet in COO form and returns its manifest.
#
# Ordering comes for free from a transpose rather than an explicit sort: a
# dgCMatrix is CSC, so t(m) has the original genes as its columns, and walking its
# @p/@i gives (gene ascending, cell ascending) directly. Sorting 10^8 rows would
# otherwise dominate the write.
#
sge_store_write_matrix = function(m, path, sample, store_id, generation = 1L){
  if(!inherits(m, 'dgCMatrix')){
    m = methods::as(methods::as(methods::as(m, 'dMatrix'), 'generalMatrix'), 'CsparseMatrix')
  }
  genes = rownames(m)
  cells = colnames(m)
  if(is.null(genes) || is.null(cells)){
    raise_err(err_code = 'error0059')
  }

  explicit_zeros = any(m@x == 0)

  tm = Matrix::t(m)   # cells x genes, CSC -> columns are genes
  coo = tibble::tibble(
    gene_idx = rep.int(seq_len(ncol(tm)), diff(tm@p)),
    cell_idx = tm@i + 1L,
    value = tm@x
  )
  coo$gene_idx = as.integer(coo$gene_idx)
  coo$cell_idx = as.integer(coo$cell_idx)

  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  # Write to a temporary file then rename. A reader must never observe a partially
  # written parquet file, and rename within a directory is atomic.
  tmp = paste0(path, '.tmp')
  arrow::write_parquet(coo, tmp)
  file.rename(tmp, path)

  new_sge_manifest(
    sample = sample, store_id = store_id, path = normalizePath(path),
    generation = generation, genes = genes, cells = cells,
    nnz = nrow(coo), value_type = 'double',
    explicit_zeros = explicit_zeros,
    fingerprint = sge_fingerprint(path))
}

sge_fingerprint = function(path){
  info = file.info(path)
  sprintf('%.0f-%s', info$size, format(info$mtime, '%Y%m%d%H%M%S'))
}


##
# @title sge_store_materialize
# @description Reconstructs a dgCMatrix from a manifest, optionally for a subset
# of genes.
#
# The match() against sorted indices is fidelity rule 9 arriving from the opposite
# direction. collect() returns rows in row-group order, which under a gene sort is
# ascending gene_idx -- so rows must be mapped back through the MANIFEST's gene
# order, never through the order the caller happened to request or the order the
# file returned.
#
sge_store_materialize = function(m, genes = NULL, cells = NULL){
  stopifnot(is_sge_manifest(m))
  if(!file.exists(m$path)){
    raise_err(err_code = 'error0060')
  }

  keep_g = if(is.null(genes)) seq_len(m$n_genes) else sort(match(intersect(genes, m$genes), m$genes))
  if(length(keep_g) == 0){
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(0L, m$n_cells),
                                dimnames = list(character(0), m$cells)))
  }

  q = arrow::read_parquet(m$path, as_data_frame = FALSE)
  tb = as.data.frame(
    dplyr::collect(dplyr::filter(q, gene_idx %in% keep_g)))

  out = Matrix::sparseMatrix(
    i = match(tb$gene_idx, keep_g),
    j = tb$cell_idx,
    x = tb$value,
    dims = c(length(keep_g), m$n_cells),
    dimnames = list(m$genes[keep_g], m$cells))

  if(!is.null(cells)){
    keep_c = intersect(cells, m$cells)
    out = out[, keep_c, drop = FALSE]
  }
  out
}


##
#' @title STlist_to_disk: Move an STlist's count matrices to on-disk parquet
#' @description Converts an in-memory STlist into one backed by parquet files, so
#' data larger than memory can be held and read in slices.
#' @details
#' Opt-in and reversible. `@counts` and `@tr_counts` elements are replaced by
#' manifests describing per-sample parquet files; every other slot is untouched.
#' All package functions continue to work, because slot reads go through internal
#' accessors that materialize on demand -- so this trades memory for I/O rather
#' than changing results.
#'
#' Matrices are stored in COO (long) form, one file per sample per slot, sorted by
#' gene then cell. `STlist_to_memory()` reverses the conversion.
#'
#' Note the store is a directory of files, not part of the R object. An STlist
#' saved with `saveRDS()` records paths, so moving or deleting `root` breaks it;
#' `validate_store()` reports that rather than failing obscurely later.
#'
#' @param x an STlist
#' @param root directory to write the store into. Created if absent
#' @param overwrite logical, whether to overwrite an existing store at `root`
#' @return an STlist whose count slots are on-disk backed
#'
#' @export
#
STlist_to_disk = function(x = NULL, root = NULL, overwrite = FALSE){
  assert_stlist(x)
  if(is.null(root) || !is.character(root) || length(root) != 1){
    raise_err(err_code = 'error0061')
  }
  if(dir.exists(file.path(root, 'samples')) && !overwrite){
    raise_err(err_code = 'error0062')
  }
  dir.create(root, showWarnings = FALSE, recursive = TRUE)

  # arrow runs its own CPU thread pool sized from detectCores(). R CMD check sets
  # _R_CHECK_LIMIT_CORES_ to catch exactly that, and count_cores() already clamps
  # for it, so arrow is clamped the same way rather than being a second, unguarded
  # source of threads (fidelity rule 12).
  sge_clamp_arrow_threads()

  samples = names(x@counts)
  for(k in seq_along(samples)){
    s = samples[k]
    sid = sge_store_id(s, k)
    x@counts[[s]] = sge_store_write_matrix(
      sample_counts(x, s), file.path(root, 'samples', sid, 'counts.parquet'), s, sid)
    if(!rlang::is_empty(x@tr_counts) && !is.null(x@tr_counts[[s]])){
      x@tr_counts[[s]] = sge_store_write_matrix(
        sample_tr_counts(x, s), file.path(root, 'samples', sid, 'tr_counts.parquet'), s, sid)
    }
  }

  # sample_order is stored explicitly so fidelity rule 7 is assertable rather than
  # inferred. portable=FALSE marks a store under tempdir(), which will not survive
  # the session -- so saveRDS() of the STlist would produce a broken object.
  x@misc[['store']] = list(
    root = normalizePath(root),
    sample_order = samples,
    schema_ver = 1L,
    portable = !startsWith(normalizePath(root), normalizePath(tempdir()))
  )
  x
}


##
#' @title STlist_to_memory: Bring an on-disk STlist back into memory
#' @description Reverses `STlist_to_disk()`, materializing every backed matrix.
#' @param x an STlist, backed or not
#' @return an STlist with in-memory count matrices
#'
#' @export
#
STlist_to_memory = function(x = NULL){
  assert_stlist(x)
  for(s in names(x@counts)){
    if(is_sge_manifest(x@counts[[s]])) x@counts[[s]] = sge_store_materialize(x@counts[[s]])
    if(!rlang::is_empty(x@tr_counts) && is_sge_manifest(x@tr_counts[[s]])){
      x@tr_counts[[s]] = sge_store_materialize(x@tr_counts[[s]])
    }
  }
  x@misc[['store']] = NULL
  x
}


##
#' @title is_backed: Is an STlist's expression data stored on disk?
#' @param x an STlist
#' @return logical
#' @export
#
is_backed = function(x = NULL){
  assert_stlist(x)
  length(x@counts) > 0 && is_sge_manifest(x@counts[[1]])
}


##
#' @title validate_store: Check an on-disk STlist against its manifests
#' @description Verifies that a backed STlist's files exist and still match what
#' its manifests record, reading parquet footers only.
#' @details
#' Checks sample order against the recorded `sample_order` (fidelity rule 7), that
#' every file exists and its size/mtime fingerprint is unchanged, that dimnames
#' lengths agree with the recorded dimensions, and that `@tr_counts` genes and
#' cells are subsets of `@counts`' in the same relative order -- the invariant
#' `filter_data()` assumes. With `deep=TRUE`, also reads the data to check index
#' bounds and duplicate coordinates.
#'
#' @param x a backed STlist
#' @param deep logical, whether to read data pages as well as metadata
#' @return invisibly TRUE, or stops with the first problem found
#'
#' @export
#
validate_store = function(x = NULL, deep = FALSE){
  assert_stlist(x)
  if(!is_backed(x)){
    raise_err(err_code = 'error0063')
  }
  store = x@misc[['store']]

  if(!identical(names(x@counts), store$sample_order)){
    raise_err(err_code = 'error0064')
  }

  for(s in names(x@counts)){
    for(slot in c('counts', 'tr_counts')){
      m = if(slot == 'counts') x@counts[[s]] else x@tr_counts[[s]]
      if(is.null(m) || !is_sge_manifest(m)) next

      if(!file.exists(m$path)) raise_err(err_code = 'error0060')
      if(!identical(sge_fingerprint(m$path), m$fingerprint)) raise_err(err_code = 'error0065')
      if(length(m$genes) != m$n_genes || length(m$cells) != m$n_cells){
        raise_err(err_code = 'error0066')
      }
      if(m$gene_dup) raise_err(err_code = 'error0067')

      if(deep){
        tb <- as.data.frame(arrow::read_parquet(m$path))
        if(nrow(tb) > 0){
          if(max(tb$gene_idx) > m$n_genes || min(tb$gene_idx) < 1L ||
             max(tb$cell_idx) > m$n_cells || min(tb$cell_idx) < 1L){
            raise_err(err_code = 'error0068')
          }
          if(anyDuplicated(paste(tb$gene_idx, tb$cell_idx))) raise_err(err_code = 'error0068')
        }
      }
    }

    # Cross-slot: tr_counts must be a same-order subset of counts.
    if(!rlang::is_empty(x@tr_counts) && is_sge_manifest(x@tr_counts[[s]])){
      cg = x@counts[[s]]$genes; tg = x@tr_counts[[s]]$genes
      if(!identical(tg, cg[cg %in% tg])) raise_err(err_code = 'error0069')
      cc = x@counts[[s]]$cells; tc = x@tr_counts[[s]]$cells
      if(!identical(tc, cc[cc %in% tc])) raise_err(err_code = 'error0069')
    }
  }
  invisible(TRUE)
}


##
# @title sge_clamp_arrow_threads
# @description Keeps arrow's own CPU thread pool inside the same limit
# count_cores() respects. arrow sizes its pool from detectCores() independently of
# anything in this package, so under _R_CHECK_LIMIT_CORES_ it is a second,
# unguarded source of threads that R CMD check would flag (fidelity rule 12). Also
# honors the spatialGE.cores override, so forcing serial really is serial.
#
sge_clamp_arrow_threads = function(){
  limit = Sys.getenv('_R_CHECK_LIMIT_CORES_', unset = '')
  override = sge_cores_override()
  target = NULL
  if(nzchar(limit) && !identical(tolower(limit), 'false')) target = 2L
  if(!is.null(override)) target = min(c(target, max(1L, override)))
  if(!is.null(target)){
    try(arrow::set_cpu_count(target), silent = TRUE)
  }
  invisible(NULL)
}
