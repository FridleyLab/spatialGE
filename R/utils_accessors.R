##
# @title Internal non-exported STlist slot accessors
# @description Uniform accessors for the sample-major STlist slots
# (`counts`, `spatial_meta`, `tr_counts`, `gene_meta`), used at the call
# sites touched while extracting Phase 7's named per-sample workers. This
# is the seam a future change to slot layout would go through -- e.g.
# making `gene_krige` sample-major like every other slot -- without having
# to touch every call site again. Not a validation layer: does not check
# that `sample` exists, matching how `x@slot[[sample]]` already behaved at
# every site these replace.
# @param x an STlist.
# @param sample a sample name or numeric index into the slot.
#
# When the optional parquet backing is in use (see R/store_parquet.R) a slot
# element is an sge_manifest rather than a dgCMatrix, so the two matrix accessors
# materialize on demand. `genes=` is forwarded to the store so a gene-subset read
# touches only the row groups it needs; in memory it is an ordinary subset. This is
# the whole reason the accessor layer exists.
sample_counts = function(x, sample, genes=NULL){
  m = x@counts[[sample]]
  if(is_sge_manifest(m)) return(sge_store_materialize(m, genes=genes))
  if(is.null(genes)) return(m)
  m[rownames(m) %in% genes, , drop=FALSE]
}

sample_tr_counts = function(x, sample, genes=NULL){
  m = x@tr_counts[[sample]]
  if(is_sge_manifest(m)) return(sge_store_materialize(m, genes=genes))
  if(is.null(genes)) return(m)
  m[rownames(m) %in% genes, , drop=FALSE]
}

# The materializing counterpart to sample_slot_mtx() below: resolves the slot by
# name and returns values, so a call site that already branches on the plotting
# functions' `data_type` ('tr'/'raw') does not branch a second time on which of
# the two accessors above to call. `genes=` is forwarded, so the subset is built
# at its final size rather than taken off a fully materialized matrix.
sample_slot_counts = function(x, sample, slot=c('counts', 'tr_counts'), genes=NULL){
  slot = match.arg(slot)
  if(slot == 'counts') sample_counts(x, sample, genes=genes) else sample_tr_counts(x, sample, genes=genes)
}

sample_spatial_meta = function(x, sample) x@spatial_meta[[sample]]
sample_gene_meta = function(x, sample) x@gene_meta[[sample]]


##
# Structural accessors: dimensions and dimnames WITHOUT touching values.
#
# These exist to separate "how big is this / what is it called" from "give me
# the numbers", because that split is what an out-of-core backing needs. Under
# the planned parquet store a slot element is a manifest recording `genes`,
# `cells` and the dims rather than a matrix, so `nrow(x@counts[[i]])` stops
# working while `names(x@counts)` keeps working -- the slot stays a named list
# in user order, one element per sample. That is why only the dim/dimname reads
# are routed here and the many `names(x@counts)` call sites are deliberately
# left alone: converting them would be churn with no backend benefit.
#
# The `slot` argument is named to match resolve_samples()'s, not the exported
# plotting functions' `data_type` (which means 'raw'/'tr' and is a different
# vocabulary).
#
# NOTE for whoever implements the store: these three cover the *reads*. A
# dozen sites in R/filter_data_stages.R combine a dimname read with a
# subset-assign (`x@counts[[i]][, mask]`), which needs a store-level filter
# operation rather than an accessor, and is not something this layer solves.
# Grep for `sample_genes(`/`sample_cells(` adjacency to find them.
#
# @param x an STlist.
# @param sample a sample name or numeric index.
# @param slot which sample-major counts slot to describe.
# @return sample_dim: integer c(n_genes, n_cells). sample_genes/sample_cells:
# a character vector, or NULL when the matrix carries no dimnames.
#
sample_slot_mtx = function(x, sample, slot=c('counts', 'tr_counts')){
  slot = match.arg(slot)
  if(slot == 'counts') x@counts[[sample]] else x@tr_counts[[sample]]
}

sample_dim = function(x, sample, slot=c('counts', 'tr_counts')){
  m = sample_slot_mtx(x, sample, slot)
  # A manifest answers this with zero I/O -- base::nrow() would return NULL on it.
  if(is_sge_manifest(m)) return(manifest_dim(m))
  c(base::nrow(m), base::ncol(m))
}

sample_n_genes = function(x, sample, slot=c('counts', 'tr_counts')){
  m = sample_slot_mtx(x, sample, slot)
  if(is_sge_manifest(m)) return(manifest_dim(m)[1])
  base::nrow(m)
}

sample_n_cells = function(x, sample, slot=c('counts', 'tr_counts')){
  m = sample_slot_mtx(x, sample, slot)
  if(is_sge_manifest(m)) return(manifest_dim(m)[2])
  base::ncol(m)
}

sample_genes = function(x, sample, slot=c('counts', 'tr_counts')){
  m = sample_slot_mtx(x, sample, slot)
  if(is_sge_manifest(m)) return(manifest_genes(m))
  base::rownames(m)
}

sample_cells = function(x, sample, slot=c('counts', 'tr_counts')){
  m = sample_slot_mtx(x, sample, slot)
  if(is_sge_manifest(m)) return(manifest_cells(m))
  base::colnames(m)
}
