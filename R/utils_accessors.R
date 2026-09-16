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
sample_counts = function(x, sample) x@counts[[sample]]
sample_spatial_meta = function(x, sample) x@spatial_meta[[sample]]
sample_tr_counts = function(x, sample) x@tr_counts[[sample]]
sample_gene_meta = function(x, sample) x@gene_meta[[sample]]
