## Larger synthetic fixtures for dev/profiling/, NOT for the test suite.
##
## Deliberately a separate file from make_fixtures.R, which
## tests/testthat/helper-fixtures.R sources on every test session. Nothing here
## is ever sourced by the suite -- these build objects in the thousands-to-tens-
## of-thousands of spots range and take seconds to minutes.
##
## WHY THIS EXISTS
##
## dev/profiling/FINDINGS.md records that the standard fixture (3 samples, 120
## spots, ~200 genes) puts almost every rprof sampling.time under 0.05s, too
## few samples to attribute cost. The consequence was not just imprecision: the
## worst scaling defect in the package was invisible. create_listw_from_dist()
## (R/utils.R:210) -- SThet()'s DEFAULT weighting -- builds a COMPLETE graph,
## every spot neighbouring every other, as a list of named vectors, with an
## O(n) string comparison inside an O(n) loop. At 120 spots that is 14k
## weights and looks fine. At 4,000 spots it is 16M weights plus 16M character
## names. At Xenium scale it cannot run at all.
##
## Asymptotics are only visible if the fixture reaches the range where they
## dominate, so these presets bracket it: `small` is where the current fixture
## sits, `visium` is where the quadratic term takes over, `xenium` is where it
## becomes infeasible. Benchmarking across all three is what turns "this looks
## quadratic" into a measured curve.

FIXTURE_SCALE_SEED <- 20260117

## Coordinates only. create_listw_from_dist(), create_listw_from_knn() and
## stgradient_dists_one_sample() consume nothing but coordinates, so the
## coordinate-only generator scales to sizes where a counts matrix of matching
## width would not fit -- which is exactly the range those three need probing
## at. Returns a plain data.frame with spatial_meta's column contract
## (libname/ypos/xpos), ready for the internals that take it directly.
fixture_scale_coords <- function(n_spots, sample = "scale1", seed = FIXTURE_SCALE_SEED){
  set.seed(seed)
  side <- ceiling(sqrt(n_spots))
  grid <- expand.grid(xpos = seq_len(side), ypos = seq_len(side))
  grid <- grid[seq_len(n_spots), , drop = FALSE]
  data.frame(
    libname = paste0(sample, "_spot", sprintf("%06d", seq_len(n_spots))),
    ypos = grid$ypos,
    xpos = grid$xpos,
    stringsAsFactors = FALSE
  )
}

## A sparse counts matrix at scale. lambda is the Poisson mean: 4 approximates
## a Visium spot (most genes detected), 0.3 approximates a Xenium cell (~74%
## structural zeros), which matters because a COO/parquet store's size and
## every sparse-aware aggregate scale with nnz, not with n_genes * n_cells.
##
## Built straight as a dgCMatrix via triplets rather than densified and
## converted -- a dense 50,000 x 400 integer intermediate is 80 MB and is the
## allocation this file exists to help remove, so generating one here would be
## self-defeating.
fixture_scale_counts <- function(n_spots, n_genes, lambda = 4, sample = "scale1",
                                 seed = FIXTURE_SCALE_SEED){
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  barcodes <- paste0(sample, "_spot", sprintf("%06d", seq_len(n_spots)))

  ## Draw the nonzero count directly rather than drawing n_genes * n_spots
  ## Poisson values and discarding the zeros.
  p_nonzero <- 1 - stats::dpois(0, lambda)
  n_nz <- stats::rbinom(1, size = as.double(n_genes) * n_spots, prob = p_nonzero)

  i <- sample.int(n_genes, n_nz, replace = TRUE)
  j <- sample.int(n_spots, n_nz, replace = TRUE)
  ## Truncated-at-zero Poisson for the stored values.
  v <- stats::qpois(stats::runif(n_nz, stats::dpois(0, lambda), 1), lambda)

  m <- Matrix::sparseMatrix(i = i, j = j, x = as.double(v),
                            dims = c(n_genes, n_spots),
                            dimnames = list(genes, barcodes))
  methods::as(m, "CsparseMatrix")
}

## The presets. n_genes is held low relative to real panels on purpose: the
## quadratic costs these fixtures probe are all in n_spots, and inflating
## n_genes only slows generation.
FIXTURE_SCALE_PRESETS <- list(
  ## Where the standard test fixture sits -- the control point, so a benchmark
  ## table shows the curve rather than a single number.
  small  = list(n_spots = 120L,   n_genes = 200L, lambda = 4.0),
  ## A Visium section. 16M pairwise distances: create_listw_from_dist()'s
  ## complete graph is catastrophic here and this is the preset that proves it.
  visium = list(n_spots = 4096L,  n_genes = 500L, lambda = 4.0),
  ## A small Xenium/CosMx run. 2.5e9 pairwise distances -- the complete-graph
  ## listw cannot be built at all, so this preset is the post-fix target rather
  ## than a before/after comparison.
  xenium = list(n_spots = 50000L, n_genes = 200L, lambda = 0.3)
)

fixture_scale_preset <- function(name = c("small", "visium", "xenium")){
  name <- match.arg(name)
  p <- FIXTURE_SCALE_PRESETS[[name]]
  list(
    name = name,
    n_spots = p$n_spots,
    n_genes = p$n_genes,
    coords = fixture_scale_coords(p$n_spots, sample = name),
    counts = fixture_scale_counts(p$n_spots, p$n_genes, lambda = p$lambda, sample = name)
  )
}

## An STlist at scale, assembled from the pieces above.
##
## Bypasses STlist() and its readers on purpose. The readers' own cost is
## measured separately (read_one_* units in bench_units.R); routing through
## them here would mean writing and re-parsing hundreds of MB of MatrixMarket
## text per benchmark run, and would put file I/O in the timing of units that
## do no I/O. The slots written are exactly the four STlist() populates for a
## coordinates+counts input, in the same order, so the object satisfies the
## slot contract in dev/ARCHITECTURE.md.
##
## Deliberately does NOT go through process_lists(): that would re-derive
## total_counts/total_genes by densifying, which is the cost under study.
fixture_scale_stlist <- function(preset = "visium", n_samples = 2L){
  presets <- lapply(seq_len(n_samples), function(k){
    p <- FIXTURE_SCALE_PRESETS[[preset]]
    nm <- paste0(preset, k)
    list(name = nm,
         coords = fixture_scale_coords(p$n_spots, sample = nm,
                                       seed = FIXTURE_SCALE_SEED + k),
         counts = fixture_scale_counts(p$n_spots, p$n_genes, lambda = p$lambda,
                                       sample = nm, seed = FIXTURE_SCALE_SEED + k))
  })
  names(presets) <- vapply(presets, function(e) e$name, character(1))

  counts <- lapply(presets, function(e) e$counts)
  spatial_meta <- lapply(presets, function(e){
    d <- e$coords
    d$total_counts <- as.double(Matrix::colSums(counts[[e$name]]))
    d$total_genes <- as.double(Matrix::colSums(counts[[e$name]] != 0))
    d
  })

  methods::new("STlist",
               counts = counts,
               spatial_meta = spatial_meta,
               gene_meta = stats::setNames(
                 lapply(counts, function(m) tibble::tibble(gene = rownames(m))),
                 names(counts)),
               sample_meta = tibble::tibble(sample_name = names(counts)),
               tr_counts = list(),
               gene_krige = list(),
               misc = list(platform = "visium"))
}
