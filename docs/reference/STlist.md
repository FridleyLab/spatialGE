# STlist: Creation of STlist objects for spatial transcriptomics analysis

Creates an STlist object from one or multiple spatial transcriptomic
samples.

## Usage

``` r
STlist(
  rnacounts = NULL,
  spotcoords = NULL,
  samples = NULL,
  cores = NULL,
  verbose = TRUE
)
```

## Arguments

- rnacounts:

  the count data which can be provided in one of these formats:

  - File paths to comma- or tab-delimited files containing raw gene
    counts, one file for each spatial sample. The first column contains
    gene names and subsequent columns contain the expression data for
    each cell/spot. Duplicate gene names will be modified using
    `make.unique`. Requires `spotcoords` and `samples`.

  - File paths to Visium output directories (one per spatial sample).
    The directory must follow the structure resulting from
    `spaceranger count`. The directory contains the `.h5` and `spatial`
    sub-directory. If no `.h5` file is available, sparse matrices (MEX)
    from `spaceranger count`. In that case a second sub-directory called
    `filtered_feature_bc_matrix` should contain contain the
    `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` files. The
    `spatial` sub-directory minimally contains the coordinates
    (`tissue_positions_list.csv`), and optionally the high resolution
    PNG image and accompanying scaling factors
    (`scalefactors_json.json`). Requires `samples`. \#'

  - File paths to Xenium output directories (one per spatial sample).
    The directory must follow the structure resulting from the
    `xeniumranger` pipeline. The directory contains the `.h5` or sparse
    matrices (MEX). In that case a second sub-directory called
    `cell_feature_matrix` should contain contain the `barcodes.tsv.gz`,
    `features.tsv.gz`, and `matrix.mtx.gz` files. The coordinates must
    be available in the `cells.parquet`. Requires `samples`.

  - The `exprMat` file for each slide of a CosMx-SMI output. The file
    must contain the "fov" and "cell_ID" columns. The `STlist` function
    will separate data from each FOV, since analysis in spatialGE is
    conducted at the FOV level. Requires `samples` and `spotcoords`.

  - Seurat object (V4). A Seurat V4 object produced via
    [`Seurat::Load10X_Spatial`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).
    Multiple samples are allowed as long as they are stored as "slices"
    in the Seurat object. Does not require `samples` as sample names are
    taken from `names(seurat_obj@images)`

  - A named list of data frames with raw gene counts (one data frame per
    spatial sample). Requires `spotcoords`. Argument `samples` only
    needed when a file path to sample metadata is provided.

- spotcoords:

  the cell/spot coordinates. Not required if inputs are Visium or Xenium
  (spaceranger or xeniumranger outputs).

  - File paths to comma- or tab-delimited files containing cell/spot
    coordinates, one for each spatial sample. The files must contain
    three columns: cell/spot IDs, Y positions, and X positions. The
    cell/spot IDs must match the column names for each cells/spots
    (columns) in the gene count files. Requires `samples` and
    `rnacounts`.

  - The `metadata` file for each slide of a CosMx-SMI output. The file
    must contain the "fov", "cell_ID", "CenterX_local_px", and
    "CenterY_local_px" columns. The `STlist` function will separate data
    from each FOV, since analysis in spatialGE is conducted at the FOV
    level. Requires `samples` and `rnacounts`.

  - A named list of data frames with cell/spot coordinates. The list
    names must match list names of the gene counts list

- samples:

  the sample names/IDs and (optionally) metadata associated with each
  spatial sample. The following options are available for `samples`:

  - A vector with sample names, which will be used to match gene the
    counts and cell/spot coordinates file paths. A sample name must not
    match file paths for two different samples. For example, instead of
    using "tissue1" and "tissue12", use "tissue01" and "tissue12".

  - A path to a file containing a table with metadata. This file is a
    comma- or tab-separated table with one sample per row and sample
    names/IDs in the first column. Subsequent columns may contain
    variables associated with each spatial sample

- cores:

  integer indicating the number of cores to use during parallelization.
  If NULL, the function uses half of the available cores at a maximum.
  The parallelization uses
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) and
  works only in Unix systems

- verbose:

  logical, whether to print text to console

## Value

an STlist object containing the counts and coordinates, and optionally
the sample metadata, which can be used for downstream analysis with
`spatialGE`

## Details

Objects of the S4 class STlist are the starting point of analyses in
**`spatialGE`**. The STlist contains data from one or multiple samples
(i.e., tissue slices), and results from most `spatialGE`'s functions are
stored within the object.

- Raw gene counts and spatial coordinates. Gene count data have genes in
  rows and sampling units (e.g., cells, spots) in columns. Spatial
  coordinates have sampling units in rows and three columns: sample unit
  IDs, Y position, and X position.

- Visium outputs from *Space Ranger*. The Visium directory must have the
  directory structure resulting from `spaceranger count`, with either a
  count matrix represented in MEX files or a h5 file. The directory
  should also contain a `spatial` sub-directory, with the spatial
  coordinates (`tissue_positions_list.csv`), and optionally the high
  resolution tissue image and scaling factor file
  `scalefactors_json.json`.

- Xenium outputs from *Xenium Ranger*. The Xenium directory must have
  the directory structure resulting from the `xeniumranger` pipeline,
  with either a cell-feature matrix represented in MEX files or a h5
  file. The directory should also contain a parquet file, with the
  spatial coordinates (`cells.parquet`).

- CosMx-SMI outputs. Two files are required to process SMI outputs: The
  `exprMat` and `metadata` files. Both files must contain the "fov" and
  "cell_ID" columns. In addition, the `metadata` files must contain the
  "CenterX_local_px" and "CenterY_local_px" columns.

- Seurat object (V4). A Seurat V4 object produced via
  [`Seurat::Load10X_Spatial`](https://satijalab.org/seurat/reference/Load10X_Spatial.html).

Optionally, the user can input a path to a file containing a table of
sample-level metadata (e.g., clinical outcomes, tissue type, age). This
sample metadata file should contain sample IDs in the first column
partially matching the file names of the count/coordinate file paths or
Visium directories. *Note:* The sample ID of a given sample cannot be a
substring of the sample ID of another sample. For example, instead of
using "tissue1" and "tissue12", use "tissue01" and "tissue12".

The function uses parallelization if run in a Unix system. Windows users
will experience longer times depending on the number of samples.

## Examples

``` r
# \donttest{
# Using included melanoma example (Thrane et al.)
# Download example data set from spatialGE_Data
thrane_tmp = tempdir()
unlink(thrane_tmp, recursive=TRUE)
dir.create(thrane_tmp)
lk='https://github.com/FridleyLab/spatialGE_Data/raw/refs/heads/main/melanoma_thrane.zip?download='
tryCatch({ # In case data is not available from network
  download.file(lk, destfile=paste0(thrane_tmp, '/', 'melanoma_thrane.zip'), mode='wb')
  zip_tmp = list.files(thrane_tmp, pattern='melanoma_thrane.zip$', full.names=TRUE)
  unzip(zipfile=zip_tmp, exdir=thrane_tmp)
  # Generate the file paths to be passed to the STlist function
  count_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='counts')
  coord_files <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                            full.names=TRUE, pattern='mapping')
  clin_file <- list.files(paste0(thrane_tmp, '/melanoma_thrane'),
                          full.names=TRUE, pattern='clinical')
  # Create STlist
  library('spatialGE')
  melanoma <- STlist(rnacounts=count_files,
                     spotcoords=coord_files,
                     samples=clin_file)
  melanoma
}, error = function(e) {
  message("Could not run example. Are you connected to the internet?")
  return(NULL)
})
#> Warning: URL 'https://media.githubusercontent.com/media/FridleyLab/spatialGE_Data/refs/heads/main/melanoma_thrane.zip?download=true': status was 'SSL connect error'
#> Could not run example. Are you connected to the internet?
#> NULL
# }
```
