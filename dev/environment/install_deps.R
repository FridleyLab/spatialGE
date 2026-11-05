## Install every Imports / LinkingTo / Suggests dependency of spatialGE
## into the currently active R (expected: the spatialge-dev conda env).
## Needs network access. Run with:
##   conda run -n spatialge-dev Rscript dev/environment/install_deps.R
## or, from an interactive shell with that env activated:
##   Rscript dev/environment/install_deps.R

repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos)

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

## Version pins, as pkg = "version". Currently empty.
##
## ggplot2 used to be pinned to 3.5.1 here, on the theory that 4.x's S7 rewrite
## would change the ggplot_build() layer data the golden plot comparisons diff.
## That pin is gone because it was both stale and actively harmful: the golden
## baseline is captured under ggplot2 4.x, so 4.x is the correct comparison
## basis, and ggrepel (>= 0.9.8) and ggpubr (>= 1.0.0) both declare
## ggplot2 (>= 3.5.2), so installing 3.5.1 breaks package loading outright.
## Don't reinstate a ggplot2 pin without re-capturing tests/golden/ first.
pinned <- list()

cran_pkgs <- c(
  "arrow", "concaveman", "data.table", "dplyr", "ggforce", "ggpolypath",
  "ggrepel", "gstat", "hdf5r", "jpeg", "jsonlite", "khroma", "magrittr",
  "Matrix", "MASS", "png", "RColorBrewer", "Rcpp", "readr", "readxl",
  "rlang", "scales", "sctransform", "sfsmisc", "sf", "sp", "spaMM",
  "spdep", "stringr", "tibble", "tidyr", "uwot", "wordspace",
  "RcppEigen", "RcppProgress",
  "curl", "geoR", "ggpubr", "httr", "janitor", "kableExtra", "knitr",
  "msigdbr", "progress", "rmarkdown", "rpart", "R.utils", "spatstat",
  "tidyverse", "testthat", "bench", "profvis", "covr", "waldo", "devtools",
  "roxygen2", "writexl"
)

bioc_pkgs <- c(
  "BiocParallel", "ComplexHeatmap", "DelayedArray", "DelayedMatrixStats",
  "dynamicTreeCut", "EBImage", "GSVA"
)

install_pinned_or_latest <- function(pkgs) {
  for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) next
    if (!is.null(pinned[[pkg]])) {
      remotes::install_version(pkg, version = pinned[[pkg]], upgrade = "never")
    } else {
      install.packages(pkg)
    }
  }
}

install_pinned_or_latest(cran_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)

## scSpatialSIM and SeuratObject are not on CRAN/Bioconductor.
if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  remotes::install_github("mojaveazure/seurat-object", upgrade = "never")
}
if (!requireNamespace("scSpatialSIM", quietly = TRUE)) {
  remotes::install_github("FridleyLab/scSpatialSIM", upgrade = "never")
}

missing_after <- Filter(
  function(p) !requireNamespace(p, quietly = TRUE),
  c(cran_pkgs, bioc_pkgs, "SeuratObject", "scSpatialSIM")
)
if (length(missing_after) > 0) {
  message("Still missing after install: ", paste(missing_after, collapse = ", "))
}

dir.create("dev/environment", showWarnings = FALSE, recursive = TRUE)
sink("dev/environment/RECORDED.md")
tryCatch({
  cat("# Recorded build environment\n\n")
  cat("## sessionInfo()\n\n```\n")
  print(sessionInfo())
  cat("```\n\n## extSoftVersion()\n\n```\n")
  print(extSoftVersion())
  cat("```\n\n## BLAS / LAPACK\n\n```\n")
  cat(La_library(), "\n")
  cat(extSoftVersion()["BLAS"], "\n")
  cat("```\n")
}, finally = sink())

message("Done. See dev/environment/RECORDED.md for the recorded environment.")
