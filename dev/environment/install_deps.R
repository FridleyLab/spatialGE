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

## ggplot2 is pinned below 4.0: the 4.x S7 rewrite changes ggplot_build()
## layer-data columns that the golden-baseline plot comparisons (see
## dev/baseline/capture_baseline.R) depend on, and STdiff_volcano.R's
## element_rect(fill = NULL) call is a known 4.x compatibility risk.
pinned <- list(ggplot2 = "3.5.1")

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
