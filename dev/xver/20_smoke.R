## Smoke-run the unmodified baseline pipeline against ONE installed version and
## print a per-step ok/error/warning summary.
##
## This runs before any cross-version comparison, on purpose. spatialGE has only
## ever been used here via devtools::load_all(), so `main` (1.2.2) has never
## actually executed against ggplot2 4.0.3 and the current dplyr/Matrix. If main
## fails here for a reason that is not on the 15-delta list, every downstream
## "the versions differ" verdict would be measuring dependency drift instead of
## the refactor.
##
## Usage:
##   . dev/xver/env.sh
##   XVER_VERSION=main Rscript --vanilla dev/xver/20_smoke.R
##   XVER_VERSION=alex Rscript --vanilla dev/xver/20_smoke.R

source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "common.R"))

which_ver <- Sys.getenv("XVER_VERSION", "main")
spec <- switch(which_ver,
  main = list(lib = Sys.getenv("XVER_LIB_MAIN"), version = "1.2.2"),
  alex = list(lib = Sys.getenv("XVER_LIB_ALEX"), version = "2.0.0"),
  stop("XVER_VERSION must be 'main' or 'alex', got '", which_ver, "'")
)

info <- xv_load_spatialge(spec$lib, spec$version)
cat(sprintf("== smoke: %s (spatialGE %s) ==\n", which_ver, info$version))

xv_source_driver()
xv_shim_fixture_image()

scratch <- xv_scratch(paste0("smoke-", which_ver))
setwd(scratch)

t0 <- Sys.time()
res <- bl_pipeline(cores = 1, verbose = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("\n%-24s %-6s %s\n", "step", "ok", "detail"))
cat(strrep("-", 78), "\n")
n_ok <- 0L
for(nm in names(res)){
  r <- res[[nm]]
  detail <- if(isTRUE(r$ok)){
    n_ok <- n_ok + 1L
    paste(class(r$value), collapse = "/")
  } else {
    paste("ERROR:", r$value)
  }
  cat(sprintf("%-24s %-6s %s\n", nm, if(isTRUE(r$ok)) "ok" else "FAIL",
              substr(detail, 1, 110)))
}
cat(strrep("-", 78), "\n")
cat(sprintf("%d/%d steps ok, %.1fs\n", n_ok, length(res), elapsed))

## Persist so the two versions' smoke results can be diffed cheaply before the
## full harness exists.
saveRDS(
  lapply(res, function(r) list(ok = r$ok,
                                classes = if(isTRUE(r$ok)) class(r$value) else NULL,
                                msg = if(!isTRUE(r$ok)) r$value else NULL)),
  file.path(Sys.getenv("XVER_OUT"), paste0("smoke-", which_ver, ".rds")))
