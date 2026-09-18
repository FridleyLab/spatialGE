## Peak resident-set-size measurement across the whole process tree.
##
## WHY THIS IS NOT bench::mark() OR gc()
##
## The number Phase 2 is judged on is peak RSS of the parent PLUS its mclapply
## children, and neither of the usual tools can see it:
##
## - gc() reports only the calling R process's own heap. A forked child's
##   allocations are invisible to it, which is exactly the memory under study.
## - bench::mark() likewise measures the calling process.
## - /usr/bin/time -l reports max RSS for the parent; children that have been
##   reaped are folded in inconsistently across platforms, and on macOS the
##   figure excludes pages a child touched before exiting.
##
## So this samples `ps` over the process tree from a watcher, at a fixed
## interval, and takes the maximum of the summed RSS. That undercounts if a peak
## falls entirely between two samples, so the interval is deliberately short
## relative to the work being measured, and every reported figure is the max over
## several repeats.
##
## IMPORTANT for interpreting fork numbers. Fork is copy-on-write: a child's RSS
## initially shares pages with the parent and `ps` counts those pages in BOTH
## processes, so summed RSS over-reports real usage at first and converges toward
## the truth as the GC touches pages and forces genuine copies. That over-report
## is the thing Phase 2 reduces -- a worker that receives only its own sample has
## far fewer shared pages to fault -- so the metric is the right one even though
## its absolute value is not a clean "bytes allocated".
##
## Usage:
##   Rscript dev/profiling/rss_peak.R            # runs the standard comparison
##   source('dev/profiling/rss_peak.R'); rss_peak(quote(expr), cores=4)

## Sum RSS in KB over a process and all its descendants.
rss_tree_kb <- function(root_pid){
  ## One ps call for the whole table, then walk the pid/ppid graph -- calling ps
  ## per pid would itself perturb the measurement.
  out <- tryCatch(
    system2("ps", c("-eo", "pid=,ppid=,rss="), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if(length(out) == 0) return(NA_real_)

  parts <- strsplit(trimws(out), "[[:space:]]+")
  parts <- parts[vapply(parts, length, integer(1)) == 3]
  pid  <- as.integer(vapply(parts, `[`, character(1), 1))
  ppid <- as.integer(vapply(parts, `[`, character(1), 2))
  rss  <- as.numeric(vapply(parts, `[`, character(1), 3))

  keep <- root_pid
  repeat {
    children <- pid[ppid %in% keep & !(pid %in% keep)]
    if(length(children) == 0) break
    keep <- c(keep, children)
  }
  sum(rss[pid %in% keep], na.rm = TRUE)
}

## Runs `expr` in a FRESH Rscript child and samples that child's process tree.
## A fresh process matters: peak RSS is a high-water mark, so anything the
## current session already did would contaminate it permanently.
rss_peak <- function(expr_text, setup_text = "", interval = 0.02, reps = 3){
  script <- tempfile(fileext = ".R")
  writeLines(c(
    ## pkgload rather than devtools::load_all -- but note this barely matters.
    ## Measured floor is ~1355 MB either way: the overhead is spatialGE's own 44
    ## Imports (EBImage, GSVA, ComplexHeatmap, spaMM, sf, sp, arrow,
    ## DelayedArray, ...), not the loader. It is irreducible here, so every figure
    ## this script reports is a DELTA above that floor. Absolute peaks are
    ## meaningless for judging a change -- a 4 GB fork cost looks like a 40%
    ## increase if you quote the total.
    "suppressMessages(pkgload::load_all('.', quiet=TRUE, helpers=FALSE))",
    "suppressMessages(source('inst/testdata/make_scale_fixtures.R'))",
    setup_text,
    "invisible(suppressWarnings({",
    expr_text,
    "}))"
  ), script)

  peaks <- numeric(reps)
  for(r in seq_len(reps)){
    p <- processx_like_spawn(script)
    peak <- 0
    repeat {
      alive <- is.na(p$get_exit_status())
      cur <- rss_tree_kb(p$get_pid())
      if(!is.na(cur)) peak <- max(peak, cur)
      if(!alive) break
      Sys.sleep(interval)
    }
    peaks[r] <- peak
  }
  list(peak_mb = max(peaks) / 1024, all_mb = peaks / 1024)
}

## Minimal spawn wrapper so this works without taking a processx dependency.
processx_like_spawn <- function(script){
  pid_file <- tempfile()
  ## Launch detached, recording the child's pid so the tree walk has a root.
  system2("bash", c("-c", shQuote(sprintf(
    "Rscript %s >/dev/null 2>&1 & echo $! > %s; wait $!",
    shQuote(script), shQuote(pid_file)))), wait = FALSE)

  ## Wait for the pid file to appear.
  for(i in 1:500){
    if(file.exists(pid_file) && nzchar(readLines(pid_file, warn = FALSE)[1])) break
    Sys.sleep(0.01)
  }
  pid <- as.integer(readLines(pid_file, warn = FALSE)[1])

  list(
    get_pid = function() pid,
    get_exit_status = function(){
      alive <- system2("kill", c("-0", pid), stdout = FALSE, stderr = FALSE) == 0
      if(alive) NA_integer_ else 0L
    }
  )
}

## The standard comparison. Reports delta above the measured floor, and the
## fork cost -- peak at N cores minus peak at 1 core -- which is the number
## Phase 2 exists to reduce.
rss_report <- function(steps, cores_list = c(1, 4), reps = 3){
  floor_mb <- rss_peak("invisible(NULL)", reps = reps)$peak_mb
  cat(sprintf("floor (package load only): %.0f MB -- subtracted from every figure below\n\n", floor_mb))
  cat(sprintf("%-26s %6s %12s %12s\n", "step", "cores", "above_floor", "fork_cost"))

  for(nm in names(steps)){
    base1 <- NA_real_
    for(cores in cores_list){
      res <- rss_peak(sprintf(steps[[nm]], cores), reps = reps)
      above <- res$peak_mb - floor_mb
      if(cores == 1) base1 <- above
      cat(sprintf("%-26s %6d %12.1f %12s\n", nm, cores, above,
                  if(cores == 1) "-" else sprintf("%+.1f", above - base1)))
    }
  }
  invisible(floor_mb)
}

STEPS <- list(
  `transform_data[visium]` = "st <- fixture_scale_stlist(preset='visium', n_samples=4); transform_data(st, method='log', cores=%d)",
  `transform_data[xenium]` = "st <- fixture_scale_stlist(preset='xenium', n_samples=4); transform_data(st, method='log', cores=%d)"
)

if(sys.nframe() == 0){
  cat("Peak RSS over the process tree (parent + mclapply children)\n")
  cat("interval 20ms, max of 3 repeats\n\n")
  rss_report(STEPS)
}
