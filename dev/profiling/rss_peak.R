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
  sel <- pid %in% keep
  ## Both are reported because the SUM is badly misleading under fork and the MAX
  ## is not. A forked child shares the parent's pages copy-on-write and `ps`
  ## charges those pages to every process, so the sum counts the ~1.35 GB of
  ## loaded packages once per worker -- which is where an apparent "+5 GB for 4
  ## workers" comes from even when almost nothing was really allocated. The max
  ## over single processes is the honest proxy for what one worker costs.
  c(sum = sum(rss[sel], na.rm = TRUE),
    max = if(any(sel)) max(rss[sel], na.rm = TRUE) else NA_real_,
    n = sum(sel))
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

  peak_sum <- numeric(reps); peak_max <- numeric(reps); peak_n <- numeric(reps)
  for(r in seq_len(reps)){
    p <- processx_like_spawn(script)
    ps_ <- 0; pm <- 0; pn <- 0
    repeat {
      alive <- is.na(p$get_exit_status())
      cur <- rss_tree_kb(p$get_pid())
      if(!is.na(cur[["sum"]])){
        ps_ <- max(ps_, cur[["sum"]])
        pm <- max(pm, cur[["max"]], na.rm = TRUE)
        pn <- max(pn, cur[["n"]])
      }
      if(!alive) break
      Sys.sleep(interval)
    }
    peak_sum[r] <- ps_; peak_max[r] <- pm; peak_n[r] <- pn
  }
  list(peak_mb = max(peak_sum) / 1024,
       peak_proc_mb = max(peak_max) / 1024,
       max_procs = max(peak_n),
       all_mb = peak_sum / 1024)
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
  fl <- rss_peak("invisible(NULL)", reps = reps)
  floor_mb <- fl$peak_mb
  cat(sprintf("floor (package load only): %.0f MB\n", floor_mb))
  cat("tree_sum  = summed RSS over parent + children, MINUS floor. Inflated by\n")
  cat("            copy-on-write double counting: shared pages are charged to\n")
  cat("            every process, so this grows ~floor per extra worker even when\n")
  cat("            nothing is really allocated. Do NOT read it as memory used.\n")
  cat("worst_proc = peak RSS of the single largest process. The honest proxy for\n")
  cat("            what one worker actually costs.\n\n")
  cat(sprintf("%-26s %6s %10s %12s %6s\n", "step", "cores", "tree_sum", "worst_proc", "procs"))

  for(nm in names(steps)){
    for(cores in cores_list){
      res <- rss_peak(sprintf(steps[[nm]], cores), reps = reps)
      cat(sprintf("%-26s %6d %10.1f %12.1f %6d\n", nm, cores,
                  res$peak_mb - floor_mb, res$peak_proc_mb, res$max_procs))
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
