##
# @title sge_cores_override
# @description Checks the package-wide override for the number of parallel
# cores to use (`options(spatialGE.cores=)`, then env var `SPATIALGE_CORES`).
# Shared by resolve_cores(), sge_lapply(), and sge_bpparam() so that setting
# either one forces every parallelized call site in the package to run
# serially and in-process, with no other source changes.
# @return an integer, or NULL if no override is set.
#
sge_cores_override = function(){
  opt = getOption('spatialGE.cores', NULL)
  if(!is.null(opt)){
    return(as.integer(opt))
  }
  env = Sys.getenv('SPATIALGE_CORES', unset=NA)
  if(!is.na(env) && nzchar(env)){
    return(as.integer(env))
  }
  return(NULL)
}


##
# @title resolve_cores
# @description Resolves the number of cores for a parallelized call site,
# replacing the boilerplate that used to be duplicated at every such site:
# force 1 core on Windows, auto-detect via count_cores(n) when the user did
# not request a specific number, otherwise coerce whatever the user passed.
#
# Validation of an unusable `cores` value lives here too. It used to be the
# caller's job and was done inconsistently: of the 17 call sites, some raised
# error0032, others error0036, five raised an uncatalogued bare stop(), and
# transform_data()'s gene_meta path did not check at all -- so `cores='abc'`
# there became NA and was handed straight to sge_lapply(). Resolving it in one
# place makes the behavior uniform and removes six copies of the same guard.
#
# @param n size of the work unit, passed to count_cores() when cores=NULL.
# @param cores user-requested cores, or NULL to auto-detect. A value that
# cannot be resolved raises error0093.
# @param coerce how to coerce a user-supplied cores: 'as_integer', 'ceiling',
# or 'none' (used as-is, no coercion).
# @return the resolved core count. The caller's coercion is preserved: a
# non-integral value survives for coerce='ceiling'/'none'.
#
resolve_cores = function(n, cores=NULL, coerce=c('as_integer', 'ceiling', 'none')){
  coerce = match.arg(coerce)

  # The override is consulted FIRST, before any platform decision. It used to sit
  # after a `return(1L)` for Windows, so options(spatialGE.cores=) and
  # SPATIALGE_CORES were silently ignored there -- a developer aid that did not
  # work on the one platform where the parallel behavior differs most.
  override = sge_cores_override()
  if(!is.null(override)){
    return(check_resolved_cores(override))
  }

  # Windows has no fork(), so mclapply falls back to serial. It is no longer
  # forced to 1 core here, because sge_lapply() now has a non-forking path for it
  # (see sge_use_psock below). Returning 1L unconditionally made the whole package
  # serial on Windows regardless of what the user asked for.

  if(is.null(cores)){
    return(check_resolved_cores(count_cores(n)))
  }

  # suppressWarnings: as.integer('abc') warns "NAs introduced by coercion" and
  # then the check below raises a catalogued error that actually explains what to
  # do. Leaving the warning in place printed both.
  resolved = suppressWarnings(switch(coerce,
        as_integer=as.integer(cores),
        ceiling=ceiling(cores),
        none=cores))

  check_resolved_cores(resolved)
}


##
# @title check_resolved_cores
# @description Single guard for the result of cores resolution: it must be one
# finite, positive number. Kept separate so every return path in resolve_cores()
# passes through it.
#
# Validates only -- the value is returned unchanged. coerce='ceiling' and
# 'none' deliberately preserve a non-integral value (2.5 stays 2.5), and
# as.integer() here would silently truncate it to 2, breaking the fidelity rule
# that non-mclapply consumers depend on. This is also why validation cannot live
# in sge_lapply(): it never sees these call sites.
#
# @param resolved the candidate core count.
# @return `resolved`, unchanged, if usable; otherwise raises error0093.
#
check_resolved_cores = function(resolved){
  ok = length(resolved) == 1 && !is.na(resolved) &&
       is.finite(suppressWarnings(as.numeric(resolved))) &&
       suppressWarnings(as.numeric(resolved)) >= 1
  if(!ok){
    raise_err(err_code='error0093')
  }
  resolved
}


##
# @title sge_lapply
# @description Thin wrapper over parallel::mclapply. Forwards straight to
# mclapply rather than reimplementing its dispatch, so its own semantics
# (falling back to plain, in-process lapply when cores < 2 or length(X) < 2,
# mc.reset.stream(), the .check_ncores() guard, etc.) are inherited rather
# than duplicated. This is the single choke point through which every
# parallelized call site now runs, so it is also where the
# spatialGE.cores/SPATIALGE_CORES override is enforced -- including for call
# sites that never resolved a `cores` argument in the first place.
# @param cores number of cores to request. NULL forwards mclapply's own
# default (i.e. mc.cores is simply not passed).
# @param preschedule passed through as mc.preschedule.
#
##
# @title sge_use_psock
# @description Whether sge_lapply() should take its non-forking (PSOCK) path.
# TRUE on Windows, where fork() does not exist and mclapply silently degrades to
# serial. Overridable with options(spatialGE.psock=TRUE) so the path is testable
# on Unix -- PSOCK works everywhere, it is just not the better choice where fork
# is available.
#
sge_use_psock = function(){
  opt = getOption('spatialGE.psock', NULL)
  if(!is.null(opt)){
    return(isTRUE(opt))
  }
  .Platform$OS.type == 'windows'
}


##
# @title sge_lapply_psock
# @description The non-forking backend. Creates a PSOCK cluster, loads spatialGE
# in each worker, and parLapply()s over X.
#
# WHY THIS IS NOT JUST mclapply WITH A DIFFERENT NAME
#
# Forking shares the parent's memory copy-on-write, so a worker gets the whole
# session for free. PSOCK workers are fresh R processes: every element of X and
# everything FUN's environment reaches is SERIALIZED and sent per task. That makes
# the payload-slicing work load-bearing here rather than cosmetic -- a call site
# whose closure still captured the whole STlist would ship it to every worker,
# turning a free fork into a real copy. It is also why this path is used only
# where fork is unavailable.
#
# Workers need the package. Under an installed spatialGE, library() suffices;
# under devtools/pkgload during development the package is not installed, so the
# loaded namespace's path is passed through and pkgload is used instead. Without
# that, FUN's references to package internals fail on the worker with
# "could not find function".
#
sge_lapply_psock = function(X, FUN, ..., cores){
  cores = max(1L, as.integer(cores))
  if(cores < 2 || length(X) < 2){
    return(lapply(X=X, FUN=FUN, ...))
  }

  pkg_lib = tryCatch(getNamespaceInfo('spatialGE', 'path'), error=function(e) NULL)

  # Detecting a development load: pkgload/devtools place a .__DEVTOOLS__ object in
  # the namespace, which is a dependency-free marker. Testing for a DESCRIPTION
  # file does NOT work -- a dev-loaded package's namespace path is the source
  # directory, which has one, so that check reports "installed" and the worker
  # then fails with "there is no package called 'spatialGE'".
  is_dev = !is.null(tryCatch(
    get0('.__DEVTOOLS__', envir=asNamespace('spatialGE'), inherits=FALSE),
    error=function(e) NULL))

  cl = parallel::makePSOCKcluster(cores)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  parallel::clusterExport(cl, varlist=c('pkg_lib', 'is_dev'), envir=environment())
  parallel::clusterEvalQ(cl, {
    if(is_dev){
      # Development load: the package is not installed, so reproduce the
      # in-development namespace on the worker from the source path.
      suppressMessages(pkgload::load_all(pkg_lib, quiet=TRUE, helpers=FALSE))
    } else{
      suppressMessages(library('spatialGE'))
    }
    NULL
  })

  parallel::parLapply(cl=cl, X=X, fun=FUN, ...)
}


sge_lapply = function(X, FUN, ..., cores=NULL, preschedule=TRUE){
  override = sge_cores_override()
  if(!is.null(override)){
    cores = override
  }

  if(sge_use_psock()){
    # preschedule has no PSOCK analogue: parLapply chunks X evenly and there is
    # no per-element-fork mode to mirror preschedule=FALSE. Deliberately ignored
    # rather than faked, and fidelity rule 4 (per-call-site preschedule semantics)
    # therefore does not carry over to this backend.
    return(sge_lapply_psock(X=X, FUN=FUN, ..., cores=if(is.null(cores)) 1L else cores))
  }

  if(is.null(cores)){
    parallel::mclapply(X=X, FUN=FUN, ..., mc.preschedule=preschedule)
  } else{
    parallel::mclapply(X=X, FUN=FUN, ..., mc.cores=cores, mc.preschedule=preschedule)
  }
}


##
# @title sge_bpparam
# @description Returns the BiocParallel back-end to use, honoring the same
# spatialGE.cores/SPATIALGE_CORES override as resolve_cores()/sge_lapply().
# Returns SerialParam() whenever cores resolve to less than 2, so a stray
# MulticoreParam(workers=1) never forks a worker process.
#
sge_bpparam = function(cores=NULL){
  override = sge_cores_override()
  if(!is.null(override)){
    cores = override
  }
  if(is.null(cores) || cores < 2){
    BiocParallel::SerialParam()
  } else{
    BiocParallel::MulticoreParam(workers=cores)
  }
}
