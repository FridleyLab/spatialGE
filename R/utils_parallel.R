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
# @param n size of the work unit, passed to count_cores() when cores=NULL.
# @param cores user-requested cores, or NULL to auto-detect.
# @param coerce how to coerce a user-supplied cores: 'as_integer' (may
# produce NA -- the caller is responsible for checking and raising, as
# before), 'ceiling', or 'none' (used as-is, no coercion).
# @return an integer (or NA if coerce='as_integer' could not parse cores).
#
resolve_cores = function(n, cores=NULL, coerce=c('as_integer', 'ceiling', 'none')){
  coerce = match.arg(coerce)

  if(.Platform$OS.type == 'windows'){
    return(1L)
  }

  override = sge_cores_override()
  if(!is.null(override)){
    return(override)
  }

  if(is.null(cores)){
    return(count_cores(n))
  }

  switch(coerce,
        as_integer=as.integer(cores),
        ceiling=ceiling(cores),
        none=cores)
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
sge_lapply = function(X, FUN, ..., cores=NULL, preschedule=TRUE){
  override = sge_cores_override()
  if(!is.null(override)){
    cores = override
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
