##
# @title assert_stlist
# @description Validates that `x` is an STlist, replacing the boilerplate
# duplicated across many top-level functions. Some call sites only check
# `is.null(x)` and never checked the class (a pre-existing gap), and one
# checks the class but not `is.null()` explicitly (harmless, since
# `is(NULL, 'STlist')` is already FALSE) -- `check_null`/`check_class` let
# each call site keep its original condition exactly.
# @param x the object to validate.
# @param message the stop() message to use when err_code is NULL.
# @param err_code if supplied, raise_err(err_code=) is used instead of stop().
# @param check_null whether to fail on is.null(x).
# @param check_class whether to fail on !is(x, 'STlist').
#
assert_stlist = function(x, message='The input must be a STlist.', err_code=NULL,
                          check_null=TRUE, check_class=TRUE){
  # Accept a Seurat object anywhere an STlist is accepted, by coercing here.
  #
  # This is the "dispatch, then coerce at the boundary" decision from the plan's
  # Phase 7, implemented at the single validation helper rather than as 22 S4
  # generics. The reason is that the bodies read slots directly, so a genuine
  # Seurat method would mean maintaining two implementations of everything; and
  # since every entry point has to validate its input anyway, the validation
  # helper is already the boundary. An S4 layer whose methods only coerced would
  # add indirection without adding behavior.
  #
  # Consequence worth knowing: coercion runs STlist() on every call, so calling
  # several functions in sequence on the same Seurat object re-imports it each
  # time. Convert once with STlist(rnacounts=<seurat>) for any real workflow.
  if(!is.null(x) && !is(x, 'STlist') && inherits(x, 'Seurat')){
    x = STlist(rnacounts=x, verbose=FALSE)
  }

  failed = (check_null && is.null(x)) || (check_class && !is(x, 'STlist'))
  if(failed){
    if(!is.null(err_code)){
      raise_err(err_code=err_code)
    } else{
      stop(message)
    }
  }

  # Returned so call sites can pick up the coerced object: `x = assert_stlist(x)`.
  # Invisible, so the pre-existing statement-form calls stay silent.
  invisible(x)
}


##
# @title assert_transformed
# @description Validates that `x@tr_counts` has been populated (i.e.
# transform_data() has already been run), replacing the duplicated
# `rlang::is_empty(x@tr_counts)` / `length(x@tr_counts) < 1` checks.
# @param x an STlist.
# @param message the stop() message to use when err_code is NULL.
# @param err_code if supplied, raise_err(err_code=) is used instead of stop().
#
assert_transformed = function(x, message='There are no transformed counts in this STlist.', err_code=NULL){
  if(rlang::is_empty(x@tr_counts)){
    if(!is.null(err_code)){
      raise_err(err_code=err_code)
    } else{
      stop(message)
    }
  }
}


##
# @title resolve_samples
# @description Resolves a user-supplied `samples` argument (NULL, numeric
# indices, or a character vector) against the names of one slot of an
# STlist. The ~14 call sites this replaces are NOT uniform, so every axis of
# divergence is a parameter rather than a fixed choice, and each call site
# must be migrated with the parameters that reproduce its exact prior
# behavior -- including the ones that never validated character `samples`
# at all, and the ones that use a buggy substring/regex `grep()` existence
# check (against a single '|'-joined pattern) instead of an exact `%in%`
# match. Do not "fix" any of this while migrating a call site; that is
# Phase 11 territory, tracked and tested separately.
# @param x an STlist.
# @param samples NULL, a numeric vector of indices, or a character vector
# of names, exactly as accepted at the original call site.
# @param slot the STlist slot to resolve names from (e.g. 'spatial_meta',
# 'tr_counts', 'counts').
# @param numeric_handling how numeric `samples` are converted to names:
# 'na_omit' wraps in `as.vector(na.omit(...))`, 'plain' does not.
# @param char_transform how character `samples` are handled: 'exact' filters
# via `%in%`; 'none' leaves `samples` untouched (no filtering, no check);
# 'validate_all' stops/raises unless every requested name is present, but
# (matching the original) does not filter `samples` down.
# @param post_check an additional check applied to the resolved `samples`
# regardless of whether they came from the numeric or character branch:
# 'none' skips it; 'any_in' stops/raises when `length(samples) == 0` or
# none of `samples` are in the available names; 'grep' reproduces the
# known bug (`grep(paste0(samples, collapse='|'), available)`) and does
# NOT filter `samples` -- only checks that at least one loose match exists.
# @param message the stop() message to use when err_code is NULL.
# @param err_code if supplied, raise_err(err_code=) is used instead of stop().
# @return a character vector of resolved sample names (or, for
# char_transform='none'/'validate_all' and post_check='grep', the original
# `samples` largely as-is -- reproducing each site's prior behavior).
#
resolve_samples = function(x, samples=NULL, slot='spatial_meta',
                            numeric_handling=c('na_omit', 'plain'),
                            char_transform=c('exact', 'none', 'validate_all'),
                            post_check=c('none', 'any_in', 'grep'),
                            message=NULL, err_code=NULL){
  numeric_handling = match.arg(numeric_handling)
  char_transform = match.arg(char_transform)
  post_check = match.arg(post_check)

  available = names(methods::slot(x, slot))

  if(is.null(samples)){
    return(available)
  }

  raise = function(){
    if(!is.null(err_code)){
      raise_err(err_code=err_code)
    } else{
      stop(message)
    }
  }

  if(is.numeric(samples)){
    if(numeric_handling == 'na_omit'){
      samples = as.vector(na.omit(available[samples]))
    } else{
      samples = available[samples]
    }
  } else{
    if(char_transform == 'exact'){
      samples = samples[samples %in% available]
    } else if(char_transform == 'validate_all'){
      if(!all(samples %in% available)){
        raise()
      }
    }
    # char_transform == 'none': samples left untouched
  }

  if(post_check == 'any_in'){
    if(length(samples) == 0 | !any(samples %in% available)){
      raise()
    }
  } else if(post_check == 'grep'){
    if(length(grep(paste0(samples, collapse='|'), available)) == 0){
      raise()
    }
  }

  return(samples)
}
