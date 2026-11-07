## Producer: runs one scenario set against ONE installed spatialGE version and
## writes a reduced payload plus a stdout transcript.
##
## Usage:
##   XVER_VERSION=main XVER_REP=1 XVER_SCENARIOS=S0,S1 \
##     R_LIBS=$XVER_LIB_MAIN Rscript --vanilla dev/xver/produce.R
##
## Everything version-specific stops here. The payload contains no S4 objects,
## environments or closures (asserted by xv_assert_plain), so the comparator can
## run with spatialGE absent from its library path and a leak fails loudly rather
## than silently deserializing against the wrong class table.

source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "common.R"))
source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "reduce.R"))

which_ver <- Sys.getenv("XVER_VERSION", "")
rep_n     <- as.integer(Sys.getenv("XVER_REP", "1"))
want      <- Sys.getenv("XVER_SCENARIOS", "S0")

spec <- switch(which_ver,
  main = list(lib = Sys.getenv("XVER_LIB_MAIN"), version = "1.2.2"),
  alex = list(lib = Sys.getenv("XVER_LIB_ALEX"), version = "2.0.0"),
  stop("XVER_VERSION must be 'main' or 'alex'"))

info <- xv_load_spatialge(spec$lib, spec$version)
xv_source_driver()
xv_shim_fixture_image()
source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "scenarios.R"))

tag <- sprintf("%s-rep%d", which_ver, rep_n)
scratch <- xv_scratch(tag)
setwd(scratch)

## S12's download is cached outside the per-rep scratch so two reps and two
## versions share one copy -- re-downloading would be four downloads of the same
## bytes and would make the two versions differ on nothing but timing.
thrane_cache <- file.path(Sys.getenv("XVER_ROOT"), "thrane")

requested <- trimws(strsplit(want, ",")[[1]])
registry <- list(
  S0  = function() xv_s0(cores = 1),
  S0b = function() xv_s0(cores = 3),
  S1  = function() xv_s1(file.path(scratch, "formats")),
  S2  = function() xv_s2(),
  S3  = function() xv_s3(),
  S4  = function() xv_s4(),
  S5  = function() xv_s5(),
  S6  = function() xv_s6(),
  S7  = function() xv_s7(),
  S8  = function() xv_s8(),
  S9  = function() xv_s9(),
  S10 = function() xv_s10(),
  S11 = function() xv_s11(),
  S12 = function() xv_s12(thrane_cache),
  S13 = function() xv_s13(),
  S14 = function() xv_s14(),
  S15 = function() xv_s15(file.path(scratch, "h5probe"))
)
unknown <- setdiff(requested, names(registry))
if(length(unknown)) stop("unknown scenario(s): ", paste(unknown, collapse = ", "))

cat(sprintf("== produce: %s spatialGE %s rep%d scenarios=%s ==\n",
            which_ver, info$version, rep_n, paste(requested, collapse = ",")))

payload <- list(
  header = list(
    version_tag = which_ver,
    package_version = info$version,
    lib = "<LIB>",
    rep = rep_n,
    scenarios = requested,
    r_version = as.character(getRversion()),
    ## The comparator refuses to compare payloads whose dependency sets differ:
    ## spatialGE has only ever been used here via load_all(), so main has never
    ## run against these versions and an unguarded diff could be dependency drift.
    dep_versions = local({
      csv <- file.path(Sys.getenv("XVER_OUT"), "dep-versions.csv")
      if(file.exists(csv)) as.list(setNames(read.csv(csv)$version, read.csv(csv)$package)) else list()
    }),
    collate = Sys.getlocale("LC_COLLATE"),
    mc_cores = getOption("mc.cores")
  ),
  scenarios = list()
)

for(sc in requested){
  cat(sprintf("\n-- %s --\n", sc))
  t0 <- Sys.time()
  recs <- tryCatch(registry[[sc]](), error = function(e){
    ## A scenario that cannot even set itself up is recorded as such. It must not
    ## silently vanish: an absent scenario reads as "no difference found".
    list(.scenario_setup_error = list(ok = FALSE, value = conditionMessage(e),
                                       warnings = character(0), messages = character(0)))
  })
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  reduced <- list()
  for(nm in names(recs)){
    rec <- recs[[nm]]
    if(is.null(rec$warnings)) rec$warnings <- character(0)
    if(is.null(rec$messages)) rec$messages <- character(0)
    raw <- tryCatch(xv_reduce_step(rec, canon = FALSE),
                    error = function(e) list(.reduce_error = conditionMessage(e)))
    can <- tryCatch(xv_reduce_step(rec, canon = TRUE),
                    error = function(e) list(.reduce_error = conditionMessage(e)))
    reduced[[nm]] <- list(raw = raw, canon = can)
    flag <- if(!is.null(raw$.reduce_error)) "REDUCE-ERR"
            else if(isTRUE(raw$ok)) "ok" else "err"
    cat(sprintf("   %-34s %-11s %s\n", nm, flag,
                if(isTRUE(raw$ok)) "" else substr(paste(rec$value, collapse=" "), 1, 70)))
  }
  payload$scenarios[[sc]] <- list(steps = reduced, elapsed = elapsed)
}

## The leak assertion. Without it the comparator's isolation is an intention
## rather than a guarantee.
xv_assert_plain(payload)

fp <- file.path(Sys.getenv("XVER_OUT"), sprintf("payload-%s.rds", tag))
saveRDS(payload, fp, version = 3)
cat(sprintf("\nwrote %s (%.1f MB)\n", fp, file.size(fp) / 1e6))
