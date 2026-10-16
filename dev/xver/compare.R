## The comparator. Runs with spatialGE ABSENT from the library path: every
## payload is plain data, so an accidental S4 or environment leak fails loudly
## here instead of silently deserializing against the wrong class table.
##
## Usage:  Rscript --vanilla dev/xver/compare.R
##
## Four passes plus two gates:
##   Pass 0   rep1 vs rep2 within each version   -> the noise allowlist
##   Pass A   main vs alex, as produced          -> everything, including order
##   Pass B   main vs alex, canonicalized        -> value differences only
##   Pass B'  Pass B at relative tolerance 1e-12 -> material value differences
##   gate C   name-set equality per matrix/frame -> makes Pass B a keyed comparison
##   gate P   cluster partition equality         -> guards the downstream cascade

source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "reduce.R"))
source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "classify.R"))

if("spatialGE" %in% loadedNamespaces()){
  stop("spatialGE is loaded in the comparator; run it with a library path that excludes both version libraries")
}

OUT <- Sys.getenv("XVER_OUT")

xv_read <- function(tag){
  fp <- file.path(OUT, sprintf("payload-%s.rds", tag))
  if(!file.exists(fp)) return(NULL)
  readRDS(fp)
}

## ---- leaf comparison ------------------------------------------------------

## Returns NULL when equal, or a record describing how the two leaves differ.
## tolerance NULL means exact.
##
## When a tolerance is given it is applied SCALE-RELATIVE: max|a-b| divided by the
## magnitude of the vector as a whole, not elementwise. Elementwise relative
## difference is the wrong question for anything with a tail. A kernel density
## curve evaluated on the same data in a different summation order agrees to ~1e-16
## in absolute terms, but its tail values are ~1e-20, so an elementwise ratio there
## reports 0.11 -- or 1.0 -- and calls a numerically identical curve a regression.
xv_leaf_diff <- function(a, b, tolerance = NULL){
  if(is.null(a) && is.null(b)) return(NULL)
  if(xor(is.null(a), is.null(b))) return(list(kind = "null_vs_value"))
  if(length(a) != length(b)){
    return(list(kind = "length", left = length(a), right = length(b)))
  }
  if(is.numeric(a) && is.numeric(b)){
    na <- is.na(a); nb <- is.na(b)
    if(!identical(na, nb)) return(list(kind = "na_pattern", n = sum(na != nb)))
    ok <- !na
    if(!any(ok)) return(NULL)
    av <- a[ok]; bv <- b[ok]
    if(identical(av, bv)) return(NULL)
    d <- abs(av - bv)
    elem <- pmax(abs(av), abs(bv))
    rel <- ifelse(elem > 0, d / elem, d)
    ## Vector-scale relative difference: how big is the worst disagreement
    ## compared with the magnitude of the data it sits in.
    vscale <- max(elem)
    scaled <- if(vscale > 0) max(d) / vscale else max(d)
    if(!is.null(tolerance) && scaled <= tolerance) return(NULL)
    if(max(d) == 0) return(NULL)
    return(list(kind = "numeric", n = sum(d > 0), max_rel = max(rel),
                max_abs = max(d), max_scaled = scaled))
  }
  ## Character/logical/integer and anything else: exact.
  if(identical(a, b)) return(NULL)
  if(is.character(a) && is.character(b)){
    n <- sum(a != b | is.na(a) != is.na(b))
    ## A pure permutation of the same multiset is worth naming: it is the
    ## signature of the ordering fix rather than a content change.
    same_set <- setequal(a, b) && identical(sort(table(a)), sort(table(b)))
    return(list(kind = if(same_set) "permutation" else "character",
                n = n,
                example_left = paste(utils::head(a[a != b], 3), collapse = ","),
                example_right = paste(utils::head(b[a != b], 3), collapse = ",")))
  }
  list(kind = "value", n = sum(a != b, na.rm = TRUE))
}

## ---- flat-dict comparison -------------------------------------------------

xv_flat_payload <- function(payload, mode){
  out <- list()
  for(sc in names(payload$scenarios)){
    steps <- payload$scenarios[[sc]]$steps
    for(st in names(steps)){
      node <- steps[[st]][[mode]]
      fl <- xv_flatten_list(node)
      names(fl) <- paste(sc, st, names(fl), sep = "/")
      out <- c(out, fl)
    }
  }
  out
}

xv_diff_flat <- function(a, b, tolerance = NULL){
  ka <- names(a); kb <- names(b)
  diffs <- list()
  for(k in setdiff(ka, kb)) diffs[[k]] <- list(kind = "only_main")
  for(k in setdiff(kb, ka)) diffs[[k]] <- list(kind = "only_alex")
  for(k in intersect(ka, kb)){
    d <- xv_leaf_diff(a[[k]], b[[k]], tolerance)
    if(!is.null(d)) diffs[[k]] <- d
  }
  diffs
}

xv_scenario_of <- function(path) sub("/.*$", "", path)

## ---- Pass 0: noise floor --------------------------------------------------

cat("== Pass 0: self-reproducibility within each version ==\n")
noise <- character(0)
pass0 <- list()
REP_A <- Sys.getenv("XVER_REP_A", "1"); REP_B <- Sys.getenv("XVER_REP_B", "2")
for(v in c("main", "alex")){
  p1 <- xv_read(sprintf("%s-rep%s", v, REP_A)); p2 <- xv_read(sprintf("%s-rep%s", v, REP_B))
  if(is.null(p1) || is.null(p2)){
    cat(sprintf("  %-5s SKIPPED (need both reps)\n", v)); next
  }
  d <- xv_diff_flat(xv_flat_payload(p1, "raw"), xv_flat_payload(p2, "raw"))
  pass0[[v]] <- d
  noise <- union(noise, names(d))
  cat(sprintf("  %-5s %d non-reproducible leaf path(s)\n", v, length(d)))
  ## Never a silent allowlist: an unprinted noise mask is how a harness rots.
  for(k in utils::head(names(d), 40)) cat(sprintf("      %s  [%s]\n", k, d[[k]]$kind))
  if(length(d) > 40) cat(sprintf("      ... and %d more\n", length(d) - 40))
}

## ---- load the two versions ------------------------------------------------

## Which rep pair to compare, and which pair forms the noise floor. Defaults to
## rep1 vs rep1 with rep1/rep2 as the floor; the extra-scenario batch uses rep9/rep10.
CMP_REP <- Sys.getenv("XVER_CMP_REP", "1")
pm <- xv_read(sprintf("main-rep%s", CMP_REP)); pa <- xv_read(sprintf("alex-rep%s", CMP_REP))
if(is.null(pm) || is.null(pa))
  stop(sprintf("need payload-main-rep%s.rds and payload-alex-rep%s.rds", CMP_REP, CMP_REP))

## Dependency-set gate. Both sides must have run against the same dependency
## versions or a difference could be dependency drift rather than the refactor.
dm <- pm$header$dep_versions; da <- pa$header$dep_versions
drift <- names(dm)[vapply(names(dm), function(k) !identical(dm[[k]], da[[k]]), logical(1))]
if(length(drift)) stop("dependency versions differ between payloads: ", paste(drift, collapse = ", "))
if(!identical(pm$header$collate, pa$header$collate)) stop("LC_COLLATE differs between payloads")
cat(sprintf("\ndependency set identical (%d packages), collate=%s\n",
            length(dm), pm$header$collate))

common_sc <- intersect(names(pm$scenarios), names(pa$scenarios))
cat(sprintf("scenarios compared: %s\n", paste(common_sc, collapse = ", ")))

## ---- Passes A, B, B' ------------------------------------------------------

fa_raw <- xv_flat_payload(pm, "raw");   fb_raw <- xv_flat_payload(pa, "raw")
fa_can <- xv_flat_payload(pm, "canon"); fb_can <- xv_flat_payload(pa, "canon")

pass_a  <- xv_diff_flat(fa_raw, fb_raw)
pass_b  <- xv_diff_flat(fa_can, fb_can)
pass_bp <- xv_diff_flat(fa_can, fb_can, tolerance = 1e-12)

cat(sprintf("\n== Pass A (raw)      %d differing leaf path(s) of %d\n", length(pass_a), length(fa_raw)))
cat(sprintf("== Pass B (canon)    %d differing leaf path(s) of %d\n", length(pass_b), length(fa_can)))
cat(sprintf("== Pass B' (tol1e-12) %d differing leaf path(s)\n", length(pass_bp)))

## ---- gate C: name-set equality -------------------------------------------

## Canonical comparison sorts both sides by name, so comparing element-wise IS a
## keyed comparison -- but only if the two sides carry the SAME names. Without
## this gate a genuinely dropped or duplicated spot or gene could be absorbed by
## the sort instead of reported.
cat("\n== gate C: name-set equality (spots and genes) ==\n")
gateC <- list()
name_keys <- grep("/(rownames|colnames)$", intersect(names(fa_can), names(fb_can)), value = TRUE)
for(k in name_keys){
  a <- fa_can[[k]]; b <- fb_can[[k]]
  if(!is.character(a) || !is.character(b)) next
  if(!setequal(a, b)){
    gateC[[k]] <- list(only_main = setdiff(a, b), only_alex = setdiff(b, a))
  }
}
## Scoped, not global. A name-set mismatch invalidates the order-only reasoning for
## the OBJECT it belongs to, not for the whole run: treating it as a global switch
## would reclassify every genuinely order-only path in every other scenario as a
## regression, which buries the one real finding under ten thousand false ones.
gateC_prefixes <- unique(sub("/(rownames|colnames)$", "", names(gateC)))
if(!length(gateC)){
  cat(sprintf("  ok  %d name vector(s) carry identical sets on both versions\n", length(name_keys)))
} else {
  for(k in names(gateC)){
    cat(sprintf("  MISMATCH %s\n    only main: %s\n    only alex: %s\n", k,
                paste(utils::head(gateC[[k]]$only_main, 5), collapse = ","),
                paste(utils::head(gateC[[k]]$only_alex, 5), collapse = ",")))
  }
  cat(sprintf("  -> order-only reasoning suspended for %d object prefix(es)\n",
              length(gateC_prefixes)))
}

## ---- gate P: cluster partition -------------------------------------------

## STclust's labels permute under a spot reordering because hclust breaks ties by
## input order. Labels are COLUMNS, so a canonical sort would absorb the
## relabeling and report clean -- while STenrich/STgradient/STdiff downstream were
## handed different annot=/domain=/ref= values and answered different questions.
## Compare the partition (the set of spot-sets) rather than the labels.
cat("\n== gate P: cluster partition equality ==\n")
xv_partitions <- function(payload, mode){
  out <- list()
  for(sc in names(payload$scenarios)){
    steps <- payload$scenarios[[sc]]$steps
    for(st in names(steps)){
      node <- steps[[st]][[mode]]
      sm <- node$value$spatial_meta
      if(is.null(sm)) next
      for(s in names(sm)){
        cols <- sm[[s]]$columns
        lib <- cols$libname
        for(cn in grep("^stclust_", names(cols), value = TRUE)){
          lab <- as.character(cols[[cn]])
          parts <- split(lib, lab)
          ## Canonical form of a partition: each block sorted, blocks ordered by
          ## their own minimum member. Independent of what the blocks are called.
          parts <- lapply(parts, sort)
          parts <- parts[order(vapply(parts, function(p) p[1], character(1)))]
          out[[paste(sc, st, s, cn, sep = "/")]] <- unname(parts)
        }
      }
    }
  }
  out
}
partm <- xv_partitions(pm, "raw"); parta <- xv_partitions(pa, "raw")
gateP <- list()
for(k in intersect(names(partm), names(parta))){
  if(!identical(partm[[k]], parta[[k]])) gateP[[k]] <- TRUE
}
cascade_scenarios <- character(0)
if(!length(partm)){
  cat("  (no stclust_ columns present)\n")
} else if(!length(gateP)){
  cat(sprintf("  ok  %d clustering(s) have identical partitions on both versions\n", length(partm)))
  ## Labels may still differ; that is ORDER_INDUCED_VALUE, not a partition change.
  lbl <- 0L
  for(k in intersect(names(partm), names(parta))){
    if(!identical(names(partm[[k]]), names(parta[[k]]))) lbl <- lbl + 1L
  }
} else {
  for(k in names(gateP)) cat(sprintf("  PARTITION DIFFERS %s\n", k))
  cascade_scenarios <- unique(xv_scenario_of(names(gateP)))
  cat(sprintf("  -> downstream results in scenario(s) %s are INCOMPARABLE, not order-only\n",
              paste(cascade_scenarios, collapse = ", ")))
}

## ---- gate S: selection sets -----------------------------------------------

## STclust, STdiff and STgradient all pick their working gene set with
## arrange(desc(vst.variance.standardized)) then head(topgenes). The variance
## statistic is computed by a reduction over spots, so a spot permutation moves its
## last bits -- and a 1e-16 move can flip a rank across the topgenes cut. When that
## happens the two versions computed their statistics on DIFFERENT genes, so
## comparing the result rows is meaningless rather than merely noisy.
cat("\n== gate S: gene-selection sets in result tables ==\n")
xv_selection_sets <- function(payload){
  out <- list()
  for(sc in names(payload$scenarios)){
    steps <- payload$scenarios[[sc]]$steps
    for(st in names(steps)){
      fl <- xv_flatten_list(steps[[st]][["canon"]])
      for(k in grep("/columns/gene$", names(fl), value = TRUE)){
        v <- fl[[k]]
        if(is.character(v)) out[[paste(sc, st, k, sep = "|")]] <- sort(unique(v))
      }
    }
  }
  out
}
selm <- xv_selection_sets(pm); sela <- xv_selection_sets(pa)
gateS <- list()
for(k in intersect(names(selm), names(sela))){
  if(!setequal(selm[[k]], sela[[k]])) gateS[[k]] <- list(
    only_main = setdiff(selm[[k]], sela[[k]]), only_alex = setdiff(sela[[k]], selm[[k]]))
}
## Prefixes whose downstream leaves cannot be compared row-for-row.
selection_prefixes <- character(0)
selection_steps <- character(0)
if(!length(selm)){
  cat("  (no gene columns present)\n")
} else if(!length(gateS)){
  cat(sprintf("  ok  %d gene-selection set(s) identical on both versions\n", length(selm)))
} else {
  for(k in names(gateS)){
    parts <- strsplit(k, "|", fixed = TRUE)[[1]]
    cat(sprintf("  SELECTION DIFFERS %s/%s\n    only main: %s\n    only alex: %s\n",
                parts[1], parts[2],
                paste(utils::head(gateS[[k]]$only_main, 6), collapse = ","),
                paste(utils::head(gateS[[k]]$only_alex, 6), collapse = ",")))
    ## The path prefix that shares this table, e.g. S0/STgradient/value/elements/bravo
    selection_prefixes <- c(selection_prefixes,
                            paste0(parts[1], "/", parts[2], "/", sub("/columns/gene$", "", parts[3])))
    ## A step's own warnings/messages sit OUTSIDE the result-table prefix, but they
    ## are just as much a consequence of which genes were selected: fitting a
    ## different gene set produces a different number of degenerate loess/wilcox
    ## fits, hence a different warning count. Same cascade, different path depth.
    selection_steps <<- c(selection_steps, paste0(parts[1], "/", parts[2]))
  }
}

## ---- classification -------------------------------------------------------

cat("\n== classification ==\n")
verdicts <- list()
all_paths <- union(names(pass_a), names(pass_b))
matched_expect <- character(0)

for(p in all_paths){
  sc <- xv_scenario_of(p)
  in_a  <- p %in% names(pass_a)
  in_b  <- p %in% names(pass_b)
  in_bp <- p %in% names(pass_bp)

  cls <- NULL; rec <- NULL

  if(p %in% noise){
    cls <- "NOISE"
  } else if(sc %in% cascade_scenarios && !identical(sc, "S11")) {
    cls <- "INCOMPARABLE_CASCADE"
  } else if((length(selection_prefixes) &&
             any(startsWith(p, paste0(selection_prefixes, "/")))) ||
            (length(selection_steps) &&
             p %in% c(paste0(selection_steps, "/warnings"),
                      paste0(selection_steps, "/messages")))) {
    ## The two versions selected different genes for this table, so its rows are
    ## not the same rows. Not order-only, and not evidence of a regression either.
    cls <- "SELECTION_CASCADE"
  } else {
    rec <- xv_match_expect(p, sc, main_val = fa_raw[[p]], alex_val = fb_raw[[p]])
    ## EXPECTED_CHANGE is distinct from EXPECTED_FIX on purpose: a deliberate
    ## change that is not a bug fix (a stored intermediate dropped for capacity, a
    ## new column, a message rewritten) should not be reported as if main had been
    ## broken. Both are "accounted for"; only one is a claim about correctness.
    if(!is.null(rec) && rec$class %in% c("NON_COMPARABLE", "EXPECTED_FIX",
                                         "EXPECTED_CHANGE",
                                         "UNDOCUMENTED_CHANGE")){
      cls <- rec$class; matched_expect <- c(matched_expect, rec$id)
    } else if(in_bp && !is.null(rec) && identical(rec$class, "ORDER_INDUCED_VALUE")){
      cls <- "ORDER_INDUCED_VALUE"; matched_expect <- c(matched_expect, rec$id)
    } else if(in_a && !in_b){
      ## Order-only requires BOTH that the canonical pass is clean AND that this
      ## object's name sets matched (gate C). Absence from Pass B alone could mean a
      ## canonicalization rule is masking a real difference.
      suspect <- length(gateC_prefixes) &&
                 any(startsWith(p, paste0(gateC_prefixes, "/")))
      cls <- if(suspect) "REGRESSION" else "ORDER_ONLY"
    } else if(in_b && !in_bp){
      cls <- "FLOAT_NOISE"
    } else {
      cls <- "REGRESSION"
    }
  }
  verdicts[[p]] <- list(class = cls, expect = if(is.null(rec)) NA_character_ else rec$id,
                        a = in_a, b = in_b, bp = in_bp,
                        detail = if(in_b) pass_b[[p]] else pass_a[[p]])
}

tab <- table(vapply(verdicts, function(v) v$class, character(1)))
for(nm in names(tab)) cat(sprintf("  %-22s %d\n", nm, tab[[nm]]))

## Exhaustiveness: if the classifier cannot account for every observed path, say
## so rather than quietly dropping one.
stopifnot(length(verdicts) == length(all_paths))
cat(sprintf("  %-22s %d (accounted for: %d)\n", "TOTAL", length(all_paths), sum(tab)))

## ---- MISSING_EXPECTED -----------------------------------------------------

cat("\n== inverse check: expected differences that were NOT observed ==\n")
missing <- list()
for(rec in XV_EXPECT){
  if(rec$id %in% matched_expect) next
  scen_ran <- is.na(rec$scenario) || rec$scenario %in% common_sc
  missing[[rec$id]] <- list(rec = rec,
                            verdict = if(!scen_ran) "COVERAGE_GAP" else
                                      if(isTRUE(rec$optional)) "NOT_OBSERVED" else "MISSING_EXPECTED")
}
if(!length(missing)){
  cat("  ok  every expectation was observed\n")
} else {
  for(id in names(missing)){
    m <- missing[[id]]
    cat(sprintf("  %-16s %-24s scenario=%s commit=%s\n", m$verdict, id,
                if(is.na(m$rec$scenario)) "any" else m$rec$scenario,
                if(is.na(m$rec$commit)) "-" else m$rec$commit))
  }
}

## ---- vacuous coverage -----------------------------------------------------

## A step that errors identically on both versions compares only an error string
## and therefore tests no code, yet reads as "no difference". That is exactly how
## the pseudobulk CV fix stayed unverified: count it as coverage debt, not green.
cat("\n== vacuous coverage (error on BOTH versions) ==\n")
vac <- character(0); tot <- 0L
for(sc in common_sc){
  for(st in names(pm$scenarios[[sc]]$steps)){
    if(is.null(pa$scenarios[[sc]]$steps[[st]])) next
    tot <- tot + 1L
    om <- pm$scenarios[[sc]]$steps[[st]]$raw$ok
    oa <- pa$scenarios[[sc]]$steps[[st]]$raw$ok
    if(!isTRUE(om) && !isTRUE(oa)) vac <- c(vac, paste(sc, st, sep = "/"))
  }
}
cat(sprintf("  %d of %d compared steps compared only an error string\n", length(vac), tot))
for(v in vac) cat(sprintf("      %s\n", v))

## ---- persist --------------------------------------------------------------

saveRDS(list(pass0 = pass0, noise = noise, pass_a = pass_a, pass_b = pass_b,
             pass_bp = pass_bp, gateC = gateC, gateP = gateP,
             verdicts = verdicts, missing = missing, vacuous = vac,
             scenarios = common_sc,
             headers = list(main = pm$header, alex = pa$header)),
        file.path(OUT, sprintf("comparison-rep%s.rds", CMP_REP)), version = 3)

regressions <- names(Filter(function(v) identical(v$class, "REGRESSION"), verdicts))
cat(sprintf("\n== %d REGRESSION path(s) ==\n", length(regressions)))
for(p in utils::head(regressions, 60)){
  d <- verdicts[[p]]$detail
  cat(sprintf("  %s\n     kind=%s%s\n", p, d$kind,
              if(!is.null(d$max_rel)) sprintf(" max_rel=%.3g n=%d", d$max_rel, d$n)
              else if(!is.null(d$n)) sprintf(" n=%s", d$n) else ""))
}
if(length(regressions) > 60) cat(sprintf("  ... and %d more\n", length(regressions) - 60))
cat("\nwrote ", file.path(OUT, sprintf("comparison-rep%s.rds", CMP_REP)), "\n", sep = "")
