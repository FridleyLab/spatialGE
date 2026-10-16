#!/usr/bin/env Rscript
#
# Static argument check for the package vignettes.
#
# Three of the four vignettes are listed in .Rbuildignore, so R CMD check only
# ever builds `basic_functions_vignette.Rmd`. The other three can therefore rot
# silently: a chunk that calls an exported function with an argument name that no
# longer exists is never executed by anything in CI, and nothing parses it either.
#
# This script extracts every R chunk from each .Rmd by hand, parses the combined
# code (catching syntax errors), and walks the AST looking for calls to exported
# spatialGE functions whose named arguments are not in that function's formals.
#
# Deliberately allowed: R performs *partial* matching of named arguments, so
# `compare_SThet(gene=)` resolves to the formal `genes=` and is legal -- it runs
# without error as long as only one formal starts with "gene". Such names are not
# reported. (They are still worth avoiding in documentation, but they are not
# broken.) What this catches is a name that matches no formal at all, e.g. after
# an argument is renamed or removed.
#
# Chunk extraction scans for the fences rather than using knitr::purl(), because
# purl() *evaluates* chunk options: these vignettes use `eval=net_check` (a
# variable defined only when the network is available), and under purl() those
# chunks error out so their code never reaches the output -- which silently hides
# exactly the code this check exists to inspect.
#
# Nothing here is evaluated, so it needs no network and no data. It can run on
# every push.
#
# Usage: Rscript dev/ci/check_vignette_args.R   (from the package root)
# Exits non-zero if any problem is found.
#
# Usage: Rscript dev/ci/check_vignette_args.R   (from the package root)
# Exits non-zero if any problem is found.

suppressMessages(pkgload::load_all(".", compile = FALSE, quiet = TRUE))

exported <- getNamespaceExports("spatialGE")
ns <- asNamespace("spatialGE")

vignettes <- list.files("vignettes", pattern = "[.]Rmd$", full.names = TRUE)
if (length(vignettes) == 0L) {
  stop("No vignettes found; run this from the package root.")
}

problems <- character()

record <- function(file, msg) {
  problems <<- c(problems, sprintf("%s: %s", basename(file), msg))
}

# Resolve a call head to an exported spatialGE function name, if it is one.
head_name <- function(hd) {
  if (is.name(hd)) return(as.character(hd))
  # pkg::fn / pkg:::fn
  if (is.call(hd) && length(hd) == 3L &&
      as.character(hd[[1]]) %in% c("::", ":::") &&
      identical(as.character(hd[[2]]), "spatialGE")) {
    return(as.character(hd[[3]]))
  }
  NULL
}

# Recursively inspect an expression for calls to exported functions.
walk <- function(e, file) {
  if (!is.call(e)) return(invisible(NULL))

  nm <- head_name(e[[1]])
  if (!is.null(nm) && nzchar(nm) && nm %in% exported) {
    fn <- tryCatch(get(nm, envir = ns), error = function(...) NULL)
    if (is.function(fn)) {
      fmls <- names(formals(fn))
      has_dots <- "..." %in% fmls
      # names(e)[-1] is "" for positional arguments; only check named ones.
      supplied <- names(e)[-1]
      supplied <- supplied[!is.na(supplied) & nzchar(supplied)]
      for (a in supplied) {
        if (a %in% fmls) next
        # R allows unambiguous partial matching of argument names.
        if (length(fmls[startsWith(fmls, a)]) > 0L) next
        if (has_dots) next
        record(file, sprintf(
          "%s(): no argument named '%s'. Formals are: %s",
          nm, a, paste(fmls, collapse = ", ")))
      }
    }
  }

  for (i in seq_along(e)) {
    if (!is.null(e[[i]])) walk(e[[i]], file)
  }
  invisible(NULL)
}

# Return the R source of every ```{r ...} chunk in an .Rmd, as one string.
extract_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  n <- length(lines)
  out <- character()
  i <- 1L
  while (i <= n) {
    if (!grepl("^\\s*```+\\s*\\{r", lines[i])) {
      i <- i + 1L
      next
    }
    # A chunk header may span lines; it ends at the first line closing with '}'.
    j <- i
    while (j <= n && !grepl("\\}\\s*$", lines[j])) j <- j + 1L
    if (j > n) break
    code_start <- j + 1L
    k <- code_start
    while (k <= n && !grepl("^\\s*```+\\s*$", lines[k])) k <- k + 1L
    if (k > code_start) out <- c(out, lines[code_start:(k - 1L)])
    i <- k + 1L
  }
  paste(out, collapse = "\n")
}

for (v in vignettes) {
  code <- extract_chunks(v)
  if (!nzchar(trimws(code))) next

  exprs <- tryCatch(parse(text = code), error = function(err) {
    record(v, paste("contains R that does not parse:", conditionMessage(err)))
    NULL
  })
  if (is.null(exprs)) next

  for (e in exprs) walk(e, v)
}

if (length(problems) > 0L) {
  cat("Vignette argument check FAILED:\n\n")
  cat(paste0("  - ", problems, collapse = "\n"), "\n", sep = "")
  cat("\nFix the vignette (or update the check) before releasing.\n")
  quit(status = 1L)
}

cat(sprintf("Vignette argument check passed (%d vignettes, no unknown arguments).\n",
            length(vignettes)))
