## Single-version order-sensitivity probe.
##
## The cross-version comparison finds that SThet(k=) and STgradient(distsumm='avg')
## differ materially between main and alex-refactor. That alone does not say WHY:
## a difference could be the ordering fix exposing pre-existing order sensitivity,
## or it could be a regression the refactor introduced.
##
## This distinguishes them without involving two versions at all. It loads ONE
## build, permutes the spot order inside an STlist, and re-runs. If the output moves
## by the same magnitude, the code path was always order-sensitive and the
## cross-version difference is a consequence of delta 1 (f278bd0) rather than a
## regression. If it does not move, the difference came from somewhere else and is a
## regression.
##
## Measured result on the fixture:
##   SThet(k=4)                     Moran's I moves up to 1.8e-2   <- order-sensitive
##   SThet() default distance listw Moran's I moves up to 1.5e-15  <- float noise only
##   STgradient(distsumm='avg')     Spearman r moves up to 4.3e-3  <- order-sensitive
##
## Usage: Rscript --vanilla dev/xver/30_order_probe.R

source(file.path(Sys.getenv("XVER_REPO"), "dev", "xver", "common.R"))

xv_load_spatialge(Sys.getenv("XVER_LIB_ALEX"), "2.0.0")
xv_source_driver(); xv_shim_fixture_image()
setwd(xv_scratch("orderprobe"))

base <- fixture_base_data(); fx <- fixture_write_list_dfs(base)
st <- STlist(rnacounts=fx$rnacounts, spotcoords=fx$spotcoords, cores=1, verbose=FALSE)
st <- transform_data(st, method="log", cores=1)

# Permute spots within each sample, keeping the same data.
perm <- function(x){
  set.seed(99)
  for(s in names(x@counts)){
    o <- sample.int(ncol(x@counts[[s]]))
    x@counts[[s]] <- x@counts[[s]][, o, drop=FALSE]
    if(length(x@tr_counts)) x@tr_counts[[s]] <- x@tr_counts[[s]][, o, drop=FALSE]
    d <- x@spatial_meta[[s]]
    d2 <- d[match(colnames(x@counts[[s]]), d$libname), , drop=FALSE]
    ## spatial_meta is a plain data.frame, so base row-subsetting leaves
    ## non-contiguous row names that column_to_rownames() later rejects -- the
    ## footgun dev/ARCHITECTURE.md documents, and the mechanism of the STenrich
    ## domain= bug. Reset them.
    rownames(d2) <- NULL
    x@spatial_meta[[s]] <- d2
  }
  x
}
stp <- perm(st)

cmp <- function(a, b, lab){
  ka <- sort(a$gene); kb <- sort(b$gene)
  cat(sprintf("  %-28s same gene set: %s\n", lab, identical(ka,kb)))
  m <- merge(a, b, by="gene", suffixes=c(".1",".2"))
  for(cn in setdiff(names(a), "gene")){
    c1 <- m[[paste0(cn,".1")]]; c2 <- m[[paste0(cn,".2")]]
    if(!is.numeric(c1)) next
    ok <- !is.na(c1)&!is.na(c2)
    if(!any(ok)) next
    cat(sprintf("      %-26s max_abs=%.3g\n", cn, max(abs(c1[ok]-c2[ok]))))
  }
}

cat("== SThet(k=4): unpermuted vs permuted, SAME version ==\n")
h1 <- SThet(st,  genes=c("Gene1","GradGene1","DomainGeneA"), method="moran", k=4, cores=1, verbose=FALSE)
h2 <- SThet(stp, genes=c("Gene1","GradGene1","DomainGeneA"), method="moran", k=4, cores=1, verbose=FALSE)
for(s in names(h1@gene_meta)){
  a <- as.data.frame(h1@gene_meta[[s]]); b <- as.data.frame(h2@gene_meta[[s]])
  a <- a[!is.na(a$moran_i), c("gene","moran_i")]; b <- b[!is.na(b$moran_i), c("gene","moran_i")]
  cmp(a, b, paste0("sample ", s))
}

cat("\n== SThet default (no k): unpermuted vs permuted ==\n")
d1 <- SThet(st,  genes=c("Gene1","GradGene1","DomainGeneA"), method="moran", cores=1, verbose=FALSE)
d2 <- SThet(stp, genes=c("Gene1","GradGene1","DomainGeneA"), method="moran", cores=1, verbose=FALSE)
for(s in names(d1@gene_meta)){
  a <- as.data.frame(d1@gene_meta[[s]]); b <- as.data.frame(d2@gene_meta[[s]])
  a <- a[!is.na(a$moran_i), c("gene","moran_i")]; b <- b[!is.na(b$moran_i), c("gene","moran_i")]
  cmp(a, b, paste0("sample ", s))
}

cat("\n== STgradient(distsumm='avg'): unpermuted vs permuted ==\n")
addann <- function(x){ for(s in names(x@spatial_meta)){ d <- x@spatial_meta[[s]]
  d[["xv_domain"]] <- ifelse(d[["xpos"]] <= 6, "A","B"); x@spatial_meta[[s]] <- d }; x }
g1 <- STgradient(addann(st),  annot="xv_domain", ref="A", distsumm="avg", topgenes=50, cores=1, verbose=FALSE)
g2 <- STgradient(addann(stp), annot="xv_domain", ref="A", distsumm="avg", topgenes=50, cores=1, verbose=FALSE)
for(s in names(g1)) cmp(g1[[s]][,c("gene","avg_spearman_r")], g2[[s]][,c("gene","avg_spearman_r")], paste0("sample ", s))
