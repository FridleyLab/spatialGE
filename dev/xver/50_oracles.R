## Directional oracles.
##
## The comparison establishes THAT two versions differ. That is not the claim this
## exercise needs to support. The claim is that main was WRONG and alex-refactor is
## RIGHT, and a diff cannot say which side is which.
##
## So each oracle recomputes the correct answer independently of both packages,
## straight from the input data, and asserts which version matches it. These are
## falsifiable predictions: if main's output matched the correct formula, or alex's
## did not, the oracle would say so.
##
## Measured results:
##
##   delta 9 (e7917bc, pseudobulk_heatmap coefficient of variation)
##     @misc$scaled_pbulk_mtx is IDENTICAL on both versions, so the oracle is valid
##     rather than comparing two different inputs. Given that matrix:
##       alex's displayed genes == top-k by sd/mean  (the true CV)      -> TRUE
##       main's displayed genes == top-k by mean/sd  (the inverted one) -> TRUE
##       overlap between the two selections: 0 of 5 (fixture), 0 of 10 (Thrane)
##     So the two versions select entirely disjoint gene sets, and main's set is
##     exactly the least-variable genes. Every 1.2.x heatmap showed the opposite of
##     what its documentation promises.
##
##   delta 1 (f278bd0, the silent no-op coordinate sort)
##     main: xpos is NOT non-decreasing in any sample
##     alex: xpos IS non-decreasing in all samples
##     So the sort genuinely did nothing before and genuinely applies now.
##
## Usage: Rscript --vanilla dev/xver/50_oracles.R

pm1 <- readRDS(file.path(Sys.getenv("XVER_OUT"),"payload-main-rep1.rds"))
pa1 <- readRDS(file.path(Sys.getenv("XVER_OUT"),"payload-alex-rep1.rds"))
pm9 <- readRDS(file.path(Sys.getenv("XVER_OUT"),"payload-main-rep9.rds"))
pa9 <- readRDS(file.path(Sys.getenv("XVER_OUT"),"payload-alex-rep9.rds"))

cat("=== ORACLE delta 9: pseudobulk_heatmap CV direction ===\n")
oracle <- function(pm, pa, sc, step, mtxstep, k){
  gm <- pm$scenarios[[sc]]$steps[[step]]$raw$value$heatmaps[[1]]$genes_display_order
  ga <- pa$scenarios[[sc]]$steps[[step]]$raw$value$heatmaps[[1]]$genes_display_order
  sm <- pm$scenarios[[sc]]$steps[[mtxstep]]$raw$value
  sa <- pa$scenarios[[sc]]$steps[[mtxstep]]$raw$value
  # scaled_pbulk_mtx payload: sparse or dense matrix reduction
  getm <- function(p){
    if(!is.null(p$values)) matrix(p$values, nrow=p$dim[1], dimnames=list(p$rownames,p$colnames))
    else { m <- Matrix::sparseMatrix(i=p$i, p=p$p, x=p$x, dims=p$dim, index1=FALSE,
             dimnames=list(p$rownames,p$colnames)); as.matrix(m) }
  }
  Mm <- getm(sm); Ma <- getm(sa)
  cat(sprintf("  scaled_pbulk_mtx identical across versions: %s (dim %s)\n",
              isTRUE(all.equal(Mm, Ma, tolerance=0)), paste(dim(Mm), collapse="x")))
  # genes are columns of scaled_pbulk_mtx (samples x genes); heatmap transposes
  X <- if(nrow(Mm) < ncol(Mm)) t(Mm) else Mm   # genes as rows
  cv_true <- apply(X, 1, function(v) sd(v)/ (mean(v)+1e-20))
  cv_inv  <- apply(X, 1, function(v) (mean(v)+1e-20)/sd(v))
  top_true <- names(sort(cv_true, decreasing=TRUE))[1:k]
  top_inv  <- names(sort(cv_inv,  decreasing=TRUE))[1:k]
  cat(sprintf("  alex genes == top-%d by TRUE CV (sd/mean): %s\n", k, setequal(ga, top_true)))
  cat(sprintf("  main genes == top-%d by INVERTED  (mean/sd): %s\n", k, setequal(gm, top_inv)))
  cat(sprintf("  overlap between the two selections: %d of %d\n", length(intersect(gm,ga)), k))
}
oracle(pm1, pa1, "S4",  "s4_heatmap",  "s4_scaled_mtx", 5)
cat("\n  -- same oracle on real Thrane data --\n")
oracle(pm9, pa9, "S12", "s12_heatmap", "s12_scaled_mtx", 10)

cat("\n=== ORACLE delta 1: coords sorted by the third coords column? ===\n")
for(nm in c("S0")){
  sm <- pm1$scenarios[[nm]]$steps$STlist$raw$value$spatial_meta
  sa <- pa1$scenarios[[nm]]$steps$STlist$raw$value$spatial_meta
  for(s in names(sm)){
    xm <- sm[[s]]$columns$xpos; xa <- sa[[s]]$columns$xpos
    cat(sprintf("  %-8s main xpos non-decreasing: %-5s | alex: %-5s\n", s,
        !is.unsorted(xm), !is.unsorted(xa)))
  }
}
