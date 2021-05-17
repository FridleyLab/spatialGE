##
#' @title cluster_STspot
#' @description Perform {INSERT CLUSTER ALGORITHM HERE :|} of spots within ST array.
#' @details
#' yadda yadda yadda
#'
#' @param x, an STList.
#' @param who, an integer indicating the spatial array to be analyzed.
#' @export
#
#
cluster_STspot <- function(x=NULL, who=NULL) {

  require('magrittr')
  require('ggplot2')
  require('NMF')
  require('fastICA')

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # TEMPORARY VARIABLE... Intended to be part of a loop later.
  who=1
  method='brunet'
  seed=123465
  #seed ='ica'

  # if(!is.null(clinvar)){
  #   clinvar_vals <- as.character(x@clinical[[clinvar]])
  #   clinvar_vals <- tibble::as_tibble_col(clinvar_vals, column_name='var_vals')
  # }
  #
  # # Create data frame to store "bulk counts".
  # bulkcounts_df <- x@counts[[1]][, 1]
  #
  # # Loop through count matrices in STList
  # for(i in 1:length(x@counts)){
  #
  #   bulk_expr <- rowSums(x@counts[[i]][, -1])
  #   bulk_expr <- tibble::as_tibble_col(bulk_expr,
  #                                      column_name=paste0("st_", i))
  #   bulk_expr <- bulk_expr %>%
  #     tibble::add_column(x@counts[[i]][, 1], .before=1)
  #
  #   bulkcounts_df <- dplyr::inner_join(bulkcounts_df,
  #                                      bulk_expr, by=c('gene'))
  #   }
  #
  # norm_factors <- edgeR::calcNormFactors(bulkcounts_df[, -1], method='TMM',
  #                                        lib.size=NULL)
  #
  # # Create new object to store the size-corrected counts.
  # bulkvoom_df <- c()
  #
  # # Divide each count value by their respective column (sample) normalization
  # # factor.
  # for(raw_col in names(bulkcounts_df[, -1])){
  #   bulkvoom_df <- dplyr::bind_cols(bulkvoom_df, tibble::as_tibble(
  #     bulkcounts_df[raw_col] / norm_factors[raw_col]))
  # }
  #
  # bulkvoom <- limma::voom(bulkcounts_df[, -1], design=NULL,
  #                         lib.size=colSums(bulkcounts_df[, -1]),
  #                         normalize.method='none', plot=F)
  #
  # bulkvoom_df <- tibble::as_tibble(bulkvoom$E, .name_repair='minimal')
  # bulkvoom_df <- bulkvoom_df %>%
  #   tibble::add_column(bulkcounts_df[, 1], .before=1)

  #bulkvoom_df <- x@voom_counts[[who]]
  bulkvoom_df <- x@counts[[who]]

  spotlibs <- c('gene')
  for(voomcol in 2:ncol(bulkvoom_df[, -1])){
    if(sum(x@counts[[who]][, voomcol] > 0) > 500){
      spotlibs <- append(spotlibs, names(bulkvoom_df)[voomcol])
    }
  }

  bulkvoom_df <- x@voom_counts[[who]]

  bulkvoom_subset_df <- bulkvoom_df[, names(bulkvoom_df) %in% spotlibs]

  variable_genes <- x@gene_stdev[[who]][order(x@gene_stdev[[who]]$gene_stdevs, decreasing=T), 1][1:500, ]

  bulkvoom_subset_df <- bulkvoom_subset_df[unlist(bulkvoom_subset_df$gene) %in% unlist(variable_genes), ]

  bulkvoom_mx <- t(as.matrix(bulkvoom_subset_df[, -1]))

  ranks <- 5

  res <- nmf(bulkvoom_mx, 2:ranks, nrun=30, seed=seed, .options='v', method=method)


 # BIC <- mclust::mclustBIC(bulkvoom_mx, G=c(1:20))

  #mod1 <- mclust::Mclust(bulkvoom_mx, x=BIC)


  #factoextra::fviz_mclust(mod1, "BIC", palette = "jco")

clustplots <- list()

for(rank in 2:ranks){

W_mx <- res$fit[[as.character(rank)]]@fit@W
colnames(W_mx) <- paste0('r', rep(1:rank))

hard_rank <- colnames(W_mx)[max.col(W_mx,ties.method="first")]
names(hard_rank) <- rownames(W_mx)


  clusters <- x@coords[[who]][x@coords[[who]]$X1 %in% names(hard_rank), ]

  clusters <- dplyr::bind_cols(clusters,
                               tibble::as_tibble_col(hard_rank, column_name='cluster'))



  p1 <-   ggplot(clusters) +
    geom_point(data=x@coords[[who]], aes(x=X2, y=X3), shape=1) +
    geom_point(aes(x=X2, y=X3, col=as.factor(cluster))) +
    labs(col='Ranks') +
    theme_classic()


clustplots[[rank]] <- p1

}

pdf('~/Desktop/new_nmf_clusters_Subj3.pdf', width=9, height=5.5)
gridExtra::grid.arrange(clustplots[[2]], clustplots[[3]], clustplots[[4]], clustplots[[5]])
dev.off()

pdf('~/Desktop/nmf_rankEval_Subj3.pdf')
plot(res)
dev.off()


}


