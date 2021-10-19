##
# This function returns a list of plots with three QC stats: The total number of
# reads, number of genes, and percent of mtDNA genes per each spot.
#
#
plot_QCstats = function(x=NULL, qc='percent_mt', plot_who=NULL, saveplot=NULL){

  require('ggplot2')

  if(is.null(plot_who)){
    plot_who = 1:length(x@counts)
  }

  qcplot = list()
  for(i in plot_who){
    total_reads = colSums(x@counts[[i]][, -1])

    if(qc == 'total_reads'){
      df_current = dplyr::bind_cols('sample'=rep(names(x@counts)[i], ncol(x@counts[[i]][, -1])),
                                    'total_reads'=total_reads)
      qcplot[[names(x@counts)[i]]] = ggplot(df_current, aes(x=sample, y=total_reads)) +
        geom_violin() +
        geom_point(position=position_jitter(seed=1, width=0.2), size=0.1, alpha=0.5) +
        theme(axis.title.x=element_blank())
    }

    if(qc == 'total_genes'){
      nonzero_df = x@counts[[i]][, -1] != 0
      total_genes = colSums(nonzero_df)
      df_current = dplyr::bind_cols('sample'=rep(names(x@counts)[i], ncol(x@counts[[i]][, -1])),
                                    'total_genes'=total_genes)
      qcplot[[names(x@counts)[i]]] = ggplot(df_current, aes(x=sample, y=total_genes)) +
        geom_violin() +
        geom_point(position=position_jitter(seed=1, width=0.2), size=0.1, alpha=0.5) +
        theme(axis.title.x=element_blank())
    }

    if(qc == 'percent_mt'){
      total_mtct = colSums(x@counts[[i]][grep("^MT-", x@counts[[i]][[1]]), -1])
      df_current = dplyr::bind_cols('sample'=rep(names(x@counts)[i], ncol(x@counts[[i]][, -1])),
                                  'total_mt'=total_mtct,
                                  'percent_mt'=total_mtct/total_reads)
      qcplot[[names(x@counts)[i]]] = ggplot(df_current, aes(x=sample, y=percent_mt)) +
        geom_violin() +
        geom_point(position=position_jitter(seed=1, width=0.2), size=0.1, alpha=0.5) +
        theme(axis.title.x=element_blank())
    }

  }

  mfrow_n = n2mfrow(length(plot_who))
  qc_pgrid = ggpubr::ggarrange(plotlist=qcplot, ncol=mfrow_n[2], nrow=mfrow_n[1])

  if(is.null(saveplot)){
    print(qc_pgrid)
  } else{
    ggsave(saveplot, plot=qc_pgrid)
  }
}

