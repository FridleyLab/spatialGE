require('geoR')
require('RColorBrewer')
plot_gene_krige <- function(x=NULL, gene=NULL, krige_type='ord'){

  moran_est <- round(as.vector(x@gene_het[[gene]]$morans_I$estimate[[1]]), 2)
  geary_est <- round(as.vector(x@gene_het[[gene]]$gearys_C$estimate[[1]]), 2)
  getis_est <- round(as.vector(x@gene_het[[gene]]$getis_ord_Gi$estimate[[1]]), 2)

  if(krige_type == 'ord'){
    title <- paste0(gene, " - ordinary kriging voom-norm counts\nMorans I=",
                    moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                    getis_est)
    assign('stlist', x@gene_krige[[gene]]$ord, envir=.GlobalEnv)
  }
  else if(krige_type == 'univ'){
    title <- paste0(gene, " - universal kriging voom-norm counts\nMorans I=",
                    moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                    getis_est)
    assign('stlist', x@gene_krige[[gene]]$univ, envir=.GlobalEnv)
  }

  assign("titlekrige", title, envir=.GlobalEnv)

  image(stlist, col=rev(brewer.pal(9, "Spectral")),
        xlab="X Position", ylab="Y Position", main=titlekrige#,
#        y.leg=c(min(x@coords[,3]), (min(x@coords[,3]) + 1)),
#        x.leg=c((max(x@coords[,2])-1), (max(x@coords[,2])+3))#,
#        scale.vals=c(round(min(stlist$predict) + (min(stlist$predict) * 0.05), 2),
#                     round(max(stlist$predict) - (max(stlist$predict) * 0.05), 2))
  )

  rm(titlekrige, envir=.GlobalEnv)
  rm(stlist, envir=.GlobalEnv)
}
