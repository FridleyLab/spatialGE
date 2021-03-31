require('geoR')
require('RColorBrewer')
plot_cell_krige <- function(x=NULL, cell=NULL, krige_type='ord'){

  moran_est <- round(as.vector(x@cell_het[[cell]]$morans_I$estimate[[1]]), 2)
  geary_est <- round(as.vector(x@cell_het[[cell]]$gearys_C$estimate[[1]]), 2)
  getis_est <- round(as.vector(x@cell_het[[cell]]$getis_ord_Gi$estimate[[1]]), 2)

  if(krige_type == 'ord'){
    title <- paste0(cell, " - ordinary kriging ", x@cell_deconv$deconv_method,
                    " deconvolution\nMorans I=", moran_est,
                    "  Gearys C=", geary_est, "  GetisOrd Gi=", getis_est)
    assign('stlist', x@cell_krige[[cell]]$ord, envir=.GlobalEnv)
  }
  else if(krige_type == 'univ'){
    title <- paste0(cell, " - universal kriging ", x@cell_deconv$deconv_method,
                    " deconvolution\nMorans I=", moran_est,
                    "  Gearys C=", geary_est, "  GetisOrd Gi=", getis_est)
    assign('stlist', x@cell_krige[[cell]]$univ, envir=.GlobalEnv)
  }

  assign("titlekrige", title, envir = .GlobalEnv)

  if((max(stlist$predict) - (max(stlist$predict) * 0.1)) >= 0.001){
    image(stlist, col=rev(brewer.pal(9, "Spectral")),
          xlab="X Position", ylab="Y Position", main=titlekrige,
          y.leg=c(min(x@coords[,3]), (min(x@coords[,3]) + 1)),
          x.leg=c((max(x@coords[,2])-5), max(x@coords[,2])),
          scale.vals=c(round(min(stlist$predict) + (min(stlist$predict) * 0.1), 3),
                       round(max(stlist$predict) - (max(stlist$predict) * 0.1), 3))
    )
  }else{
    image(stlist, col=rev(brewer.pal(9, "Spectral")),
          xlab="X Position", ylab="Y Position", main=titlekrige,
          y.leg=c(min(x@coords[,3]), (min(x@coords[,3]) + 1)),
          x.leg=c((max(x@coords[,2])-5), max(x@coords[,2])))
  }

  rm(titlekrige, envir=.GlobalEnv)
  rm(stlist, envir=.GlobalEnv)
}
