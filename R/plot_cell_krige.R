##
# This function produces a kriging plot for a series of cell names and subjects.
#
# @param x, an STList with kriging objects for the cells selected.
# @param genes, a vector of cell names (one or several) to plot.
# @param plot_who, a vector of subject indexes as ordered within the STList, to
# plot cells from. If NULL, will plot for all subjects.
# @color_pal, a scheme from 'khroma'.
# @param saveplot, a file path where quilt plots will be saved. If NULL, plots
# are printed to console
#
#
# require('RColorBrewer')
require('geoR')
plot_cell_krige <-
  function(x=NULL, cells=NULL, krige_type='ord', plot_who=NULL,
           color_pal='sunset', saveplot=NULL){

    # Test that a cell name was entered.
    if (is.null(cells)) {
      stop("Please, enter one or more cell names to plot.")
    }

    # Test if no specific subject plot was requested.
    if (is.null(plot_who)) {
      plot_who <- c(1:length(x@counts))
    }

    # Prepare color palette
    col_fn <- khroma::colour(color_pal)
    color_pal <- as.vector(col_fn(9))

    # Loop through each of ythe subjects.
    for (i in plot_who) {

      # Loop though genes to plot.
      for (cell in cells) {

        assign("color_palGb", color_pal, envir=.GlobalEnv)

        if(cell != 'stroma_score'){
          moran_est <- round(as.vector(x@cell_het[[cell]][[i]]$morans_I$estimate[[1]]), 2)
          geary_est <- round(as.vector(x@cell_het[[cell]][[i]]$gearys_C$estimate[[1]]), 2)
          getis_est <- round(as.vector(x@cell_het[[cell]][[i]]$getis_ord_Gi$estimate[[1]]), 4)
        }

        if(krige_type == 'ord'){
          if(cell != 'stroma_score'){
            title <- paste0(cell, " - ordinary kriging ", x@cell_deconv$deconv_method,
                            " deconvolution\nMorans I=", moran_est,
                            "  Gearys C=", geary_est, "  GetisOrd Gi=", getis_est)
          }else{
            title <- paste0(cell, " - ordinary kriging ", x@cell_deconv$deconv_method,
                            " deconvolution")
          }
          assign('stlist', x@cell_krige[[cell]][[i]]$ord, envir=.GlobalEnv)
        }
        else if(krige_type == 'univ'){
          if(cell != 'stroma_score'){
            title <- paste0(cell, " - universal kriging ", x@cell_deconv$deconv_method,
                            " deconvolution\nMorans I=", moran_est,
                            "  Gearys C=", geary_est, "  GetisOrd Gi=", getis_est)
          }else{
            title <- paste0(cell, " - universal kriging ", x@cell_deconv$deconv_method,
                            " deconvolution")
          }
          assign('stlist', x@cell_krige[[cell]][[i]]$univ, envir=.GlobalEnv)
        }

        assign("titlekrige", title, envir=.GlobalEnv)

        ymin <- (min(x@coords[[i]][,3]) - 4)
        ymax <- (max(x@coords[[i]][,3]) + 1)

        assign('yminGb', ymin, envir=.GlobalEnv)
        assign('ymaxGb', ymax, envir=.GlobalEnv)

        image(stlist, col=color_palGb,
              locations=x@prediction_grid[[i]],
              borders=x@prediction_border[[i]],
              ylim=c(yminGb, ymaxGb),
              xlab="X Position", ylab="Y Position", main=titlekrige,
              y.leg=c((min(x@coords[[i]][,3]) - 3), (min(x@coords[[i]][,3]) - 2)),
              x.leg=c((min(x@coords[[i]][,2])), (max(x@coords[[i]][,2])))
              #          scale.vals=c(round(min(stlist$predict) + (min(stlist$predict) * 0.05), 2),
              #                     round(max(stlist$predict) - (max(stlist$predict) * 0.05), 2))
        )

        # if((max(stlist$predict) - (max(stlist$predict) * 0.1)) >= 0.001){
        #   image(stlist, col=rev(brewer.pal(9, "Spectral")),
        #         xlab="X Position", ylab="Y Position", main=titlekrige,
        #         y.leg=c(min(x@coords[,3]), (min(x@coords[,3]) + 1)),
        #         x.leg=c((max(x@coords[,2])-5), max(x@coords[,2])),
        #         scale.vals=c(round(min(stlist$predict) + (min(stlist$predict) * 0.1), 3),
        #                      round(max(stlist$predict) - (max(stlist$predict) * 0.1), 3))
        #   )
        # }else{
        #   image(stlist, col=rev(brewer.pal(9, "Spectral")),
        #         xlab="X Position", ylab="Y Position", main=titlekrige,
        #         y.leg=c(min(x@coords[,3]), (min(x@coords[,3]) + 1)),
        #         x.leg=c((max(x@coords[,2])-5), max(x@coords[,2])))
        # }

        rm(list=c(
          'color_palGb',
          'titlekrige',
          'stlist',
          'yminGb',
          'ymaxGb'),
          envir=.GlobalEnv)
      }
    }
  }
