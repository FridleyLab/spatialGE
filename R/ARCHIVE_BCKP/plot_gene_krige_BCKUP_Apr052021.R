##
# This function produces a kriging plot for a series of HUGO gene names and subjects.
#
# @param x, an STList with kriging objects for the genes selected.
# @param genes, a vector of gene names (one or several) to plot.
# @param plot_who, a vector of subject indexes as ordered within the STList, to
# plot genes from. If NULL, will plot for all subjects.
# @color_pal, a scheme from 'khroma'.
# @param saveplot, a file path where quilt plots will be saved. If NULL, plots
# are printed to console
#
#
# require('geoR')
# require('RColorBrewer')
plot_gene_krige <-
  function(x=NULL, genes=NULL, krige_type='ord', plot_who=NULL,
           color_pal='sunset', saveplot=NULL){

    # Test that a gene name was entered.
    if (is.null(genes)) {
      stop("Please, enter one or more genes to plot.")
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
      for (gene in genes) {

          if(length(x@gene_krige[[gene]]) > i){
            if(is.null(x@gene_krige[[gene]][[i]])){
              cat(paste0(gene, " kriging for subject ", i, " is not prsent in STList\n"))
              next
            }
          }else{
            cat(paste0(gene, " kriging for subject ", i, " is not prsent in STList\n"))
            next
          }

        assign("color_palGb", color_pal, envir=.GlobalEnv)

        moran_est <- round(as.vector(x@gene_het[[gene]][[i]]$morans_I$estimate[[1]]), 2)
        geary_est <- round(as.vector(x@gene_het[[gene]][[i]]$gearys_C$estimate[[1]]), 2)
        getis_est <- round(as.vector(x@gene_het[[gene]][[i]]$getis_ord_Gi$estimate[[1]]), 4)

        if(krige_type == 'ord'){
          title <- paste0(gene, ", subj ", i, " - ordinary kriging voom-norm counts\nMorans I=",
                          moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                          getis_est)
          assign('stlist', x@gene_krige[[gene]][[i]]$ord, envir=.GlobalEnv)
        }
        else if(krige_type == 'univ'){
          title <- paste0(gene, ", subj ", i, " - universal kriging voom-norm counts\nMorans I=",
                          moran_est, "  Gearys C=", geary_est, "  GetisOrd Gi=",
                          getis_est)
          assign('stlist', x@gene_krige[[gene]][[i]]$univ, envir=.GlobalEnv)
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
