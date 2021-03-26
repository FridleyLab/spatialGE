# This function performs spatial interpolation of normalized gene counts. This
# function takes a STList and a list of gene names or the token 'top' for the 10
# genes with the highest variation. It also calculates spatial heterogeneity
# measures for the genes. The function can perform ordinary or universal kriging.
# The result can be plotted using the plot_gene_krige() function.
#
# @param x, a STList with normalized counts
# @plot_who, a vector of HUGO names or 'top'. If 'top', kriging for the 10 genes
# with highest standard deviation is estimated.
# @univ, a logical stating whether or not to perform universal or ordinary kriging.
# @res, a number to adjust the resolution of the plot. Fractions of 1 lead to
# more resolution.
# @return x, a STList including an spatial interpolation object.
#
#
# require('concaveman')
# require('geoR')
# require('RColorBrewer')

gene_krige <- function(x=NULL, genes='top', univ=F, res=0.1, who=NULL){

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # If genes='top', get names of 10 genes with the highest standard deviation.
  if(length(genes) == 1){
    if(genes == 'top'){
      genes <- x@gene_stdev[[i]]$gene[order(x@gene_stdev[[i]]$gene_stdevs, decreasing=T)][1:10]
    }
  }

  # Test if no specific subject plot was requested.
  if (is.null(who)) {
    who <- c(1:length(x@counts))
  }

  # Test if voom normalized counts are available.
  if (is_empty(stlists@voom_counts)) {
    stop(paste("There are not normalized matrices in this STList."))
  }

  # Loop through each normalized count matrix.
  for (i in who) {

    # Loop through genes.
    for(gene in genes){

      # Test that a normalized matrix for a gene is present.
      if(!any(x@voom_counts[[i]]$gene == gene)){
        cat(paste(gene, "is not a gene in the normalized count matrix from
                  subject", i, "."))
        next
      }

      # Test if slot for gene kriging is already present. Else, create it.
      if(!is.null(x@gene_krige[[gene]])){
        if(length(x@gene_krige[[gene]]) < i){
          x@gene_krige[[gene]][[i]] <- list(ord=NULL,
                                            univ=NULL)
        }
      }else{
        x@gene_krige[[gene]] <- list()
        x@gene_krige[[gene]][[i]] <- list(ord=NULL,
                                          univ=NULL)
      }

      # Extract expression data for a goven gene.
      gene_expr <- x@voom_counts[[i]][x@voom_counts[[i]]$gene == gene, -1]

      # Transpose expression data to turn it into a column. Then turn library
      # names into a column and assign column names (first row).
      gene_expr <- as.data.frame(t(gene_expr))
      gene_expr <- gene_expr %>% rownames_to_column(., var='position')
      #colnames(gene_expr) <- gene_expr[1,]
      colnames(gene_expr)[2] <- 'gene_expr'
      #gene_expr <- as.data.frame(gene_expr[-1,])

      # Sort expression data using the order in the mapping file. Then add
      # coordinates data to expression data frame.
      gene_expr <- gene_expr[match(x@coords[[i]][[1]], gene_expr[[1]]), ]
      gene_geo_df <- cbind(x@coords[[i]][2:3], as.numeric(gene_expr[[2]]))

      # Create concave hull to use as delimiter of sampled area. Needs to be
      # done before converting data frame to spatial object.
      conc_hull <- concaveman(as.matrix(gene_geo_df[1:2]))

      # Create geodata object from expression and coordinate data
      gene_geo <- as.geodata(gene_geo_df, coords.col=c(1,2), data.col=3)

      # Create a grid finer than the sampled locations to predict locations.
      gene_geo_grid <-expand.grid(
        seq((min(x@coords[[i]][[2]])-1), (max(x@coords[[i]][[2]])+1), by=res),
        seq((min(x@coords[[i]][[3]])-1), (max(x@coords[[i]][[3]])+1), by=res)
      )

      # Add concave hull to geodata.
      gene_geo$borders <- conc_hull

      # OC <- output.control(simulations=TRUE, n.pred=10,
      #                      quantile=c(0.1, 0.25, 0.5, 0.75, 0.9),
      #                      threshold = 350)

      # Create controls for either ordinary or universal kriging, and perform
      # estimation.
      if(univ == F){
        # NOTE: Need to check how to decide on init.cov.pars
        gene_geo_lhood <- likfit(gene_geo, trend='cte', ini.cov.pars=c(1, 0.15))

        # Specify control (and output) parameters for ordinary kriging.
        KC <- krige.control(obj.model=gene_geo_lhood)

        # Perform ordinary kriging.
        gene_krig <- krige.conv(gene_geo, locations=gene_geo_grid, krige=KC)

        x@gene_krige[[gene]][[i]][['ord']] <- gene_krig

      }else if(univ == T){
        # NOTE: Need to use regression analysis of variogram to get values for
        # nugget.
        gene_geo_lhood <- likfit(gene_geo, trend='cte',
                                 ini.cov.pars=c(1000, 500), nug=100)

        # Specify control (and output) parameters for universal kriging.
        KC <- krige.control(type.krige="OK", obj.m=gene_geo_lhood,
                            trend.d="cte",
                            trend.l="cte")

        # Perform universal kriging.
        gene_krig <- krige.conv(gene_geo, locations=gene_geo_grid, krige=KC
                                #output = OC
        )

        x@gene_krige[[gene]][[i]][['univ']] <- gene_krig
      }

      # Calculate spatial heterogeneity statistics.
      x <- gene_moran_I(x, genes=gene, who=i)
      x <- gene_geary_C(x, genes=gene, who=i)
      x <- gene_getis_Gi(x, genes=gene, who=i)

    }

  }

  return(x)

}
