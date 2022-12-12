##
# @title plot_spatial_geneset: Gene set expression at each spot from an spatial array
# @description Plot averaged transformed expression of gene within a gene set (quilt plot).
# @details
# This function produces a quilt plot for a series of gene sets and spatial
# arrays within an STList.
#
# @param x an STList
# @param genes a vector of one or more gene names to plot.
# @param plot_who a vector of numbers indicating the spatial arrays to plot
# genes from. Numbers follow the order of `names(x@counts)`. If NULL, will plot
# all spatial arrays.
# @param color_pal, a color scheme from 'khroma' or RColorBrewer.
# @param visium logical, whether or not the samples are from a Visium experiment.
# @param ptsize a number specifying the size of the points. Passed to `size` aesthetic.
# @return a list with plots.
#
# @examples
# # In this example, melanoma is an STList.
# # qplots <- plot_gene_quilt(melanoma, genes=gene_sets, samples=2, visium=F)
#
# @export
#
# @importFrom methods as is new
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
plot_spatial_geneset = function(x=NULL, genes=NULL, samples=NULL, color_pal='BuRd', visium=T, ptsize=NULL){
  # May implement later
  image = F

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Define which samples to plot
  if(is.null(samples)){
    samples = names(x@spatial_meta)
  } else{
    if(is.numeric(samples)){
      samples = names(x@spatial_meta)[samples]
    }
    if(length(grep(paste0(samples, collapse='|'), names(x@spatial_meta))) == 0){
      stop('The requested samples are not present in the STList spatial metadata.')
    }
  }

  # Test if no transformed counts are available
  if(rlang::is_empty(x@tr_counts)){
    stop('There are no transformed counts in the STlist.')
  }

  avg_expr = list()
  # Loop through samples
  for(i in samples){
    counts = expandSparse(x@tr_counts[[i]])
    # Loop through gene sets
    avg_expr[[i]] = list()
    for(set in names(genes)){
      if(any(rownames(counts) %in% genes[[set]])){
        # Get averaged expression
        avg_expr[[i]][[set]] = as.data.frame(t(data.frame(colMeans(counts[rownames(counts) %in% genes[[set]], ]))))
        rownames(avg_expr[[i]][[set]]) = set
      }
    }
    # Make data frame from lots of averaged expression
    avg_expr[[i]] = dplyr::bind_rows(avg_expr[[i]])

    # Merge coordinate data
    avg_expr[[i]] = t(avg_expr[[i]]) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var='libname') %>%
      dplyr::left_join(x@spatial_meta[[i]] %>%
                         dplyr::select(libname, xpos, ypos), ., by='libname')
  }

  # Store maximum expression value to standardize color legend.
  maxvalue <- c()
  minvalue <- c()
  for(i in names(avg_expr)){
    for(set in names(genes)){
      # Test if gene name exists in normalized count matrix.
      if(any(colnames(avg_expr[[i]]) == set)){
        # Find maximum expression value for each spatial array.
        values_tmp <- unlist(avg_expr[[i]][[set]])
        maxvalue <- append(maxvalue, max(values_tmp))
        minvalue <- append(minvalue, min(values_tmp))
        rm(values_tmp) # Clean environment
      }
    }
  }
  # Find maximum value among selected spatial arrays.
  maxvalue <- max(maxvalue)
  minvalue <- min(minvalue)

  # Define size of points
  if(is.null(ptsize)){
    ptsize = 0.5
  }

  # Create list of plots.
  qp_list <- list()
  # Loop through each normalized count matrix.
  for (i in names(avg_expr)) {
    # Loop though genes to plot.
    for (set in colnames(avg_expr[[i]] %>% dplyr::select(-libname, -xpos, -ypos))){
      # Extract relevant data frame
      df_tmp = avg_expr[[i]] %>%
        dplyr::select(xpos, ypos, values:=!!set)

      # The color palette function in khroma is created by quilt_p() function.
      qp <- quilt_p(data_f=df_tmp, leg_name='averaged\ngene_expr', color_pal=color_pal,
                    title_name=paste0(set, "\n", "sample: ", i),
                    minvalue=minvalue, maxvalue=maxvalue, visium=visium, ptsize=ptsize)

      # Append plot to list.
      qp_list[[paste0(set, "_", i)]] = qp

      rm(df_tmp) # Clean environment
    }

    if(image && !is.null(x@misc[['sp_images']][[i]])){
      img_obj = grid::rasterGrob(x@misc[['sp_images']][[i]])
      qp_list[[paste0('image', names(counts[i]))]] = ggplot() +
        annotation_custom(img_obj)
    }
  }

  return(qp_list)
}

