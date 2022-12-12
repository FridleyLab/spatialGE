##
# @title plot_gene_expression: Gene expression at each spot from an spatial array
# @description Plot raw or transformed gene expression levels at each spot (quilt plot)
# within a spatially-resolved transcriptomic array.
# @details
# This function produces a quilt plot for a series of gene names and spatial
# arrays within an STList. The function can also plot tumor/stroma classifications
# from ESTIMATE deconvolution results.
#
# @param x, an STList
# @param genes, a vector of one or more gene names to plot.
# @param plot_who, a vector of numbers indicating the spatial arrays to plot
# genes from. Numbers follow the order of `names(x@counts)`. If NULL, will plot
# all spatial arrays.
# @param color_pal, a color scheme from 'khroma' or RColorBrewer.
# @param data_type, one of 'tr' or 'raw', to plot transformed or raw counts
# respectively.
# @param image logical, whether to print the image stored for the spatial arrays
# @param visium, logical, whether or not the samples are from a Visium experiment.
# @param ptsize, a number specifying the size of the points. Passed to `size` aesthetic.
# @return a list with plots.
#
# @examples
# # In this example, melanoma is an STList.
# # qplots <- plot_gene_quilt(melanoma, genes='CD74', samples=2, visium=F)
#
# @export
#
# @importFrom methods as is new
#' @import ggplot2
#' @importFrom magrittr %>%
#
#
plot_spatial_expression = function(x=NULL, genes=NULL, samples=NULL, color_pal='BuRd', data_type='tr', image=F, visium=T, ptsize=NULL){

  # Test that a gene name was entered.
  if (is.null(genes)) {
    stop("Please, enter one or more genes to plot.")
  }

  # Remove duplicated genes entered by user.
  genes = unique(genes)

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

  # Select appropriate slot to take counts from
  if(data_type == 'tr'){
    counts = x@tr_counts[samples]
  }else if(data_type == 'raw'){
    counts = x@counts[samples]
    # Expand sparse matrices
    # for(i in samples){
    #   counts[[i]] = expandSparse(counts[[i]])
    #   counts[[i]] = tibble::as_tibble(tibble::rownames_to_column(counts[[i]], var='gene'))
    # }
  } else(
    stop('Please, select one of transformed (tr) or raw (raw) counts')
  )

  # Extract genes and trasnpose expression data (to get genes in columns)
  # Also join coordinate data
  for(i in samples){
    if(any(rownames(counts[[i]]) %in% genes)){
      counts[[i]] = counts[[i]] %>%
        expandSparse() %>%
        tibble::rownames_to_column(var='gene') %>%
        dplyr::filter(gene %in% genes) %>%
        tibble::column_to_rownames(var='gene') %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var='libname') %>%
        dplyr::left_join(x@spatial_meta[[i]] %>%
                           dplyr::select(libname, xpos, ypos), ., by='libname')
    } else{
      # Remove samples for which none of the requested genes are available
      counts = counts[grep(i, names(counts), value=T, invert=T)]
      warning(paste0('Not all requested genes are available for sample ', i, '.'))
    }
  }

  # Test if requested data slot is available.
  if(rlang::is_empty(counts)) {
    stop("Data was not found in the specified slot of this STList.")
  }

  # Store maximum expression value to standardize color legend.
  maxvalue <- c()
  minvalue <- c()
  for(i in names(counts)){
    for(gene in genes){
      # Test if gene name exists in normalized count matrix.
      if(any(colnames(counts[[i]]) == gene)){
        # Find maximum expression value for each spatial array.
        values_tmp <- unlist(counts[[i]][[gene]])
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
  for (i in names(counts)) {
    # Loop though genes to plot.
    for (gene in colnames(counts[[i]] %>% dplyr::select(-libname, -xpos, -ypos))){
      # Extract relevant data frame
      df_tmp = counts[[i]] %>%
        #dplyr::rename(values:=!!gene)
        dplyr::select(xpos, ypos, values:=!!gene)

      # Set legend title
      if(data_type == 'tr'){
        qlegname = paste0(x@misc[['transform']], "\nexpr")
      } else{
        qlegname = 'raw\nexpr'
      }

      # The color palette function in khroma is created by quilt_p() function.
      qp <- quilt_p(data_f=df_tmp, leg_name=qlegname, color_pal=color_pal,
                    title_name=paste0(gene, "\n", "sample: ", i),
                    minvalue=minvalue, maxvalue=maxvalue, visium=visium, ptsize=ptsize)

      # Append plot to list.
      qp_list[[paste0(gene, "_", i)]] = qp

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

