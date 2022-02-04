##
#' @title STclust: Detect clusters on ST arrays
#' @description Perform unsupervised spatially-informed clustering on the spots of an
#' spatial array.
#' @details
#' The function takes an STList and calculates euclidean distances between spots based on the x,y
#' spatial locations and expression of the top variable genes (`Seurat::FindVariableFeatures`).
#' Weighted averages between the two distances are calculated to perform hierarchical clustering. The
#' user can define how much weight the spatial distance should have (values
#' between 0.05 to 0.25 work reasonably well). The method aims to detect compartments
#' within a tissue.
#'
#' @param x an STList with normalized expression data.
#' @param weights a double (0-1) indicating the weight to be applied to spatial distances.
#' Defaults to w=0.025.
#' @param dist, the distrance measurement to be used. Defaults to 'euclidean'. Other
#' options as provided by `wordspace::dist.matrix`.
#' @param linkage the linkage method applied to hierarchical clustering. Passed to
#' `hclust` and defaults to 'ward.D'.
#' @param ks the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied.
#' @param topgenes the number of genes with highest spot-to-spot expression variation. The
#' variance is calculated via `Seurat::FindVariableFeatures`.
#' @param deepSplit, a logical or integer (1-4), to be passed to `cutreeDynamic` and
#' control cluster resolution.
#' @return x, the STList with cluster assignments.
#'
#' @examples
#' # In this example, melanoma is an STList.
#' # Using Dynamic Tree Cuts:
#' # melanoma <- STclust(melanoma, ks='dtc', weights=0.1)
#'
#' # Using a range of ks:
#' # melanoma <- STclust(melanoma, ks=c(2:6), weights=0.1)
#'
#' @export
#
#
STclust = function(x=NULL, weights=0.025, dist='euclidean', linkage='ward.D', ks='dtc', topgenes=2000, deepSplit=F) {

  require('magrittr')

  # Clustering method to use. Set because other methods will be supported in future versions
  clmethod = 'hclust'

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STList')){
    stop("The input must be a STList.")
  }

  # Check to ensure number of ks is acceptable
  if(is.numeric(ks)){
    if(any(ks < 2)){
      stop('Refusing to generate < 2 clusters.')
    }
  }

  grp_list = list()
  for(i in 1:length(x@tr_counts)){

    # Find tops variable genes using Seurat approach. In the past, instead of Seurat, genes with the highest stdev were used
    if(any(colnames(x@gene_var[[i]]) == 'vst.variance.standardized')){
      x@gene_var[[i]] = x@gene_var[[i]][, !grepl('vst.variance.standardized', colnames(x@gene_var[[i]]))]
    }
    x@gene_var[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
      tibble::rownames_to_column(var='gene') %>%
      dplyr::select('gene', 'vst.variance.standardized') %>%
      dplyr::left_join(x@gene_var[[i]], ., by='gene')
    topgenenames = x@gene_var[[i]][order(x@gene_var[[i]]$vst.variance.standardized, decreasing=T), ]
    topgenenames = topgenenames[['gene']][1:topgenes]

    # Get transformed counts
    trcounts_df = x@tr_counts[[i]][x@tr_counts[[i]]$gene %in% topgenenames, ]

    # Convert counts and coordinate data to matrices
    A = t(as.matrix(trcounts_df[, -1]))
    rownames(A) = colnames(trcounts_df[, -1])
    B = as.matrix(x@coords[[i]][,-1])
    rownames(B) = x@coords[[i]][[1]]

    # Get distance matrices
    dA <- wordspace::dist.matrix(A, method=dist)
    dB <- dist(B, upper=T, diag=T)
    dAm <- as.matrix(dA)
    dBm <- as.matrix(dB)

    # Scale matrices
    dAm <- dAm/max(dAm)
    dBm <- dBm/max(dBm)

    for(w in 1:length(weights)){
      weight_d <- weights[w]
      weight_g <- 1-weight_d

      # Create vector of weights for Reduce
      weightList <- c(weight_g, weight_d)
      dMxs <- list(dAm, dBm)

      # Apply wieight element-wise
      weight_m <- Reduce(`+`,Map(`*`, dMxs, weightList))
      # Convert result to matrix
      weight_m <- matrix(weight_m, nrow(dAm))
      weight_m_d <- as.dist(weight_m)

      if(clmethod == 'hclust'){
        # Performe hclust
        hierclusters <- hclust(weight_m_d, method=linkage)

        # Apply dtc or split to k
        if(is.character(ks)){
          if(ks == 'dtc'){
            x@misc[['STclust_cuttype']] <- 'dtc'
            grp_df <- dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weight_m, deepSplit=deepSplit, verbose=F)
            grp_df <- tibble::tibble(colnames(dAm), as.factor(grp_df))
            names(grp_df) <- c('libname', 'WCluster')

            grp_df$WCluster[grp_df$WCluster == 0] <- NA

            grp_df <- dplyr::full_join(x@coords[[i]], grp_df, by='libname')

            grp_list[[paste0(names(x@tr_counts[i]), "_spw", weights[w])]] <- grp_df
          }
        }else if(is.numeric(ks)){
          x@misc[['STclust_cuttype']] <- 'multiK'
          for(k in ks){
            singlek <-  cutree(hierclusters, k=k)
            singlek <- tibble::tibble(colnames(dAm), as.factor(singlek))
            names(singlek) <- c('libname', 'WCluster')
            singlek <- dplyr::full_join(x@coords[[i]], singlek, by='libname')
            grp_list[[paste0(names(x@tr_counts[i]), "_spw", weights[w], '_k', k)]] <- singlek
          }
        } else{
          stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
        }
      }
      else{
        clusterlouvain = bluster::clusterRows(weight_m,
                                              BLUSPARAM=bluster::NNGraphParam(
                                                shared=T,
                                                cluster.fun="louvain", k=nn))
        x@misc[['STclust_cuttype']] <- 'louvain'
        grp_df <- tibble::tibble(colnames(dAm), as.factor(clusterlouvain))
        names(grp_df) <- c('libname', 'WCluster')

        #grp_df$WCluster[grp_df$WCluster == 0] <- NA

        grp_df <- dplyr::full_join(x@coords[[i]], grp_df, by='libname')

        grp_list[[paste0(names(x@tr_counts[i]), "_spw", weights[w])]] <- grp_df
      }
    }
  }
  x@st_clusters <- grp_list

  return(x)
}
