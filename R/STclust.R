##
#' @title STclust: Detect clusters of spots/cells
#' @description Perform unsupervised spatially-informed clustering on the spots/cells of a
#' ST sample
#' @details
#' The function takes an STlist and calculates euclidean distances between cells or spots
#' based on the x,y spatial locations, and the expression of the top variable genes
#' (`Seurat::FindVariableFeatures`). The resulting distances are weighted by
#' applying 1-`ws` to the gene expression distances and `ws` to the spatial distances.
#' Hierarchical clustering is performed on the sum of the weighted distance matrices.
#' The `STclust` method allows for identification of tissue niches/domains that are
#' spatially cohesive.
#'
#' @param x an STlist with normalized expression data
#' @param ws a double (0-1) indicating the weight to be applied to spatial distances.
#' Defaults to 0.025
#' @param dist_metric the distance metric to be used. Defaults to 'euclidean'. Other
#' options are the same as in `wordspace::dist.matrix`
#' @param linkage the linkage method applied to hierarchical clustering. Passed to
#' `hclust` and defaults to 'ward.D'
#' @param ks the range of k values to assess. Defaults to `dtc`, meaning `cutreeDynamic`
#' is applied
#' @param topgenes the number of genes with highest spot-to-spot expression variation. The
#' variance is calculated via `Seurat::FindVariableFeatures`.
#' @param deepSplit a logical or integer (1-4), to be passed to `cutreeDynamic` and
#' control cluster resolution
#' @return an STlist with cluster assignments
#'
#' @examples
#' # Using included melanoma example (Thrane et al.)
#' library('spatialGE')
#' data_files <- list.files(system.file("extdata", package="spatialGE"), recursive=T, full.names=T)
#' count_files <- grep("counts", data_files, value=T)
#' coord_files <- grep("mapping", data_files, value=T)
#' clin_file <- grep("thrane_clinical", data_files, value=T)
#' melanoma <- STlist(rnacounts=count_files[c(1,3)], spotcoords=coord_files[c(1,3)], samples=clin_file)
#' melanoma <- transform_data(melanoma, method='log')
#' melanoma <- STclust(melanoma, ws=c(0, 0.025))
#' STplot(melanoma, ws=0.025, samples='ST_mel1_rep1', ptsize=1)
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom methods as is new
#' @importFrom stats as.dist complete.cases cutree dist hclust prcomp sd
#
#
STclust = function(x=NULL, ws=0.025, dist_metric='euclidean', linkage='ward.D2', ks='dtc', topgenes=2000, deepSplit=F){
  # Record time
  zero_t = Sys.time()
  verbose = 1L
  if(verbose > 0L){
    cat(crayon::green(paste0('STclust started.\n')))
  }

  # Clustering method to use. Set because other methods will be supported in future versions
  clmethod = 'hclust'

  ws = as.double(ws)
  topgenes = as.integer(topgenes)

  # Do not allow weights higher than 1
  if(any(ws < 0) | any(ws > 1)){
    stop('Please select a spatial weight between 0 and 1.')
  }

  # Test if an STList has been input.
  if(is.null(x) | !is(x, 'STlist')){
    stop("The input must be a STlist.")
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
    if(any(colnames(x@gene_meta[[i]]) == 'vst.variance.standardized')){
      x@gene_meta[[i]] = x@gene_meta[[i]][, !grepl('vst.variance.standardized', colnames(x@gene_meta[[i]]))]
    }

    #x@gene_meta[[i]] = Seurat::FindVariableFeatures(x@counts[[i]], verbose=F) %>%
    x@gene_meta[[i]] = Seurat_FindVariableFeatures(x@counts[[i]]) %>%
      tibble::rownames_to_column(var='gene') %>%
      dplyr::select('gene', 'vst.variance.standardized') %>%
      dplyr::left_join(x@gene_meta[[i]], ., by='gene')

    topgenenames = x@gene_meta[[i]] %>%
      dplyr::arrange(desc(vst.variance.standardized)) %>%
      dplyr::slice_head(n=topgenes) %>%
      dplyr::select(gene) %>%
      unlist() %>%
      as.vector()

    # Get transformed counts
    trcounts_df = x@tr_counts[[i]][rownames(x@tr_counts[[i]]) %in% topgenenames, ]

    # Convert counts and coordinate data to matrices
    A = t(as.matrix(trcounts_df))
    B = x@spatial_meta[[i]] %>%
      dplyr::select(libname, xpos, ypos) %>%
      tibble::column_to_rownames(var='libname')
    B = as.matrix(B[match(rownames(B), rownames(A)), ])

    # Get distance matrices
    dA = wordspace::dist.matrix(A, method=dist_metric)
    dB = dist(B, upper=T, diag=T)
    dAm = as.matrix(dA)
    dBm = as.matrix(dB)

    # Scale matrices
    dAm = dAm/max(dAm)
    dBm = dBm/max(dBm)

    for(w in 1:length(ws)){
      weight_d = ws[w]
      weight_g = 1-weight_d

      # Create vector of weights for Reduce
      weightList = c(weight_g, weight_d)
      dMxs = list(dAm, dBm)

      # Apply weight element-wise
      weight_m = Reduce('+', Map('*', dMxs, weightList))
      # Convert result to matrix
      weight_m_d = as.dist(weight_m)

      if(clmethod == 'hclust'){
        # Run hierarchical clustering
        hierclusters = hclust(weight_m_d, method=linkage)

        # Apply dtc or split to k
        if(is.character(ks)){
          if(ks == 'dtc'){
            # Set list element names based on weight and deepSplit
            if(is.logical(deepSplit)){
              if(deepSplit){
                dspl = 'True'
              } else{
                dspl = 'False'
              }
            } else{
              dspl = deepSplit
            }
            col_name = paste0('stclust_spw', ws[w], '_dspl', dspl)

            grp_df <- dynamicTreeCut::cutreeDynamic(hierclusters, method='hybrid', distM=weight_m, deepSplit=deepSplit, verbose=F)
            grp_df <- tibble::tibble(libname=colnames(dAm), !!col_name:=as.factor(grp_df))

            grp_df[[col_name]][grp_df[[col_name]] == 0] = NA
            if(any(colnames(x@spatial_meta[[i]]) == col_name)){
              x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
                dplyr::select(-!!col_name)
            }
            x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
              dplyr::full_join(., grp_df, by='libname')
          }
        }else if(is.numeric(ks)){
          for(k in ks){
            # Set list element names based on weight and k value
            col_name = paste0('stclust_spw', ws[w], '_k', k)
            singlek <- cutree(hierclusters, k=k)
            singlek <- tibble::tibble(libname=colnames(dAm), !!col_name:=as.factor(singlek))

            if(any(colnames(x@spatial_meta[[i]]) == col_name)){
              x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
                dplyr::select(-!!col_name)
            }
            x@spatial_meta[[i]] = x@spatial_meta[[i]] %>%
              dplyr::full_join(., singlek, by='libname')
          }
        } else{
          stop('Enter a valid number of k values to evaluate or \'dtc\' to apply cutreeDynamic.')
        }
      }
      else{
        stop('Currently, only spatially-informed hierarchical clustering is supported.')
      }
    }
  }

  # Print time
  end_t = difftime(Sys.time(), zero_t, units='min')
  if(verbose > 0L){
    cat(crayon::green(paste0('STclust completed in ', round(end_t, 2), ' min.\n')))
  }

  return(x)
}
