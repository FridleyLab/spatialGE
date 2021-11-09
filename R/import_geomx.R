##
# @title importGeomx: Creates an STList from GeoMx-DSP outputs
# @description Reads .dcc and .pkc files, as well as sample annotation files with
# spatial coordinates and returns an STList for downstream analysis with spatialGE.
# @details
# The function takes as an argument the path to a directory containig .dcc files,
# a path to a .pkc file, and a sample metadata file containing spatial or clinical
# information for each ROI. The arguments `slide_col`, `id_col`, `x_col`, and `y_col` indicate the
# column names in the metadata file from which slide IDs, ROI IDs, y coordinates, and x coordinates
# should be taken.
#
# @param dcc Path to directory containing `.dcc` files.
# @param pkc File path to the `.pkc` file.
# @param annots File path to comma- or tab-delimited, or an Excel spreadsheet with
# coordinate (x and y) and clinical metadata, with one row per ROI
# @param slide_col, String with name of column in `annots` with IDs for each slide
# @param id_col, String with name of column in `annots` with ROI IDs.
# @param x_col, String with name of column in `annots` with x coordinates
# @param y_col, String with name of column in `annots` with y coordinates
# @return x, an STList with the GeoMx counts and coordinates.
#
#
import_Geomx = function(dcc=NULL, pkc=NULL, annots=NULL, slide_col=NULL, id_col=NULL, x_col=NULL, y_col=NULL){
  # Get filepaths for each DCC file
  dcc_fps = list.files(dcc, full.names=T, pattern='.dcc$', recursive=T)
  # Create list with un-processed text for each DCC
  dcc_list = lapply(dcc_fps, function(i){
    counts = readLines(i)
    startLine = grep('<Code_Summary>', counts)
    endLine = grep('</Code_Summary>', counts)
    roiName = grep('^ID,', counts, value=T)
    roiName = gsub('ID,', '', roiName)
    counts = counts[-c(endLine:length(counts))]
    counts = append(roiName, counts[-c(1:startLine)])
  })

  # Process text for each DCC, and create data frame with probe name and counts
  roiName = c()
  for(i in 1:length(dcc_list)){
    roiName = append(roiName, dcc_list[[i]][1])
    dcc_list[[i]] = dcc_list[[i]][-1]
    dcc_list[[i]] = tibble::as.tibble(dcc_list[[i]])
    dcc_list[[i]] = tidyr::separate(dcc_list[[i]], 1, into=c('probe', roiName[i]), sep=',')
    dcc_list[[i]][[2]] = as.numeric(dcc_list[[i]][[2]])
  }
  names(dcc_list) = roiName

  # Read metadata file and get coordinate information
  if(grepl('.xls', annots)){
    annots_df = readxl::read_excel(annots)
  } else{
    delim = readLines(annots, n=1)
    if(grepl(',', delim)){
      delim = ','
    } else if(grepl('\t', delim)){
      delim = '\t'
    } else{
      stop('Annotation file is not comma- or tab-delimited.')
    }
    annots_df = readr::read_delim(annots, delim=delim, col_types=readr::cols())
  }
  # Remove rows without slide ID
  annots_df = annots_df[!is.na(annots_df[[slide_col]]), ]

  # Extract probe ID and paired gene names
  pkcs = rjson::fromJSON(file=pkc)
  pkcs_df = tibble::tibble()
  for(i in 1:length(pkcs$Targets)){
    pkcs_df = dplyr::bind_rows(pkcs_df,
                               c(probe=pkcs$Targets[[i]]$Probes[[1]]$RTS_ID, gene=pkcs$Targets[[i]]$DisplayName))
                               #c(probe=pkcs$Targets[[i]]$Probes[[1]]$RTS_ID, gene=pkcs$Targets[[i]]$Probes[[1]]$DisplayName))
  }

  # Organize ROI counts within data frames for each slide
  return_list = list()
  return_list[['counts']] = list()
  return_list[['coords']] = list()
  for(slide in unique(annots_df[[slide_col]])){
    slide_annots_df = annots_df[annots_df[[slide_col]] == slide, ]
    return_list[['coords']][[slide]] = slide_annots_df[, c(id_col, y_col, x_col)]
    colnames(return_list[['coords']][[slide]]) = c('spot', 'ypos', 'xpos')

    slide_dcc_list = dcc_list[unique(slide_annots_df[[id_col]])]

    if(length(slide_dcc_list) == 1){
      return_list[['counts']][[slide]] = slide_dcc_list[[1]]
      return_list[['counts']][[slide]] = dplyr::left_join(pkcs_df, return_list[['counts']][[slide]],  by='probe')
    } else{
      return_list[['counts']][[slide]] = slide_dcc_list[[1]]
      for(roi in 2:length(slide_dcc_list)){
        return_list[['counts']][[slide]] = dplyr::full_join(return_list[['counts']][[slide]], slide_dcc_list[[roi]], by='probe')
      }
      return_list[['counts']][[slide]] = dplyr::left_join(pkcs_df, return_list[['counts']][[slide]], by='probe')
    }
  }

  # Process duplicates genes using RTS probe IDs
  dup_rts = c()
  for(i in 1:length(names(return_list[['counts']]))){
    dup_genes_Mask = duplicated(return_list[['counts']][[i]]$gene)
    dup_rts = append(dup_rts, return_list[['counts']][[i]]$probe[dup_genes_Mask])
  }
  dup_rts = unique(dup_rts)
  for(i in 1:length(names(return_list[['counts']]))){
    # for(rts in 1:length(dup_rts)){
    #   if(length(dup_rts[rts]) !=0){
    #     geneToChange = return_list[['counts']][[i]]$gene[return_list[['counts']][[i]]$probe == dup_rts[rts]]
    #     if(!is.null(geneToChange)){
    #       return_list[['counts']][[i]]$gene[return_list[['counts']][[i]]$probe == dup_rts[rts]] = paste0(geneToChange, '_d', rts)
    #     }
    #   }
    # }
    return_list[['counts']][[i]] = return_list[['counts']][[i]][, -1]
  }

  return(return_list)
}
