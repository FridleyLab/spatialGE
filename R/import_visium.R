##
# @title importVisium: Reads Visium outputs and produce an STList.
# @description Reads files in the output folder of a Visium run, and returns an
# STList for downstream analysis with spatialGE.
# @details
# The function takes as an argument the path to an 'outs'  folder of a Visium run.
# It reads the data and converts it into an STList, which can be used in downstrean
# analysis with spatialGE. Optionally, the function can also output the count and
# coordinates files.
#
# @param features_fp File path to the features.tsv.gz file.
# @param barcodes_fp File path to the barcodes.tsv.gz file.
# @param counts_fp File path to the matrix.mtx.gz file.
# @param coords_fp File path to the tissue_positions_list.csv file.
# @param filterMT, logical, whether or not to filter mtDNA genes. Filters out gene names
# beginning with 'MT-'. Defaults to TRUE.
# @param savefiles, logical, whether or not to save the counts and coordinate files.
# The files are saved as 'visium_counts.txt' and visium_coords.txt'. Defaults to TRUE.
# @param stlist, logical, whether or not return an STList.
# @return x, an STList with the Visium counts and coordinates.
#
#
import_Visium <- function(features_fp=NULL, barcodes_fp=NULL, counts_fp=NULL, coords_fp=NULL, filterMT=F, savefiles=F, stlist=F){
  require(data.table)
  #features_df <- readr::read_delim(features_fp, delim="\t", col_names=F, col_types=readr::cols(), progress=F)
  features_df = fread(features_fp, header = F, check.names =F)
  names(features_df) <- c('emsb', 'gene', 'dtype')
  features_df <- tibble::add_column(features_df, feat_n=as.character(1:nrow(features_df)), .before=1)
  setkey(features_df, "feat_n")

  #barcodes_df <- readr::read_delim(barcodes_fp, delim="\t", col_names=F, col_types=readr::cols(), progress=F)
  barcodes_df = fread(barcodes_fp, header = F, check.names = F)
  barcodes_df <- tibble::add_column(barcodes_df, spot_n=as.character(1:nrow(barcodes_df)), .before=1)
  names(barcodes_df) <- c('spot_n', 'barcode')
  setkey(barcodes_df, "barcode")

  #coords_df <- readr::read_delim(coords_fp, delim=",", col_names=F, col_types=readr::cols(), progress=F)
  coords_df = fread(coords_fp, header=F, check.names=F)
  # spot_name <- c()
  # for(i in 1:nrow(coords_df)){
  #   spot_name <- append(spot_name, paste0('y', coords_df$V3[i], 'x', coords_df$V4[i]))
  # }
  # coords_df$spotname <- spot_name
  coords_df$spotname = paste0('y', coords_df[[3]], 'x', coords_df[[4]])
  names(coords_df) <- c('barcode', 'intissue', 'array_row', 'array_col', 'pxlcol', 'pxlrow', 'spotname')
  setkey(coords_df, "barcode")

  #counts_df <- readr::read_delim(counts_fp, delim=" ", col_names=F, skip=3, col_types="ccd", progress=F)
  counts_df = fread(counts_fp, header=F, check.names = F, sep=" ", skip = 3)
  names(counts_df) <- c('feat_n', 'spot_n', 'counts')
  setkey(counts_df, "spot_n")

  counts_all_df = coords_df[barcodes_df, on = "barcode"] %>% mutate(spot_n = as.integer(spot_n))
  setkey(counts_all_df, "spot_n")
  counts_all_df = counts_all_df[counts_df, on = "spot_n"]
  setkey(counts_all_df, "feat_n")
  counts_all_df = counts_all_df[features_df %>% mutate(feat_n = as.integer(feat_n)), on = "feat_n"]


  counts_all_df = dplyr::left_join(coords_df, barcodes_df, by='barcode')
  counts_all_df <- dplyr::left_join(counts_all_df %>% mutate(spot_n = as.integer(spot_n)), counts_df, by='spot_n')
  counts_all_df <- dplyr::full_join(counts_all_df, features_df %>% mutate(feat_n = as.integer(feat_n)), by='feat_n')

  counts_all_df$spot_n[which(is.na(counts_all_df$spot_n))] = 'otherBCs'
  counts_all_df$spotname[which(is.na(counts_all_df$spotname))] = 'otherBCs'
  counts_all_df$emsb[which(is.na(counts_all_df$emsb))] = 'noGene_'
  counts_all_df$counts[is.na(counts_all_df$counts)] = 0

  rawcounts_df <- tidyr::pivot_wider(counts_all_df, id_cols=c('emsb', 'gene'), names_from=spotname, values_from=counts, values_fill=0)
  rawcounts_df = rawcounts_df[!grepl('noGene_', rawcounts_df$emsb), ]
  #rawcounts_df <- rawcounts_df[, -1]
  rawcounts_df = rawcounts_df[, !grepl('otherBCs', colnames(rawcounts_df))]

  # keep_idx <- c()
  # for(gene in dup_genes){
  #   dup_idx <- grep(paste0("^", gene, "$"), rawcounts_df$gene)
  #
  #   dup_rows <- rowSums(rawcounts_df[dup_idx, -1])
  #   dup_rows_mask <- dup_rows >= max(dup_rows)
  #
  #   if(sum(dup_rows_mask) > 1){
  #       eq_dups <- which(dup_rows_mask)
  #       dup_rows_mask[eq_dups[2:length(eq_dups)]] <- FALSE
  #   }
  #
  #   keep_idx <- append(keep_idx, dup_idx[!dup_rows_mask])
  # }

  if(filterMT){
    #rawcounts_df <- rawcounts_df[-c(keep_idx), ]
    rawcounts_df <- rawcounts_df[!grepl("^MT-", rawcounts_df$gene), ]
  }

  spotcoords_df <- counts_all_df[, c('spotname', 'array_row', 'array_col')]
  spotcoords_df <- spotcoords_df[!duplicated(spotcoords_df$spotname), ]
  spotcoords_df = spotcoords_df[!grepl('otherBCs', spotcoords_df$spotname), ]

  zeroSpots = colnames(rawcounts_df[, -c(1, 2)])[colSums(rawcounts_df[, -c(1, 2)]) == 0]
  rawcounts_df = rawcounts_df[, !(colnames(rawcounts_df) %in% zeroSpots)]
  spotcoords_df = spotcoords_df[(spotcoords_df$spotname %in% colnames(rawcounts_df[, -c(1, 2)])), ]

  if(savefiles){
    write.table(rawcounts_df, file='./visium_counts.txt', sep="\t", row.names=F, quote=F)
    write.table(spotcoords_df, file='./visium_coords.txt', sep="\t", row.names=F, quote=F)
  }

  tmp_ctsdf <- tempfile(fileext = ".txt", pattern='spatialGEctsdf_')
  tmp_cdsdf <- tempfile(fileext = ".txt", pattern='spatialGEcdsdf_')
  write.table(rawcounts_df, file=tmp_ctsdf, sep="\t", row.names=F, quote=F)
  write.table(spotcoords_df, file=tmp_cdsdf, sep="\t", row.names=F, quote=F)

  tmp_ctsfp <- tempfile(fileext = ".txt", pattern='spatialGEctsfp_')
  tmp_cdsfp <- tempfile(fileext = ".txt", pattern='spatialGEcdsfp_')
  write(tmp_ctsdf, file=tmp_ctsfp)
  write(tmp_cdsdf, file=tmp_cdsfp)

  if(stlist){
    #TEMPORARY: To avoid pushing to Github and install over and over...
    source('~/OneDrive - Moffitt Cancer Center/SPATIAL_TRANSCRIPTOMICS/code/spatialGEdev/R/STList.R')
    #x <- STList(countfiles = tmp_ctsfp, coordfiles = tmp_cdsfp)
    #return(x)
  } else{
    visium_list = list(rawcounts=rawcounts_df, coords=spotcoords_df)
    return(visium_list)
  }

}
