##
# @title import_xenium_mex: Reads Xenium outputs in MEX format (Xenium Ranger)
# @description Reads Market Exchange (MEX) format files in the output directory of a
# Xenium run, and returns data for creation of a STList.
# @details
# The function takes as an argument the paths to the features, barcodes, and matrix
# outputs of a Xenium run, as well as the spot coordinates in the `spatial` directory.
# Then, reads the files and returns a list with sparse count matrix (genes x spots) and
# a dataframe with spot coordinates. Zero-count genes across spots are removed.
#
# @param features_fp File path to the features.tsv.gz file.
# @param barcodes_fp File path to the barcodes.tsv.gz file.
# @param counts_fp File path to the matrix.mtx.gz file.
# @param coords_fp File path to the tissue_positions_list.csv file.
# @return a list with two objects: Xenium counts in `dgCMatrix` sparse format and
# a data frames with spot coordinates.
#
#
import_visium = function(features_fp=NULL, barcodes_fp=NULL, counts_fp=NULL, coords_fp=NULL){
  # Read in feature data
  features_df = data.table::fread(features_fp, header = F, check.names =F) %>%
    dplyr::rename("emsb" = 1,
                  "gene" = 2,
                  "dtype" = 3) %>%
    dplyr::select(-dtype) %>%
    dplyr::mutate(feat_n = as.character(seq(nrow(.)))) %>%
    dplyr::relocate(feat_n, .before = 1)

  # Read in barcode data
  barcodes_df = data.table::fread(barcodes_fp, header = F, check.names = F) %>%
    dplyr::mutate(spot_n = as.character(seq(nrow(.)))) %>%
    dplyr::relocate(spot_n, .before = 1) %>%
    dplyr::rename("barcode" = 2)

  # Read in coordinate data
  coords_df = data.table::fread(coords_fp, header=F, check.names=F) %>%
    dplyr::rename('barcode' = 1, 'intissue' = 2, 'array_row' = 3,
                  'array_col' = 4, 'imagerow' = 5, 'imagecol' = 6) %>%
    #dplyr::mutate(spotname = paste0("y", array_row, "x", array_col)) %>%
    dplyr::filter(intissue == 1) %>%
    dplyr::select(barcode, imagecol,imagerow)

  # Read in count data
  counts_df = data.table::fread(counts_fp, header=F, check.names = F, sep=" ", skip = 3) %>%
    dplyr::rename("feat_n" = 1, "spot_n" = 2, "counts" = 3)

  # Merge files together
  counts_all_df = dplyr::inner_join(coords_df, barcodes_df, by='barcode')
  counts_all_df <- dplyr::inner_join(counts_all_df %>% dplyr::mutate(spot_n = as.integer(spot_n)), counts_df, by='spot_n')
  counts_all_df <- dplyr::inner_join(counts_all_df, features_df %>% dplyr::mutate(feat_n = as.integer(feat_n)), by='feat_n')
  counts_all_df = counts_all_df %>%
    dplyr::mutate(spot_n = ifelse(is.na(spot_n), "otherBCs", spot_n),
                  #spotname = ifelse(is.na(spotname), "otherBCs", spotname),
                  emsb = ifelse(is.na(emsb), "noGene_", emsb),
                  counts = ifelse(is.na(counts), 0, counts)) %>%
    data.table::as.data.table()

  rawcounts_df = data.table::dcast.data.table(counts_all_df, emsb + gene ~ barcode, value.var = "counts", fill = 0) %>%
    dplyr::filter(emsb != "noGene_") %>%
    dplyr::select(-contains("otherBCs"), -emsb) %>%
    dplyr::mutate(gene=make.unique(gene)) %>%
    data.frame(check.names = F)

  spotcoords_df <- counts_all_df[, c('barcode', 'imagerow', 'imagecol', 'spot_n')] %>%
    dplyr::distinct(.keep_all = T) %>%
    dplyr::filter(spot_n != "otherBCs") %>%
    dplyr::select(-spot_n) %>%
    data.frame(check.names = F)

  #zeroSpots = colnames(rawcounts_df)[colSums(rawcounts_df) == 0]
  #rawcounts_df = rawcounts_df[, !(colnames(rawcounts_df) %in% zeroSpots)]

  # Remove zero-count genes
  rawcounts_df = rawcounts_df[rowSums(rawcounts_df[, -1]) > 0, ] %>%
    makeSparse(.)
  # Match column names (barcodes) between counts and coordinates
  spotcoords_df = spotcoords_df[(spotcoords_df$barcode %in% colnames(rawcounts_df)), ] %>%
    tibble::as_tibble()

  xenium_list = list(rawcounts=rawcounts_df, coords=spotcoords_df)
  return(xenium_list)
}

##
# @title import_xenium_h5: Reads Xenium HDF5 outputs (Xenium Ranger)
# @description Reads the HDF5 file and spot coordinates in the output directory
# of a Xeniu, run, and returns data for creation of a STList.
# @details
# The function takes as an argument the path a HDF5 file containing the gene counts
# and a coordinates (parquet) file. It returns a list with sparse count matrix
# (genes x spots) and a dataframe with spot coordinates. Zero-count genes across
# spots are removed.
#
# @param counts_fp Path to a Xenium HDF5 (.h5) file.
# @param coords_fp Path to the cells.parquet file.
# @return a list with two objects: Xenium counts in `dgCMatrix` sparse format and
# a data frames with spot coordinates.
#
#
import_xenium_h5 = function(counts_fp=NULL, coords_fp=NULL){
  # Read HDF5 file
  h5_read = hdf5r::H5File$new(filename=counts_fp, mode='r')

  # Warning in case there is more than one group in .h5 file or group is not called
  # `matrix`. 10X indicates their HDF5 outputs have one group named ('matrix')
  h5_group = names(h5_read)
  if(length(h5_group) > 1){
    warning('There is more than one group in \"', h5_file, '\". Reading only the first group.')
  }
  if(h5_group[1] != 'matrix'){
    warning('The .h5 file group is not called \"matrix\". Is this a Xenium .h5 output?')
  }
  rm(counts_fp)

  # Create sparse matrix from h5 file
  counts_mtx = Matrix::sparseMatrix(
    x=as.numeric(h5_read[[paste0(h5_group[1], '/data')]][]),
    i=h5_read[[paste0(h5_group[1], '/indices')]][],
    p=h5_read[[paste0(h5_group[1], '/indptr')]][],
    dims=h5_read[[paste0(h5_group[1], '/shape')]][],
    index1=F
  )
  rownames(counts_mtx) = make.unique(h5_read[[paste0(h5_group[1], '/features/name')]][])
  colnames(counts_mtx) = h5_read[[paste0(h5_group[1], '/barcodes/')]][]

  # Close file handle
  h5_read$close_all()
  rm(h5_read, h5_group)

  # Remove zero-count genes
  counts_mtx = counts_mtx[Matrix::rowSums(counts_mtx) > 0, ]

  # Read in coordinate data
  if(grepl('\\.parquet$', coords_fp)){
  spotcoords_df = arrow::read_parquet(coords_fp) %>%
    dplyr::select(barcode=cell_id, ypos=y_centroid, xpos=x_centroid) %>% ### FUTURE DEV: OTHER METADATA IS AVAILABLE... SHOULD IT BE IMPORTED TOO?
    tibble::as_tibble()
  } else if(grepl('\\.csv$', coords_fp)){
    spotcoords_df = data.table::fread(coords_fp, header=F, check.names=F)
    if(all(!(c('cell_id', 'y_centroid', 'x_centroid') %in% colnames(spotcoords_df)))){
      raise_err(err_code='error0040')
    }
    spotcoords_df = spotcoords_df %>%
      dplyr::select(barcode=cell_id, ypos=y_centroid, xpos=x_centroid) %>% ### FUTURE DEV: OTHER METADATA IS AVAILABLE... SHOULD IT BE IMPORTED TOO?
      tibble::as_tibble()
  } else{
    raise_err(err_code='error0039')
  }

  # Warning if number of spot coordinate do not match number of spots in counts
  # Subset to spots with coordinate information (to possibly deal with raw_feature_bc_matrix.h5)
  if(nrow(spotcoords_df) != ncol(counts_mtx)){
    warning('Number of cells do not match between count and coordinate data. Subsetting to cells with coordinates.')
    counts_mtx = counts_mtx[ , colnames(counts_mtx) %in% spotcoords_df[['barcode']] ]

    # Check that spots are available after subsetting
    if(ncol(counts_mtx) == 0){
      stop('Cells in expression data and coordinate data do not match.')
    }
  }

  xenium_list = list(rawcounts=counts_mtx, coords=spotcoords_df)
  return(xenium_list)
}

