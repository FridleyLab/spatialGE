##
# @title import_Visium_h5: Reads Visium `space ranger` HDF5 outputs
# @description Reads the HDF5 file and spot coordinates in the output directory
# of a Visium run, and returns data for creation of a STList.
# @details
# The function takes as an argument the path a HDF5 file containing the gene counts
# and a coordinates .csv file from the `spatial` sub-directory. It returns a list
# with sparse count matrix (genes x spots) and a dataframe with spot coordinates.
# Zero-count genes across spots are removed.
#
# @param counts_fp Path to a Visium HDF5 (.h5) file within.
# @param coords_fp Path to the tissue_positions_list.csv file.
# @return x A list with two objects: Visium counts in `dgCMatrix` sparse format and
# a data frames with spot coordinates.
#
#
import_visium_h5 = function(counts_fp=NULL, coords_fp=NULL){
  # Read HDF5 file
  h5_read = hdf5r::H5File$new(filename=counts_fp, mode='r')

  # Warning in case there is more than one group in .h5 file or group is not called
  # `matrix`. 10X indicates their HDF5 outputs have one group named ('matrix')
  h5_group = names(h5_read)
  if(length(h5_group) > 1){
    warning('There is more than one group in \"', h5_file, '\". Reading only the first group.')
  }
  if(h5_group[1] != 'matrix'){
    warning('The .h5 file group is not called \"matrix\". Is this a Visium .h5 output?')
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
  spotcoords_df = data.table::fread(coords_fp, header=F, check.names=F) %>%
    dplyr::select(barcode=1, intissue=2, imagerow=5, imagecol=6) %>%
    dplyr::filter(intissue == 1) %>%
    dplyr::select(barcode, imagerow, imagecol) %>%
    tibble::as_tibble()

  # Error if number of spot coordinate do not match number of spots in counts
  if(nrow(spotcoords_df) != ncol(counts_mtx)){
    stop('Number of spots do not match between count and coordinate data.')
  }

  visium_list = list(rawcounts=counts_mtx, coords=spotcoords_df)
  return(visium_list)
}

