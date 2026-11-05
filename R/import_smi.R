##
# @title import_smi: Prepare CosMx-SMI inputs to be formated into an STlist
#
# @param counts_fp the path to a CosMx *exprMat* file
# @param coords_fp the path to a CosMx *metadata* file
# @return x a list of processed counts and coordinates
#
#
import_smi = function(counts_fp=NULL, coords_fp=NULL, slidename=NULL){
  # Read gene counts
  counts_read = data.table::fread(counts_fp)

  # Separate FOVs and create sparse matrices
  # FILTER CELLS WITH ID ZERO (0)
  fov_ls = list()
  for(fovid in unique(counts_read[['fov']])){
    fov_ls[[paste0(slidename, '_fov_', fovid)]] = counts_read[counts_read[['fov']] == fovid, ] %>%
      dplyr::filter(cell_ID != 0) %>%
      dplyr::mutate(libname=paste0('fov_', fov, '_', 'cell_',cell_ID)) %>%
      #dplyr::select(-c('fov', 'cell_ID', 'cell')) %>% # NANOSTRING BRAIN CORTEX DATA SET HAD BOTH 'cell_ID' and 'cell'
      dplyr::select(!dplyr::matches('^fov$|^cell_ID$|^cell$')) %>% # NANOSTRING BRAIN CORTEX DATA SET HAD BOTH 'cell_ID' and 'cell'
      dplyr::relocate(libname, .before=1) %>%
      tibble::column_to_rownames('libname') %>%
      t() %>% as.data.frame() %>% tibble::rownames_to_column('gene_name')

    # Remove zero-count genes
    fov_ls[[paste0(slidename, '_fov_', fovid)]]  = fov_ls[[paste0(slidename, '_fov_', fovid)]][rowSums(fov_ls[[paste0(slidename, '_fov_', fovid)]][, -1]) > 0, ]
  }
  rm(counts_read) # Clean env

  # Read in coordinate data
  spotcoords_df = data.table::fread(coords_fp)

  # Separate FOV coordinates
  fov_coord_ls = list()
  for(fovid in unique(spotcoords_df[['fov']])){
    fov_coord_ls[[paste0(slidename, '_fov_', fovid)]] = spotcoords_df[spotcoords_df[['fov']] == fovid, ] %>%
      dplyr::select(c('fov', 'cell_ID', 'CenterX_local_px', 'CenterY_local_px')) %>%
      dplyr::mutate(libname=paste0('fov_', fov, '_', 'cell_', cell_ID)) %>%
      dplyr::select(c('libname', 'CenterX_local_px', 'CenterY_local_px')) %>%
      dplyr::rename(ypos=CenterY_local_px, xpos=CenterX_local_px) %>%
      tibble::as_tibble()

    # Error if number of spot coordinate do not match number of spots in counts
    if(nrow(fov_coord_ls[[paste0(slidename, '_fov_', fovid)]]) != ncol(fov_ls[[paste0(slidename, '_fov_', fovid)]][, -1])){
      stop(paste0('Number of cells in FOV ', fovid,  ' do not match between count and coordinate data.'))
    }

  }
  rm(spotcoords_df) # Clean env

  smi_list = list(rawcounts=fov_ls, coords=fov_coord_ls)
  return(smi_list)
}
