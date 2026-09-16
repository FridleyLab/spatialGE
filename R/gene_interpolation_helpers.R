# Helpers for gene_interpolation(): per-sample and per-gene units of work for
# kriging transformed gene expression.

##
# @title interp_one_sample: Krige all requested genes for one sample
# @description The per-sample unit of work dispatched by `gene_interpolation()`'s
# `sge_lapply` call, extracted so it is directly callable/benchmarkable for one
# sample without going through the whole STlist.
#
# @param tr_counts_mat the sample's transformed count matrix.
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param genes the gene names to interpolate for this sample.
# @param sample_name the sample's name, used in progress messages.
# @param verbose either logical or an integer (0, 1, or 2); values > 1 print
# a per-gene progress message.
# @return a named list (by gene) of kriging results.
#
interp_one_sample = function(tr_counts_mat, spatial_meta_df, genes, sample_name, verbose){
  kriging_res = list()
  for(gene in genes){
    kriging_res[[gene]] = interp_krige_one_gene(tr_counts_mat, spatial_meta_df, gene, sample_name, verbose)
  }
  return(kriging_res)
}


##
# @title interp_krige_one_gene: Krige one gene's expression for one sample
# @description The per-gene unit of work called by `interp_one_sample()`,
# extracted so it is directly callable/benchmarkable for a single gene.
#
# @param tr_counts_mat the sample's transformed count matrix.
# @param spatial_meta_df the sample's spatial metadata tibble.
# @param gene the gene name to interpolate.
# @param sample_name the sample's name, used in progress messages.
# @param verbose either logical or an integer (0, 1, or 2); values > 1 print
# a progress message.
# @return a list with `krige_out` (a data frame of predicted coordinates and
# values, or NULL on error) and `success_or_not`.
#
interp_krige_one_gene = function(tr_counts_mat, spatial_meta_df, gene, sample_name, verbose){
  # Get transformed counts
  gene_expr = data.frame(
    position=names(unlist(tr_counts_mat[rownames(tr_counts_mat) == gene, ])),
    expr_val=as.vector(unlist(tr_counts_mat[rownames(tr_counts_mat) == gene, ])))

  # Match order of library names in counts data and coordinate data
  gene_expr = gene_expr[match(spatial_meta_df[[1]], gene_expr[[1]]), ]
  gene_geo_df = dplyr::bind_cols(spatial_meta_df[c('xpos', 'ypos')], gene_expr=as.numeric(gene_expr[[2]]))
  gene_geo_df = as.data.frame(gene_geo_df)
  sp::coordinates(gene_geo_df) = ~xpos+ypos

  rm(gene_expr) # Clean env

  # Create a grid with points in the middle of each coordinate
  ny = range(sp::coordinates(gene_geo_df)[, 'ypos'])
  ny = length(seq(ny[1], ny[2], by=100))
  nx = range(sp::coordinates(gene_geo_df)[, 'xpos'])
  nx = length(seq(nx[1], nx[2], by=100))
  grid_sf = sf::st_make_grid(gene_geo_df, n=c(nx, ny), what="centers")

  # Give hope to users...
  if(verbose > 1L){
    system(sprintf('echo "%s"', paste0("\tInterpolating ", gene, ' in ', sample_name, "....")))
  }

  # Perform kriging
  suppressMessages({
    krige_out = tryCatch({
      dat_vgm = gstat::variogram(gene_expr~1, gene_geo_df) # Calculates sample variogram values
      est_vgm = gstat::fit.variogram(dat_vgm, gstat::vgm('Exp')) # Fit variaogram
      krige_pred = gstat::krige(gene_expr~1, gene_geo_df, grid_sf, model=est_vgm) # Predict values
    },
    error=function(err){return(err)}
    )
  })

  # Compute surface if parameters available and save kriging values
  # Record if kriging was successful
  #if(any(class(cov_est) == 'simpleError') | !cov_est[['optimSuccess']]){
  if(any(class(krige_out) == 'simpleError')){
    result = list(krige_out=NULL,
                  success_or_not='error')
  } else{
    krige_tmp = as.data.frame(sf::st_coordinates(krige_out))
    krige_tmp[['krige']] = krige_out[['var1.pred']]
    result = list(krige_out=krige_tmp)
    rm(krige_tmp) # Clean env
    if(any(class(krige_out) == 'simpleWarning')){
      result[['success_or_not']] = 'warning'
    } else{
      result[['success_or_not']] = 'kriging_completed'
    }
  }
  rm(krige_out) # Clean env

  return(result)
}
