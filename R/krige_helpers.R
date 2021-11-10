##
# @title krige_py: Spatial interpolation of gene expression using PyKrige
# @description Performs spatial interpolation ('kriging') of transformed gene
# counts in spatially-resolved transcriptomics data using PyKrige.
#
# @param gridx
# @param gridy
# @param geo_df
# @param univ
# @return predict
#
#
krige_py = function (gridx=NULL, gridy=NULL, geo_df=NULL, univ=NULL) {

  # Requires python 3.5+
  # Prompts Miniconda installation.
  if(!reticulate::py_module_available("pykrige")){
    reticulate::py_install(c("pykrige"))
  }

  # Importing the kriging functions as their respective names.
  UniversalKriging = reticulate::import("pykrige.uk")$UniversalKriging
  OrdinaryKriging = reticulate::import("pykrige.ok")$OrdinaryKriging
  np = reticulate::import("numpy")

  # Convert to Numpy arrays
  gridx = np$asarray(gridx)
  gridy = np$asarray(gridy)
  data = np$array(t(geo_df))

  if(univ == F){
    # Computing universal kriging on numpy array
    OK = OrdinaryKriging(
      data[1,],
      data[2,],
      data[3,],
      variogram_model="spherical",
      #variogram_model="gaussian",
      exact_values=FALSE
    )

    # Calculates a kriged grid and the associated variance.
    kriged = OK$execute("grid", gridx, gridy)

    # Convert to R vector.
    predict = as.vector(t(kriged[[1]]))

    return(predict)
  } else if (univ == T) {
    # Computing universal kriging on numpy array
    UK = UniversalKriging(
      data[1,],
      data[2,],
      data[3,],
      variogram_model="spherical",
      #drift_terms="regional_linear",
      exact_values=FALSE
    )

    # Calculates a kriged grid and the associated variance.
    kriged = UK$execute("grid", gridx, gridy)

    # Convert to R vector.
    predict = as.vector(t(kriged[[1]]))

    return(predict)
  }
}

##
# @title gene_krige: Spatial interpolation of gene expression
# @description Performs spatial interpolation ('kriging') of transformed gene
# counts in spatially-resolved transcriptomics data using geoR.
#
# @param geodata
# @param locations
# @param univ
# @return res_krig$predict
#
#
krige_geor = function(geodata=NULL, locations=NULL, univ=NULL){

  if(univ == F){
    # NOTE: Need to check how to decide on init.cov.pars
    geo_lhood <- geoR::likfit(geodata=geodata, trend='cte', ini.cov.pars=c(1, 0.15), messages=F)

    # Specify control (and output) parameters for ordinary kriging.
    KC <- geoR::krige.control(obj.model=geo_lhood)
    OC <- geoR::output.control(messages=F)

    # Perform ordinary kriging.
    res_krig <- geoR::krige.conv(geodata=geodata, locations=locations, krige=KC, output=OC)

    return(res_krig$predict)

  }else if(univ == T){
    # NOTE: Need to use regression analysis of variogram to get values for
    # nugget.
    geo_lhood <- geoR::likfit(geodata=geodata, trend='cte', ini.cov.pars=c(1000, 500), nug=100, messages=F)

    # Specify control (and output) parameters for universal kriging.
    KC <- geoR::krige.control(type.krige="OK", obj.m=geo_lhood, trend.d="cte", trend.l="cte")
    OC <- geoR::output.control(messages=F)

    # Perform universal kriging.
    res_krig <- geoR::krige.conv(geodata=geodata, locations=locations, krige=KC, output=OC)

    return(res_krig$predict)
  }
}

