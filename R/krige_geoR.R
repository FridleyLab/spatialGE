##
# @title gene_krige: Spatial interpolation of gene expression
# @description Performs spatial interpolation ('kriging') of normalized gene
# counts in spatially-resolved transcriptomics data.
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
