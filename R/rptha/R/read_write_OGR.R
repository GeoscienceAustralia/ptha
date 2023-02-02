#
# In 2023 the R package 'rgdal' is retiring.
# This is widely used in "rptha" for 2 purposes:
#    readOGR, writeOGR
# Here we provide alternative interfaces to those
#

#' Read OGR supported file as an sp::SpatialXXXXX object
#'
#' This function was previously provided by rgdal::readOGR. But rgdal is retiring in 2023
#' so we use as a workaround for the rptha package.
#'
#' @param dsn Data source destination (often a folder, or for shapefiles can be "name_of_shapefile.shp")
#' @param layer data layer (e.g. for a shapefile "myfile" inside dsn, it is probably "myfile")
#' @param verbose 
#' @return An object with class of some sp::SpatialXXX object (e.g.
#' SpatialPolygonsDataFrame, SpatialLinesDataFrame, SpatialPointsDataFrame)
#' @export
readOGR<-function(dsn, layer, verbose = TRUE){
# Other arguments from rgdal::readOGR were:
#    p4s=NULL, 
#    stringsAsFactors=as.logical(NA), 
#    drop_unsupported_fields=FALSE,
#    pointDropZ=FALSE, dropNULLGeometries=TRUE,
#    useC=TRUE, disambiguateFIDs=FALSE, addCommentsToPolygons=TRUE,
#    encoding=NULL, use_iconv=FALSE, swapAxisOrder=FALSE, require_geomType = NULL,
#    integer64="no.loss", GDAL1_integer64_policy=FALSE, morphFromESRI = NULL,
#    dumpSRS = FALSE, enforce_xy = NULL, D3_if_2D3D_points=FALSE, missing_3D=0)
    require(sf)
    require(sp) 
    x = as(read_sf(dsn, layer, quiet = (!verbose)), 'Spatial')
    return(x)
}

#' Write an sp object to a spatial format
#'
#' This function was previously provided by rgdal::writeOGR. But rgdal is retiring in 2023
#' so we use as a workaround for the rptha package.
#'
writeOGR<-function(obj, dsn, layer, driver, verbose=FALSE, overwrite_layer=FALSE){
# Arguments from rgdal::writeOGR were
#      (obj, dsn, layer, driver, dataset_options = NULL,
#      layer_options=NULL, verbose = FALSE, check_exists=NULL,
#      overwrite_layer=FALSE, delete_dsn=FALSE, morphToESRI=NULL,
#      encoding=NULL, shp_edge_case_fix=FALSE, dumpSRS = FALSE)

    require(sf)
    require(sp)
    # Convert from sp::SpatialXXX to sf format    
    y = st_as_sf(obj)
    write_sf(y, dsn, layer, driver, verbose=verbose, append=(!overwrite_layer))  
}
