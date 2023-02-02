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
#' @import sf sp
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
    x = as(read_sf(dsn=dsn, layer=layer, quiet = (!verbose)), 'Spatial')
    return(x)
}

#' Write an sp object to a spatial format
#'
#' This function was previously provided by rgdal::writeOGR. But rgdal is retiring in 2023
#' so we use as a workaround for the rptha package.
#'
#' @param obj An object of class sf::SpatialXXXXX (e.g. SpatialPolygonsDataFrame)
#' @param dsn Output folder for the data
#' @param layer Layer name for the data
#' @param driver for the output data
#' @param verbose print more info
#' @param overwrite_layer If FALSE and the file exists, then this will throw an error. Will overwrite if TRUE.
#' @return Invisibly returns the object in sf format (from st_as_sf(obj) )
#' @export
#' @import sf
writeOGR<-function(obj, dsn, layer, driver, verbose=FALSE, overwrite_layer=FALSE){
# Arguments from rgdal::writeOGR were
#      (obj, dsn, layer, driver, dataset_options = NULL,
#      layer_options=NULL, verbose = FALSE, check_exists=NULL,
#      overwrite_layer=FALSE, delete_dsn=FALSE, morphToESRI=NULL,
#      encoding=NULL, shp_edge_case_fix=FALSE, dumpSRS = FALSE)

    # Convert from sp::SpatialXXX to sf format    
    y = sf::st_as_sf(obj)
    st_write(y, dsn=dsn, layer=layer, driver=driver, verbose=verbose, delete_layer=overwrite_layer)
}
