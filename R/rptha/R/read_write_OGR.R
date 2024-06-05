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
#' @param verbose print more
#' @return An object with class of some sp::SpatialXXX object (e.g.
#' SpatialPolygonsDataFrame, SpatialLinesDataFrame, SpatialPointsDataFrame)
#' @export
#' @import sf sp
readOGR<-function(dsn, layer, verbose = TRUE, stringsAsFactors=FALSE){
# Other arguments from rgdal::readOGR were:
#    p4s=NULL, 
#    stringsAsFactors=as.logical(NA), 
#    drop_unsupported_fields=FALSE,
#    pointDropZ=FALSE, dropNULLGeometries=TRUE,
#    useC=TRUE, disambiguateFIDs=FALSE, addCommentsToPolygons=TRUE,
#    encoding=NULL, use_iconv=FALSE, swapAxisOrder=FALSE, require_geomType = NULL,
#    integer64="no.loss", GDAL1_integer64_policy=FALSE, morphFromESRI = NULL,
#    dumpSRS = FALSE, enforce_xy = NULL, D3_if_2D3D_points=FALSE, missing_3D=0)

    sf_geo = read_sf(dsn=dsn, layer=layer, quiet = (!verbose), stringsAsFactors=stringsAsFactors) 

    # The sf object can include 'empty geometries, which we need to remove, or else it fails
    to_keep = which(!st_is_empty(sf_geo))
    sf_geo = sf_geo[to_keep,]
    x = as(sf_geo, 'Spatial')
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

    # Deal with shapefile limitations in field names, following the original writeOGR
    fld_names = names(obj)
    is_shpfile = (driver == 'ESRI Shapefile')
    if (is_shpfile) {
        if (any(nchar(fld_names) > 10)) {
            # Try to shorten file names while maintaining unique names
            fld_names = abbreviate(fld_names, minlength=7)
            warning("Field names abbreviated for ESRI Shapefile driver")
            if (any(nchar(fld_names) > 10)) 
                fld_names = abbreviate(fld_names, minlength=5)
        }
        # fix for dots in DBF field names 121124
        if (length(wh. <- grep("\\.", fld_names) > 0)) {
            fld_names[wh.] = gsub("\\.", "_", fld_names[wh.])
        }
    }
    if (length(fld_names) != length(unique(fld_names)))
       stop("Non-unique field names")

    names(obj) = fld_names

    # Convert from sp::SpatialXXX to sf format    
    y = sf::st_as_sf(obj)
    st_write(y, dsn=dsn, layer=layer, driver=driver, verbose=verbose, delete_layer=overwrite_layer)
}
