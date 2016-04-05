library(raster)
library(rgdal)
 
#' Contouring with gdal_contour
#'   
#' Using gdal to make a contour can be more reliable and efficient than R's
#' algorithms for this application
#'
#' @param dem raster for which contours are desired, either as a filename or a
#' RasterLayer. In the former case any format supported by gdal is ok.
#' @param contour_levels elevations for the contours
#' @param out_dsn Folder to write outputs to
#' @param contour_shp_name name for the contour shapefile (inside out_dsn)
#' @param delete_tmp_files TRUE/FALSE -- delete temp files??
#' @param verbose TRUE/FALSE -- print more or less
#' @param clean_dsn TRUE/FALSE -- delete files matching contour_shp_name
#' in out_dsn. If we dont do this, gdal will throw errors which may not be
#' detected
#' @param gdal_contour_call character. R will try to call gdal_contour by
#' passing this command (followed by various flags) to 'system'
#' @return A SpatialLinesDataFrame contour, which seems to be of higher quality
#' (and execute faster) than the native R routines
#'
gdal_contour<-function(
    dem, 
    contour_levels=0.0, 
    out_dsn=NA, 
    contour_shp_name=NA,
    delete_tmp_files=is.na(out_dsn),
    verbose=TRUE, 
    clean_dsn=TRUE,
    gdal_contour_call = 'gdal_contour '){

    # Make directory to write file in, if there is something
    if(is.na(out_dsn)){
        out_dsn = '.GCONTOUR_TMP'
    } 
    dir.create(out_dsn, showWarnings=FALSE)
    
    # Make sure the dem is a tif file 
    if(class(dem) == 'RasterLayer'){ 
        # Write the dem to a file
        if(verbose) print('Writing dem to file ...')
        tmp_raster_filename = paste0(out_dsn,'/temp_dem.tif')
        writeRaster(dem, tmp_raster_filename, driver='GTiff', overwrite=TRUE)
    }else{
        if(file.exists(dem)){
            tmp_raster_filename = dem
        }else{
            stop('ERROR: dem should either be a raster layer, or a .tif filename')
        }
    }
    
    # Make the contour layer name
    if(!is.na(contour_shp_name)){
        contour_layer = contour_shp_name
    }else{
        contour_layer = 'contour'
    }

    # Delete all files in out_dsn called contour_shape_name -- or there will be
    # problems
    if(clean_dsn){
        file.remove(dir(out_dsn, pattern=contour_layer, full.names=TRUE))
    }

    # Run the command
    gdal_contour_command = paste(gdal_contour_call, ' -a elev -fl', 
        paste(contour_levels, collapse=" "), '-nln ', contour_layer, 
        tmp_raster_filename, out_dsn)

    if(verbose){
        print('Running gdal_contour ...')
        print(gdal_contour_command)
    } 
    system(gdal_contour_command)

    # Read the file
    out = readOGR(out_dsn, layer=contour_layer, verbose=verbose)
  
    # Clean up 
    if(delete_tmp_files){
        unlink(out_dsn, recursive=TRUE)
    } 

    return(out)
}



#' Convert a SpatialLinesDataFrame containing closed contours to a
#' SpatialPolygonsDataFrame.
#' 
#' This is useful so we can use 'over' to identify lines which don't contain 
#' coastal points. 
#' 
#' @param sldf A SpatialLinesDataFrame
#' @return A SpatialPolygonsDataFrame
#'
SpatialLinesDF2Polygons<-function(sldf){

    mypolylist = lapply(sldf@lines, 
       myfun<-function(x){
            mycrds = x@Lines[[1]]@coords
            l = length(mycrds[,1])
            if(any(mycrds[1,1:2] != mycrds[l,1:2])){
                mycrds = rbind(mycrds, mycrds[1,])
            }
            return(Polygons(list(Polygon(mycrds, hole=FALSE)), ID=1))
        })

    for(i in 1:length(mypolylist)){
       mypolylist[[i]]@ID = as.character(i)
    }

    sldf_poly = SpatialPolygons(mypolylist, proj4string=CRS(proj4string(sldf)))
    sldf_poly = SpatialPolygonsDataFrame(sldf_poly, 
        data=data.frame(ID=1:length(sldf_poly)))

    return(sldf_poly)
}
