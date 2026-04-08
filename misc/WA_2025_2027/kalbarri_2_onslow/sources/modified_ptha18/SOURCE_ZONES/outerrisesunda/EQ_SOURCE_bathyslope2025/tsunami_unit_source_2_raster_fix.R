#' Convert the tsunami unit source to a z-displacement raster
#' 
#' @param tsunami_unit_source output of make_tsunami_unit_source or similar.
#' @param filename Name for output raster file. If NULL, no file is saved.
#' @param saveonly logical. If FALSE, return the raster object as well as
#' saving. If TRUE but filename is not NULL, then save the file, but return NULL.
#' @param tsunami_surface_points_lonlat matrix with 2 columns containing surface points
#' at which tsunami_unit_source$smooth_tsunami_displacement occurs. If NULL, look for
#' this in tsunami_unit_source -- however, to save memory, the latter may be set to NA. In
#' which case this argument must be provided.
#' @param res optional argument 'res' to pass to \code{rasterFromXYZ}
#' @param field_to_output either 'smoothed_vertical_displacement' (to output vertical
#' displacement with smoothing if it was applied in making the tsunami_unit_source), or
#' 'vertical_displacement' to output the vertical component of the Okada displacement vector, or 'easting_displacement'
#' to output the easting component, or 'northing_displacement' to outut the northing component
#' @return Either a raster, or NULL. Can save the raster to a file as a side effect.
#' @export
tsunami_unit_source_2_raster_fix<-function(tsunami_unit_source, filename=NULL, 
    saveonly=FALSE, tsunami_surface_points_lonlat=NULL, res=c(NA,NA), 
    field_to_output='smoothed_vertical_displacement'){

    allowed_fields_to_output = c(
        'smoothed_vertical_displacement',
        'vertical_displacement',
        'easting_displacement',
        'northing_displacement')

    # Check the field_to_output variable is legitimate
    if(! (field_to_output %in% allowed_fields_to_output)){
        msg = paste0('field_to_output was set to an illegal value: ', 
                     field_to_output, '\n It should be one of: ', 
                     paste0(allowed_fields_to_output, collapse=", "))
        stop(msg)
    }

    if(field_to_output != 'smoothed_vertical_displacement'){
        # If the tsunami_source component is NA, then we didn't store the full
        # Okada output, so cannot provide vertical/easting/northing
        if(any(is.na(tsunami_unit_source$tsunami_source)) | (length(tsunami_unit_source$tsunami_source) == 0)){
            msg = paste0(
                'tsunami_unit_source does not contain the vertical/easting/northing displacements \n (',
                'this typically reflects an effort to conserve memory, by passing the \n',
                'argument "minimal_output=TRUE" to the function make_tsunami_unit_source). \n',
                'To save these displacements you must use "minimal_output=FALSE" \n',
                'when creating the tsunami_unit_source. See ?make_tsunami_unit_source')
            stop(msg)
        }
    }

    # Get the xy output coordinates, depending on whether they were provided to
    # the function
    if(all(is.null(tsunami_surface_points_lonlat))){

        stopifnot(
            length(tsunami_unit_source$tsunami_surface_points_lonlat[,1]) == 
            length(tsunami_unit_source$smooth_tsunami_displacement))
        xy = tsunami_unit_source$tsunami_surface_points_lonlat
    }else{

        stopifnot(length(tsunami_surface_points_lonlat[,1]) == 
                  length(tsunami_unit_source$smooth_tsunami_displacement))
        xy = tsunami_surface_points_lonlat
    }

    # Make 'xyz' matrix with the data to be rasterized
    xyz = switch(field_to_output, 
        "smoothed_vertical_displacement" = 
            cbind(xy, tsunami_unit_source$smooth_tsunami_displacement),
        "vertical_displacement" = 
            cbind(xy, tsunami_unit_source$tsunami_source$zdsp),
        "easting_displacement" = 
            cbind(xy, tsunami_unit_source$tsunami_source$edsp),
        "northing_displacement" = 
            cbind(xy, tsunami_unit_source$tsunami_source$ndsp)
        )

    # Currently use local version of rasterFromXYZ to deal with some bugs in it
    # (which have been reported)
    outrast = rptha:::.local_rasterFromXYZ(xyz, crs=CRS('+init=epsg:4326'), res=res)
    #outrast = rasterFromXYZ(xyz, crs=CRS('+init=epsg:4326'), res=res)

    # Free up memory
    rm(xy, xyz); gc()

    if(!is.null(filename)){
        writeRaster(outrast, filename, driver='GTiff', 
            options='COMPRESS=DEFLATE', overwrite=TRUE)
    }

    if(saveonly){
        return()
    }else{
        return(outrast)
    }
}


