suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(sf))

#
# Define the initial part of key web-addresses where our data can be accessed
#

# Netcdf data -- this location allows remote query/subsetting
#
# We include a flag that ensures strings are long enough
# FIXME: Can reduce stringlength
max_stringlength = '[stringlength=4096]'
#
.GDATA_OPENDAP_BASE_LOCATION = paste0(max_stringlength, 
    'http://dapds00.nci.org.au/thredds/dodsC/fj6/PTHA/AustPTHA_1/')
#
# Non-netcdf data -- this location only allows download
.GDATA_HTTP_BASE_LOCATION = 'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/'

#'
#' Some files on gdata contain paths like /g/data/fj6/PTHA/AustPTHA_1/....
#' We want to replace the start of this path with the .GDATA_OPENDAP_BASE_LOCATION
#'
#' @param filepath The path of the file, starting with /g/data/fj6/....
#' @return The same filepath, with prefix suitable for remote access
#'
adjust_path_to_gdata_base_location<-function(filepath, prefix_type='opendap'){

    if(length(filepath) != 1) stop('adjust_path_to_gdata_base_location is not vectorized')

    split_path = strsplit(filepath,'/')[[1]]
    k = grep('AustPTHA_1', split_path)
    ls = length(split_path)
    if(length(k) == 0) stop('did not find AustPTHA_1 in path')

    if(k < ls){
        new_base_path = paste(split_path[(k+1):ls], collapse='/')
    }else{
        new_base_path = ''
    }

    if(prefix_type == 'opendap'){
        new_base_path = paste0(.GDATA_OPENDAP_BASE_LOCATION, new_base_path)
    }
    if(prefix_type == 'http'){
        new_base_path = paste0(.GDATA_HTTP_BASE_LOCATION, new_base_path)
    }

    return(new_base_path)

}


#
# Get base datasets
#

#' Read the unit-source-grid polygons
#'
#' Wrap in a function to avoid adding many variables to the environment
#'
#' @return unit_source_grids list containing a set of
#' SpatialPolygonsDataFrames, one for each source-zone
#'
.read_all_unit_source_grids<-function(){
    # Find names of all the shapefiles
    .all_sourcezone_unit_source_grids = Sys.glob('SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')

    # Read each one into a list
    unit_source_grids = list()
    for(i in 1:length(.all_sourcezone_unit_source_grids)){
        layer_name = gsub('.shp', '', basename(.all_sourcezone_unit_source_grids[i]))
        unit_source_grids[[layer_name]] = readOGR(.all_sourcezone_unit_source_grids[i], 
            layer=layer_name, verbose=FALSE)
        #unit_source_grids[[layer_name]] = st_read(.all_sourcezone_unit_source_grids[i], 
        #    layer=layer_name, verbose=FALSE)
        names(unit_source_grids[[i]])[1:2] = c('downdip_index', 'alongstrike_index')
        unit_source_grids[[i]]$sourcezone = rep(layer_name, len=length(unit_source_grids[[i]]$downdip_index))
    }

    return(unit_source_grids)
}
unit_source_grids = .read_all_unit_source_grids()

#' Read the hazard points
#'
#' Wrap in a function to avoid adding many variables to the environment
#'
#' @return hazard_points_spdf SpatialPointsDataFrame containing the hazard points
.read_hazard_points<-function(refresh=FALSE){

    source('R/sum_tsunami_unit_sources.R', local=TRUE)
   
    if(refresh | !file.exists('DATA/hazard_points_spdf/hazard_points_spdf.shp')){ 

        # Find a file that contains hazard points. Easiest way is to read them from a tide gauge file
        unit_source_stats_alaska = paste0(.GDATA_OPENDAP_BASE_LOCATION, 
            'SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/unit_source_statistics_alaskaaleutians.nc')
        fid = nc_open(unit_source_stats_alaska)
        tg_filename = ncvar_get(fid, 'tide_gauge_file', start=c(1, 1), count=c(4096,1))[1]
        nc_close(fid)
        # Read the hazard points
        tg_filename = adjust_path_to_gdata_base_location(tg_filename)
        hazard_points = try(get_netcdf_gauge_locations(tg_filename))

        # Make sure all columns are numeric (some read as 'array' from netcdf)
        for(i in 1:ncol(hazard_points)){
            hazard_points[,i] = as.numeric(hazard_points[,i])
        }

        if(class(hazard_points) == 'try-error'){
            stop('hazard point read failed. This may occur with slower internet connections, so you could try again')
        }

        # The data contains an numeric 'gaugeID'. It is a decimal number. The fractional
        # part 
        hp_type = match( 
            round(hazard_points$gaugeID - trunc(hazard_points$gaugeID), 1)*10,
            c(0, 1, 2, 3, 4, 5))
        hp_type_char = c('shallow', 'intermediate', 'deep', 'intermediateG', 'DART', 'gridded')[hp_type]
        hazard_points = cbind(hazard_points, data.frame(point_category=hp_type_char))

        hazard_points_spdf = SpatialPointsDataFrame(coords = hazard_points[,1:2], 
            data=hazard_points, proj4string=CRS('+init=epsg:4326'))

        # Only display a subset of points -- we could do this in a more complex way easily
        #kk = which(hp_type_char != 'intermediateG')
        #hazard_points_spdf = hazard_points_spdf[kk,]
        #browser()
        clip_points = FALSE
        if(clip_points){
            stop('need to provide point_filter_polygon')
            dartp = which(hp_type_char == 'DART')

            clip_region = readOGR(dsn='DATA/HAZARD_POINTS/point_filter_polygon', 
                layer='point_filter_polygon', verbose=FALSE)
            suppressWarnings({proj4string(clip_region) = proj4string(hazard_points_spdf)})

            clip_region_keep = which(!is.na(over(as(hazard_points_spdf, 'SpatialPoints'), clip_region)))

            hazard_points_spdf = hazard_points_spdf[c(clip_region_keep, dartp),]
        }


        writeOGR(hazard_points_spdf, dsn='DATA/hazard_points_spdf', 
            layer='hazard_points_spdf', driver='ESRI Shapefile',
            overwrite=TRUE)

    }else{

        hazard_points_spdf = readOGR('DATA/hazard_points_spdf', layer='hazard_points_spdf', 
            verbose=FALSE)
        # Fix abbreviated name
        k = which(names(hazard_points_spdf) == 'pnt_ctg')
        names(hazard_points_spdf)[k] = 'point_category'

    }

    return(hazard_points_spdf)
}
hazard_points_spdf = .read_hazard_points()


# Use this variable to determine if we have already sourced config.R
.HAVE_SOURCED_CONFIG=TRUE
