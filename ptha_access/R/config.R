suppressPackageStartupMessages(library(rgdal))

#
# Define the initial part of key web-addresses where our data can be accessed
#

# Netcdf data -- this location allows remote query/subsetting
.GDATA_OPENDAP_BASE_LOCATION = 'http://dapds00.nci.org.au/thredds/dodsC/fj6/PTHA/AustPTHA/v2017/'

# Non-netcdf data -- this location only allows download
.GDATA_HTTP_BASE_LOCATION = 'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA/v2017/'


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
        names(unit_source_grids[[i]]) = c('downdip_index', 'alongstrike_index')
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
.read_hazard_points<-function(){
    # Read as csv
    hazard_points = read.csv('DATA/HAZARD_POINTS/merged_hazard_points.csv')
    # The 3rd column contains an numeric 'ID'. It is a decimal number. The fractional
    # part 
    hp_type = round(hazard_points$ID - trunc(hazard_points$ID), 1)*10
    hp_type_char = c('shallow', 'intermediate', 'deep', 'intermediateG', 'DART')[hp_type+1]
    hazard_points = cbind(hazard_points, data.frame(point_category=hp_type_char))

    hazard_points_spdf = SpatialPointsDataFrame(coords = hazard_points[,1:2], 
        data=hazard_points[,-c(1:2)], proj4string=CRS('+init=epsg:4326'))

    # Only display a subset of points -- we could do this in a more complex way easily
    #kk = which(hp_type_char != 'intermediateG')
    #hazard_points_spdf = hazard_points_spdf[kk,]
    #browser()
    dartp = which(hp_type_char == 'DART')

    clip_region = readOGR(dsn='DATA/HAZARD_POINTS/point_filter_polygon', layer='point_filter_polygon')
    suppressWarnings({proj4string(clip_region) = proj4string(hazard_points_spdf)})

    clip_region_keep = which(!is.na(over(as(hazard_points_spdf, 'SpatialPoints'), clip_region)))

    hazard_points_spdf = hazard_points_spdf[c(clip_region_keep, dartp),]

    return(hazard_points_spdf)
}
hazard_points_spdf = .read_hazard_points()


# Use this variable to determine if we have already sourced config.R
.HAVE_SOURCED_CONFIG=TRUE
