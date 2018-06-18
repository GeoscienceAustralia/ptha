suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(mapview))
suppressPackageStartupMessages(library(leaflet))
if(!exists('config_env')){
    config_env = new.env()
    source('R/config.R', local=config_env)
}

# If the following is TRUE, then the code will re-download key datasets. Otherwise,
# it will read versions from cache, and download if they do not exist
REFRESH_MAP = config_env$REFRESH_MAP

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
    all_sourcezone_unit_source_grids = Sys.glob('SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')

    unit_source_grids_rds_file = 'SOURCE_ZONES/all_unit_source_grids.RDS'

    if( REFRESH_MAP | (length(all_sourcezone_unit_source_grids) == 0) | 
        (!file.exists(unit_source_grids_rds_file)) ){
        # Download the shapefiles
        get_shapefile_env = new.env()
        source('R/get_supporting_data.R', local=get_shapefile_env)
        get_shapefile_env$download_all_unit_source_grid_shapefiles()
        # Find names of all the shapefiles
        all_sourcezone_unit_source_grids = Sys.glob('SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp')

        # Read each one into a list
        unit_source_grids = list()
        for(i in 1:length(all_sourcezone_unit_source_grids)){
            layer_name = gsub('.shp', '', basename(all_sourcezone_unit_source_grids[i]))
            unit_source_grids[[layer_name]] = readOGR(all_sourcezone_unit_source_grids[i], 
                layer=layer_name, verbose=FALSE)
            #unit_source_grids[[layer_name]] = st_read(all_sourcezone_unit_source_grids[i], 
            #    layer=layer_name, verbose=FALSE)
            names(unit_source_grids[[i]])[1:2] = c('downdip_index', 'alongstrike_index')
            unit_source_grids[[i]]$sourcezone = rep(layer_name, len=length(unit_source_grids[[i]]$downdip_index))
        }

        saveRDS(unit_source_grids, unit_source_grids_rds_file)

    }else{

        unit_source_grids = readRDS(unit_source_grids_rds_file)
    }
   
    return(unit_source_grids)
}
unit_source_grids = .read_all_unit_source_grids()

#' Read the hazard points
#'
#' Wrap in a function to avoid adding many variables to the environment
#'
#' @param refresh If FALSE, then read a cached version of the data if it exists
#' @return hazard_points_spdf SpatialPointsDataFrame containing the hazard points
#'
.read_hazard_points<-function(refresh=REFRESH_MAP){

    hazard_points_spdf_RDS_file = 'DATA/hazard_points_spdf.RDS'
   
    if(refresh | !file.exists(hazard_points_spdf_RDS_file)){ 

        # Read the data, in either netcdf or csv format. 
        # Currently only the latter has return period information
        read_format = 'csv'
        if(read_format == 'netcdf'){
            # To read the data, we need 'get_netcdf_gauge_locations' from the following file 
            source('R/sum_tsunami_unit_sources.R', local=TRUE)

            # Find a file that contains hazard points. Easiest way is to read them from a tide gauge file
            unit_source_stats_alaska = paste0(config_env$.GDATA_OPENDAP_BASE_LOCATION, 
                'SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/unit_source_statistics_alaskaaleutians.nc')
            fid = nc_open(unit_source_stats_alaska)
            tg_filename = ncvar_get(fid, 'tide_gauge_file', start=c(1, 1), count=c(4096,1))[1]
            nc_close(fid)
            # Read the hazard points
            tg_filename = config_env$adjust_path_to_gdata_base_location(tg_filename)
            hazard_points = try(get_netcdf_gauge_locations(tg_filename))

            if(class(hazard_points) == 'try-error'){
                stop('hazard point read failed. This may occur with slower internet connections, so you could try again')
            }

            # Make sure all columns are numeric (some read as 'array' from netcdf)
            for(i in 1:ncol(hazard_points)){
                hazard_points[,i] = as.numeric(hazard_points[,i])
            }

        }else if(read_format == 'csv'){
            #
            # Read a csv version on the thredds server
            #
            tg_filename = paste0(config_env$.GDATA_HTTP_BASE_LOCATION, 
                'EVENT_RATES/tsunami_stages_at_fixed_return_periods.csv')
            tg_filename_local = paste0('DATA/tsunami_stages_at_fixed_return_periods.csv')
            download.file(tg_filename, dest= tg_filename_local)

            hazard_points = try(read.csv(tg_filename))

            if(class(hazard_points) == 'try-error'){
                stop('hazard point read failed. This may occur with slower internet connections, so you could try again')
            }

        }else{
            # 
            stop('Unknown read_format for hazard points')
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

        # Faster read with RDS format
        saveRDS(hazard_points_spdf, hazard_points_spdf_RDS_file)

    }else{
        #
        # The data should be saved locally
        #
        #hazard_points_spdf = readOGR('DATA/hazard_points_spdf', layer='hazard_points_spdf', 
        #    verbose=FALSE)
        hazard_points_spdf = readRDS(hazard_points_spdf_RDS_file)
        # Fix abbreviated name
        k = which(names(hazard_points_spdf) == 'pnt_ctg')
        names(hazard_points_spdf)[k] = 'point_category'

    }

    return(hazard_points_spdf)
}
hazard_points_spdf = .read_hazard_points()


## One method -- leaflet
#library(leaflet)
#m =  leaflet(hazard_points) %>% addTiles() %>% 
#    addMarkers(clusterOptions=markerClusterOptions(
#        maxClusterRadius=20, disableClusteringAtZoom=8, spiderfyOnMaxZoom=FALSE)
#    )

# Make an interactive map where we can view and query points. 
# 'mapview' is built off 'leaflet'
plot_interactive_map<-function(refresh_map = REFRESH_MAP){

    suppressPackageStartupMessages(library(mapview))

    if(refresh_map | !file.exists('DATA/interactive_map.RDS')){
        # Make a new map
        m = mapview(
                hazard_points_spdf, 
                clusterOptions = markerClusterOptions(maxClusterRadius=20, 
                    disableClusteringAtZoom=8, spiderfyOnMaxZoom=FALSE),
                legend=FALSE, zcol='point_category',
                homebutton=FALSE
            )

        for(i in 1:length(unit_source_grids)){
            m = m + mapview(unit_source_grids[[i]], legend=FALSE, 
                color=rainbow(50)[i], alpha.regions=0, 
                layer.name=names(unit_source_grids)[i],
                homebutton=FALSE)
                #homebutton=TRUE, layer.name=names(unit_source_grids)[i])
            #addFeatures(m, unit_source_grids[[i]])
        }

        dir.create('DATA', showWarnings=FALSE)
        saveRDS(m, 'DATA/interactive_map.RDS')

    }else{
        # Read the map from a file
        m = readRDS('DATA/interactive_map.RDS')
    }

    # Display in browser
    print(m)
    # Output if needed
    return(invisible(m))
}
local_map = plot_interactive_map()
