library(OpenStreetMap)
library(raster)

# Get the OSM data, with caching to speed it up
get_osm_backdrop_reproj<-function(lower_left, upper_right, dbox, backdrop, zoom=12){
    osm_cachefile = paste0('.osm_tile_cache_', upper_right[1], upper_right[2], lower_left[1], lower_left[2], '.RDS')
    dmax = max(dbox)
    if(!file.exists(osm_cachefile)){
        # Get a piece that is bigger than we need, to avoid gaps in plots when
        # the plot limits are extended by R (which depends on the figure size
        # too)
        osm_backdrop = try(openmap(upperLeft = c(upper_right[2] + dmax, lower_left[1]  - dmax), 
                               lowerRight= c(lower_left[2]  - dmax, upper_right[1] + dmax), 
                               type=backdrop, zoom=zoom))
        if(is(osm_backdrop, 'try-error')){
            stop('Download failed')
        }
        # Reproject it
        osm_backdrop_reproj = openproj(osm_backdrop, proj4string(CRS("+init=EPSG:4326")))
        saveRDS(osm_backdrop_reproj, file=osm_cachefile)
    }else{
        # Read from cache
        osm_backdrop_reproj = readRDS(osm_cachefile)
    }
    return(osm_backdrop_reproj)
}
