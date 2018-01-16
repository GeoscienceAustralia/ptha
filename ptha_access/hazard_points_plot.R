suppressPackageStartupMessages(library(rgdal))
if(!exists('.HAVE_SOURCED_CONFIG')) source('R/config.R', local=TRUE, chdir=FALSE)

## One method -- leaflet
#library(leaflet)
#m =  leaflet(hazard_points) %>% addTiles() %>% 
#    addMarkers(clusterOptions=markerClusterOptions(
#        maxClusterRadius=20, disableClusteringAtZoom=8, spiderfyOnMaxZoom=FALSE)
#    )


# Make an interactive map where we can view and query points. 
# 'mapview' is built off 'leaflet'
plot_interactive_map<-function(refresh_map = FALSE){

    suppressPackageStartupMessages(library(mapview))

    if(refresh_map | !file.exists('DATA/interactive_map.RDS')){
        # Mape a new map
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
