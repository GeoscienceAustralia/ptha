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
plot_interactive_map<-function(){

    #suppressPackageStartupMessages(library(leaflet))
    suppressPackageStartupMessages(library(mapview))
    m = mapview(
            hazard_points_spdf, 
            clusterOptions = markerClusterOptions(maxClusterRadius=20, 
                disableClusteringAtZoom=8, spiderfyOnMaxZoom=FALSE),
            legend=FALSE,
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

    # Display in browser
    print(m)
    #m@map
    #return(m)

}
plot_interactive_map()
