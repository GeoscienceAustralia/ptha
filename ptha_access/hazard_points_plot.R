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
            legend=FALSE
        )

    #### Quick approach, but can't quite get it to work
    ##for(i in 1:length(unit_source_grids)){
    ##    poly_data = unit_source_grids[[i]]
    ##    poly_name = names(unit_source_grids)[i]
    ##    m = garnishMap(m, addPolygons, data=poly_data, layerId = poly_name, 
    ##        group=poly_name, popup = popupTable(poly_data))
    ##    m_new = mapview(poly_data, map = m@map, layer.name = poly_name)
    ##    m_new2 = methods::new('mapview', 
    ##        object=append(m@object, m_new@object), map=m)
    ##    m = m_new2
    ##}

    for(i in 1:length(unit_source_grids)){
        m = m + mapview(unit_source_grids[[i]], legend=FALSE, 
            color=rainbow(50)[i], alpha.regions=0, 
            layer.name=names(unit_source_grids)[i],
            homebutton=FALSE)
            #homebutton=TRUE, layer.name=names(unit_source_grids)[i])
        #addFeatures(m, unit_source_grids[[i]])
    }

    # Brute force approach -- this leads to the correct source-zone names
    # appearing on the interactive map -- which gives click-zoom functionality
    #
    # Not sure how to achieve this with a loop
    #izumariana = unit_source_grids[[1]]
    #m = m + izumariana
    #kermadectonga = unit_source_grids[[2]]
    #m = m + kermadectonga
    #kurilsjapan = unit_source_grids[[3]]
    #m = m + kurilsjapan
    #newhebrides = unit_source_grids[[4]]
    #m = m + newhebrides
    #puysegur = unit_source_grids[[5]]
    #m = m + puysegur
    #solomon = unit_source_grids[[6]]
    #m = m + solomon
    #southamerica = unit_source_grids[[7]]
    #m = m + southamerica
    #sunda = unit_source_grids[[8]]
    #m = m + sunda

    # Display in browser
    print(m)
    #m@map
    #return(m)

}
plot_interactive_map()
