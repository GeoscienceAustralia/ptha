#
# We want to find the area of the largest land-mass for each country, and
# flag it in a shapefile
#
# Why? So that we don't cut out any countries when simplifying polygons

library(rgdal)
library(geosphere)

wb = readOGR('TM_WORLD_BORDERS-0.3', layer='TM_WORLD_BORDERS-0.3')

# NOTE: geosphere::areaPolygon gives different areas to the AREA attribute of
# wb. However, from my checks (e.g. googling country areas on the web)
# areaPolygon seems more accurate.  BUT, for Antarctica, there is a problem
# with areaPolygon. (Near poles / self intersection / ...??).  Regardless,
# elsewhere it is sensible.

ncountries = length(wb@data[,1])
land_area_list = list()
land_area_maxind = rep(NA, ncountries)
land_area_coord = matrix(NA, ncol=2, nrow=ncountries)

# Loop over every country, and find the area of each polygon in m^2
for(i in 1:ncountries){
    print(i)
    land_area_list[[i]] = unlist(lapply(wb[i,]@polygons[[1]]@Polygons, 
        myfun<-function(x) areaPolygon(x@coords)) )
    land_area_maxind[i] = which.max(land_area_list[[i]])
    land_area_coord[i,] = 
        wb[i,]@polygons[[1]]@Polygons[[land_area_maxind[i]]]@coords[1,]
}

land_area_sp = SpatialPointsDataFrame(land_area_coord, data=wb@data, 
    proj4string=CRS(proj4string(wb)))

writeOGR(land_area_sp,dsn='OUTPUTS/LARGEST_ISLANDS', 
         layer='LARGEST_ISLANDS',driver='ESRI Shapefile')
