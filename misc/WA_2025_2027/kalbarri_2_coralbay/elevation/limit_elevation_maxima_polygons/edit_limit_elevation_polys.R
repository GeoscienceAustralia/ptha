#
# Add a few polygons to the "limit_elevation_polys" shapefile.
#
# This is useful for polygons that were created using code.
#
library(terra)

limit_elevation_polys = vect('limit_elevation_polys/limit_elevation_polys.shp')

# Location of polygon geometries made using a coastaline contour dataset
extra_polygons = Sys.glob('/g/data/w85/tsunami/DATA/ELEVATION/WA/Australian_Coastline_50k_2024/clip_pieces/cliff_polygon*.csv')
if(!file.exists(extra_polygons[1])){
    # Try GD's home machine
    extra_polygons = Sys.glob('/media/gareth/Data2/DATA/Australian_Coastline_50k_2024/clip_pieces/cliff_polygon*.csv')
}    
extra_polys_limit = -1
extra_polys_id = 1


new_vects = vector(mode='list', length=length(extra_polygons))
for(i in 1:length(extra_polygons)){
    mat = as.matrix(read.csv(extra_polygons[i]))
    new_vects[[i]] = vect(mat, type='polygons', crs=crs(limit_elevation_polys), atts=data.frame(id=extra_polys_id, limit=extra_polys_limit))
}

new_limit_elevation_polys = rbind(
    limit_elevation_polys, do.call(rbind, new_vects))


new_limit_elevation_polys = na.omit(new_limit_elevation_polys, geom=TRUE)

dir.create('limit_elevation_polys_edited', showWarnings=FALSE)
writeVector(new_limit_elevation_polys, 'limit_elevation_polys_edited/limit_elevation_polys_edited.shp', overwrite=TRUE)
