#
# Make a shapefile defining raster extents for a "merged model".
# This should have small/highres boxes whereever any priority model has them.
#

library(sf)

# Domain shapefiles
all_shp = list(
    bunbury_busselton_alternative = "../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/bunbury_busselton_alternative_with_bunbury_floodgate_shut/domains_shapefile/domains_shapefile/domains_shapefile.shp",
    bunbury_busselton = "../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/bunbury_busselton/domains_shapefile/domains_shapefile/domains_shapefile.shp",
    greater_perth = "../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/domains_shapefile/domains_shapefile/domains_shapefile.shp",
    midwest = "../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/midwest/domains_shapefile/domains_shapefile/domains_shapefile.shp")

a_s = lapply(all_shp, st_read)

# Find the main 3 models (the bunbury....alternative model is only used in a limited way)
merge_sites = which(names(a_s) %in% c('bunbury_busselton', 'greater_perth', 'midwest'))

merge_and_clean=function(shp1, shp2){
    # Intersect the polygons
    merged_shp = st_intersection(shp1, shp2)
    # Ensure geometries are valid
    merged_shp = st_make_valid(merged_shp)
    # Get rid of "non-polygon" intersection geometries which tend to arise
    geo_class = unlist(lapply(1:nrow(merged_shp), function(i) class(st_geometry(merged_shp[i,]))[1]))
    keep = which(geo_class == 'sfc_POLYGON')
    merged_shp = merged_shp[keep,]
    return(merged_shp)
}

merge_and_clean_b=function(shp1, shp2){

    # The variable of interest was not computed on global domains
    k = which(shp1$dx > 0.9/60)
    if(length(k) > 0) shp1 = shp1[-k,]
    k = which(shp2$dx > 0.9/60)
    if(length(k) > 0) shp2 = shp2[-k,]

    # Get rid of overlapping domains
    shp1_geometry_hash = as.character(shp1$geometry)
    shp2_geometry_hash = as.character(shp2$geometry)
    to_remove = which(shp2_geometry_hash %in% shp1_geometry_hash)
    shp2 = shp2[-to_remove,]

    return(rbind(shp1, shp2))
}

merged_zones = merge_and_clean_b(a_s$bunbury_busselton, a_s$greater_perth)
merged_zones = merge_and_clean_b(merged_zones, a_s$midwest)

merged_zones$outfile = paste0('merged_raster_', 1:nrow(merged_zones), '.tif')

st_write(merged_zones, dsn='merged_model_raster_output_zones', layer='merged_model_raster_output_zones', driver='ESRI Shapefile', append=FALSE)
