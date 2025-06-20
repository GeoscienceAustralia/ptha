# Find gauges that occur within a high-resolution model domain.
#
# In practice the analysis is restricted to gauges that are inside the high-resolution parts of the model.
# There are 5 kinds of model setups:
# - 'highres_australia' -- All the high-res areas (NSW, Vic, Perth)
# - 'highres_NSW' -- Only NSW high-res areas
# - 'highres_perth' -- Only greater perth high-res area
# - 'highres_SWWA' -- More extensive WA highres area
# - 'highres_australiaSWWA' -- All the high-res areas, with SWWA instead of Perth (so more extensive in WA).
#
# The first stage of our data export scripts did not restrict the gauges to the high-res domains.
#   + When the modelling began I had less data, and it wasn't necessary.
#   + Over the years I've got more data, and now it's necessary, but can be implemented here.
# 
library(sf)

make_highres_domains<-function(){
    all_poly = list(
        # Here the names correspond to possible values of 'model_resolution_tag'
        #highres_australia = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_australia/domains_shapefile/domains_shapefile.shp'),
        highres_NSW = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_NSW/domains_shapefile/domains_shapefile.shp'),
        highres_perth = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_perth/domains_shapefile/domains_shapefile.shp'),
        #highres_SWWA = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_SWWA/domains_shapefile/domains_shapefile.shp'),
        #highres_australiaSWWA = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_australiaSWWA/domains_shapefile/domains_shapefile.shp'),
        highres_australiaWA = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_australiaWA/domains_shapefile/domains_shapefile.shp'),
        highres_NWWA = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_NWWA/domains_shapefile/domains_shapefile.shp'),
        highres_WA = read_sf('../../swals/domain_shapefiles_by_highres_regions/highres_WA/domains_shapefile/domains_shapefile.shp')
    )
    
    all_highres_poly = lapply(all_poly, function(x) x[x$dx < (1/(60*7*4)),]) # In practice this will only select highres domains
    return(all_highres_poly)
}
HIGHRES_DOMAINS = make_highres_domains()

is_gauge_in_highres_domain<-function(gauge_lon, gauge_lat, model_resolution_tag){
    stopifnot( (length(gauge_lon) == 1) & (length(gauge_lat) == 1) )
    stopifnot(model_resolution_tag %in% names(HIGHRES_DOMAINS))
    check_intersection = st_intersects(st_point(c(gauge_lon, gauge_lat)), HIGHRES_DOMAINS[[model_resolution_tag]])
    result = ( length(check_intersection[[1]]) > 0 )
    return(result)
}
