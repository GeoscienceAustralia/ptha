#
# Get arrival times at a set of "thinned" offshore and onshore points.
# The intention is to make the arrival time points easier to work with, vs using
# the rasters.
#
library(terra)

# List with names of source zones (here we do not actually use the entry values)
plot_metadata = list(
    alaskaaleutians = list(xlim = c(10, 16)),
    kermadectonga2 = list(xlim = c(0,7)),
    newhebrides2  = list(xlim=c(0,7)),
    outerrise_kermadectonga = list(xlim=c(0,7)),
    outerrisenewhebrides = list(xlim=c(0,7)),
    outerrise_puysegur = list(xlim=c(0,4)),
    puysegur2 = list(xlim=c(0,4)),
    solomon2 = list(xlim=c(0,6)),
    southamerica = list(xlim=c(11, 19))
)

store_offshore_elev = list()
store_onshore_elev = list()
    
for(nm in names(plot_metadata)){    
    print(nm)
    for(arrival_type in c('minimum', 'mean')){
        print(arrival_type)

        # Get arrival time tiles ONLY for high-res domains
        MIN_DOMAIN_INDEX = 261
        kt_arrival_times = Sys.glob(paste0('../NSW_tsunami_modelling_project_final_outputs/Arrival_times/', 
            nm, '/', paste0(arrival_type), '_arrival_time/*.tif'))
        kt_index = unlist(lapply(strsplit(kt_arrival_times, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))
        keep = which(kt_index >= MIN_DOMAIN_INDEX)
        elevation_files = Sys.glob('../NSW_tsunami_modelling_project_final_outputs/elevation_in_model/*.tif')
        elevation_index = unlist(lapply(strsplit(elevation_files, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))

        # Constraints for offshore contour
        OFFSHORE_CONTOUR_LEVEL = -20
        OFFSHORE_INCLUSION_BELOW_THIS_ELEVATION = 0

        # Constraints for onshore contour
        # FIXME: This doesn't achieve what we want {to have the 'shortest
        # onshore arrival times'}, probably a different approach is needed
        ONSHORE_CONTOUR_LEVEL = 1.3
        ONSHORE_INCLUSION_ABOVE_THIS_ELEVATION = 1.15
        ONSHORE_EXCLUSION_ABOVE_THIS_ELEVATION = 1.4

        # Ensure order matches arrival times
        elv_match = match(elevation_index, kt_index)
        elevation_files = elevation_files[elv_match]
        elevation_index = elevation_index[elv_match]
        kt_arrival_times = kt_arrival_times[keep]
        kt_index = kt_index[keep]
        elevation_files = elevation_files[keep]
        elevation_index = elevation_index[keep]

        get_arrival_times_on_contour_level<-function(arrival_time_tif, elevation_tif, contour_level){
            elev_tile = rast(elevation_tif)
            # Make contour -- maybe better to use gdal_contour to control smoothing?
            contour_tile_1 = try(as.contour(elev_tile, levels=contour_level, maxcells=1e+07))
            if(is(contour_tile_1, 'try-error')){
                # Fail gracefully when we don't have elevation values near contour_level
                return(data.frame(arrival_time=c(), x=c(), y=c()))
            }
            contour_pts = crds(as.points(contour_tile_1)) # matrix with columns lon,lat
            arrival_time = extract(rast(arrival_time_tif), contour_pts)
            elev = extract(elev_tile, contour_pts)
            output = data.frame(arrival_time = arrival_time[[1]], elev = elev[[1]], x=contour_pts[,1], y=contour_pts[,2])
            return(output)
        }

        library(parallel)
        parfun<-function(i){
            # Get arrival times around contours just offshore and just onshore
            offshore = get_arrival_times_on_contour_level(kt_arrival_times[i], elevation_files[i], contour_level=OFFSHORE_CONTOUR_LEVEL)
            onshore = get_arrival_times_on_contour_level(kt_arrival_times[i], elevation_files[i], contour_level=ONSHORE_CONTOUR_LEVEL)
            return(list(offshore=offshore, onshore=onshore))
        }

        all_at = mclapply(1:length(kt_arrival_times), parfun, mc.cores=16)

        offshore_at = do.call(rbind, lapply(all_at, function(x) x$offshore))
        k = which(offshore_at$elev < OFFSHORE_INCLUSION_BELOW_THIS_ELEVATION)
        offshore_at = offshore_at[k,]

        # For onshore points, prevent them from being too high up
        onshore_at = do.call(rbind, lapply(all_at, function(x) x$onshore))
        k = which(onshore_at$elev > ONSHORE_INCLUSION_ABOVE_THIS_ELEVATION)
        onshore_at = onshore_at[k,]
        k = which(onshore_at$elev < ONSHORE_EXCLUSION_ABOVE_THIS_ELEVATION)
        onshore_at = onshore_at[k,]

        store_offshore_elev[[paste0(nm, '_', arrival_type)]] = offshore_at
        store_onshore_elev[[paste0(nm, '_', arrival_type)]]  = onshore_at
    }
}

saveRDS(store_offshore_elev, 'offshore_arrival_time_storage.RDS')

# The onshore arrival times should have further processing.
#saveRDS(store_onshore_elev,  'onshore_arrival_time_storage.RDS')

thin_points_list<-function(points_list){

    # First check that all list entries have the same x/y points
    N = length(points_list)
    for(i in 2:N){
        same_x = all(points_list[[i]]$x == points_list[[1]]$x)
        same_y = all(points_list[[i]]$y == points_list[[1]]$y)
        if(!(same_x & same_y)) stop('Points are not identical')
    }
    
    # Make a single set of points
    lon = points_list[[1]]$x
    lat = points_list[[1]]$y
    elev = points_list[[1]]$elev
    output_elev_columns = lapply(points_list, function(x) x$arrival_time)
    names(output_elev_columns) = names(points_list)


    final_list = list(lon = lon, lat=lat, elev=elev)
    final_list = c(final_list, output_elev_columns)
    final_list = as.data.frame(final_list)

    library(GeoThinneR)
    thinned_list = thin_points(final_list, lon_col = 'lon', lat_col = 'lat', method='distance', trials=1, thin_dist = 1)

    thinned_points = final_list[thinned_list$retained[[1]],]
    
    return(thinned_points)
}

thinned_offshore_elev = thin_points_list(store_offshore_elev)
write.csv(thinned_offshore_elev, 'offshore_arrival_times_all_sources_thinned.csv', row.names=FALSE)


#thinned_onshore_elev = thin_points_list(store_onshore_elev)
#write.csv(thinned_onshore_elev, 'onshore_arrival_times_all_sources_thinned.csv', row.names=FALSE)
