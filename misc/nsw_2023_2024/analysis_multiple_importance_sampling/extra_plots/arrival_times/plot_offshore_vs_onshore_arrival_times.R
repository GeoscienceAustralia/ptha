# Show how the arrival times vary from offshore to onshore
library(terra)

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
    
for(nm in names(plot_metadata)){    
    for(arrival_type in c('minimum', 'mean')){

        # Get arrival time tiles ONLY for high-res domains
        MIN_DOMAIN_INDEX = 261
        kt_arrival_times = Sys.glob(paste0('../../NSW_tsunami_modelling_project_final_outputs_20241106/Arrival_times/', 
            nm, '/', paste0(arrival_type), '_arrival_time/*.tif'))
        kt_index = unlist(lapply(strsplit(kt_arrival_times, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))
        keep = which(kt_index >= MIN_DOMAIN_INDEX)
        elevation_files = Sys.glob('../../NSW_tsunami_modelling_project_final_outputs_20241106/elevation_in_model/*.tif')
        elevation_index = unlist(lapply(strsplit(elevation_files, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))

        OFFSHORE_CONTOUR_LEVEL = -2
        ONSHORE_CONTOUR_LEVEL = 1.3
        OFFSHORE_INCLUSION_BELOW_THIS_ELEVATION = 0
        ONSHORE_INCLUSION_ABOVE_THIS_ELEVATION = 1.3

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
        onshore_at = do.call(rbind, lapply(all_at, function(x) x$onshore))

        if(arrival_type == 'minimum') labstart = 'Minimum'
        if(arrival_type == 'mean') labstart = 'Mean'

        png(paste0('offshore_vs_onshore_', nm, '_', arrival_type, '.png'), width=9, height=8, units='in', res=300)
        # To highlight the difference between offshore and onshore, we ensure the 'offshore' points have elevation < 0,
        # while the 'onshore' points have elevation > 0.8. 
        k = which(offshore_at$elev < OFFSHORE_INCLUSION_BELOW_THIS_ELEVATION)
        plot(offshore_at$arrival_time[k]/3600, offshore_at$y[k], pch=19, cex=0.2, 
            xlim=plot_metadata[[nm]]$xlim, ylim=c(-38, -26), 
            xlab=paste0(labstart, ' arrival time (hours post earthquake)'), 
            ylab='Latitude', cex.lab=1.5, cex.axis=1.5)
        k = which(onshore_at$elev > ONSHORE_INCLUSION_ABOVE_THIS_ELEVATION)
        points(onshore_at$arrival_time[k]/3600, onshore_at$y[k], pch=19, cex=0.2, col='red')
        grid(col='orange')
        title(main=paste0(nm, ' : ', labstart, ' arrival times near the shoreline, ', 
            ' \n (includes mainland NSW with estuaries + Lord Howe & Norfolk)'), cex.main=1.6)

        source('cities.R')
        text(rep(plot_metadata[[nm]]$xlim[1] + 0.05, length(cities$lon)), cities$lat, cities$name, cex=1.4, pos=4)
        legend('topright', c(
            paste0('Offshore, elevation < ', OFFSHORE_INCLUSION_BELOW_THIS_ELEVATION, ' m'), 
            paste0('Onshore, elevation > ', ONSHORE_INCLUSION_ABOVE_THIS_ELEVATION, ' m')), 
            col=c('black', 'red'), pch=19, cex=1.4, bg='white')
        dev.off()

        store_offshore_elev[[paste0(nm, '_', arrival_type)]] = offshore_at
    }
}

saveRDS(store_offshore_elev, 'arrival_time_storage.RDS')
