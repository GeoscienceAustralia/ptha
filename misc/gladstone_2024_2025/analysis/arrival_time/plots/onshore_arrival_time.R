# Show how the arrival times vary from offshroe to onshore
library(terra)
library(parallel)

source_names = list.files('../gladstone_arrival_time_min_and_scenario_average' , full.names = FALSE)

onshore_contour_elev = 0.4
min_onshore_elev = 0.2
offshore_contour_elev = -2.
max_offshore_elev = 0

for (source_name in source_names) {
    gp_arrival_times = Sys.glob(paste0(
        '../gladstone_arrival_time_min_and_scenario_average/',
        source_name,
        '/minimum_arrival_time_domain_*.tif'
    ))
    gp_index = unlist(lapply(strsplit(gp_arrival_times, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))

    # Get arrival time tiles ONLY for high-res domains
    # MIN_DOMAIN_INDEX = 88
    # keep = which(gp_index >= MIN_DOMAIN_INDEX)

    # Specify which domains
    domains_keep = c(90, 94, 98, 103, 104, 105, 106, 107, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 129, 132)
    domains_keep = c(unlist(domains_keep), 59, 60, 65, 66, 71, 77, 78)
    keep = which(gp_index %in% domains_keep)

    elevation_files = Sys.glob('../../jatwc_to_inundation/elevation_in_model/elevation0_domain*.tif')
    elevation_index = unlist(lapply(strsplit(elevation_files, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))

    # Ensure order matches arrival times
    elv_match = match(elevation_index, gp_index)
    elevation_files = elevation_files[elv_match]
    elevation_files = rep(elevation_files, times=8)
    gp_arrival_times = gp_arrival_times[keep]
    gp_index = gp_index[keep]
    elevation_files = elevation_files[keep]


    gp_arrival_times <- gp_arrival_times[!is.na(gp_arrival_times)]
    elevation_files <- elevation_files[!is.na(elevation_files)]

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

    parfun<-function(i){
        # Get arrival times around contours
        offshore = get_arrival_times_on_contour_level(gp_arrival_times[i], elevation_files[i], contour_level=offshore_contour_elev)
        onshore = get_arrival_times_on_contour_level(gp_arrival_times[i], elevation_files[i], contour_level=onshore_contour_elev)
        return(list(offshore=offshore, onshore=onshore))
    }

    all_at = mclapply(1:length(gp_arrival_times), parfun, mc.cores=48)

    offshore_at = do.call(rbind, lapply(all_at, function(x) x$offshore))
    onshore_at = do.call(rbind, lapply(all_at, function(x) x$onshore))

    filename = paste0('offshore_vs_onshore', source_name, '.png')
    png(filename, width=6, height=8, units='in', res=300)
    # To highlight the difference between offshore and onshore, we ensure the 'offshore' points have elevation < 0,
    # while the 'onshore' points have elevation > max_offshore_elev 
    k = which(offshore_at$elev < max_offshore_elev)

    # set xlim
    min_arrival <-min(offshore_at$arrival_time[k]/3600, na.rm=TRUE)
    xlim = c(min_arrival, min_arrival + 3)
    if (source_name %in% c('kurilsjapan', 'southamerica')) {
        xlim = c(min_arrival, min_arrival + 5)
    }

    plot(offshore_at$arrival_time[k]/3600, offshore_at$y[k], pch=19, cex=0.2, 
        xlab='Minimum arrival time (hours)', ylab='Latitude', cex.lab=1.3, cex.axis=1.3, ylim=c(-24.7, -23), xlim=xlim)
    k = which(onshore_at$elev > min_onshore_elev)
    points(onshore_at$arrival_time[k]/3600, onshore_at$y[k], pch=19, cex=0.2, col='red')
    axis(2, at = seq(-24.7, -23, by = .1), labels = FALSE)

    grid(col='orange')
    main = paste0('Minimum Shoreline Arrival: ', source_name)
    title(main=main, cex.main=1.4)

    source('cities_flat.R')
    label_offset = xlim[2]
    text(rep(label_offset, length(cities$lon)), cities$lat, cities$name, pos=2)
    legend('bottom',
        c(
            paste0('Offshore, elevation < ', max_offshore_elev, ' m'),
            paste0('Onshore, elevation > ', min_onshore_elev, ' m')
        ),
        col=c('black', 'red'),
        pch=19,
        cex=1.0,
        bg='white',
        ncol = 2
    )
    dev.off()
}