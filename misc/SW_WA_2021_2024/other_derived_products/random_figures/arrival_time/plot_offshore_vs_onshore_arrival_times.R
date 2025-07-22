# Show how the arrival times vary from offshroe to onshore
library(terra)

# Get arrival time tiles ONLY for high-res domains
MIN_DOMAIN_INDEX = 293
gp_arrival_times = Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/arrival_time/sunda2/arrival_time_minimum/*.tif')
gp_index = unlist(lapply(strsplit(gp_arrival_times, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))
keep = which(gp_index >= MIN_DOMAIN_INDEX)

elevation_files = Sys.glob('../../WA_tsunami_modelling_2021_2024_definitive_outputs/model_outputs/greater_perth/elevation_in_model/*.tif')
elevation_index = unlist(lapply(strsplit(elevation_files, '_domain_'), function(x) as.numeric(gsub(".tif", "", x[2]))))

# Ensure order matches arrival times
elv_match = match(elevation_index, gp_index)
elevation_files = elevation_files[elv_match]
elevation_index = elevation_index[elv_match]
gp_arrival_times = gp_arrival_times[keep]
gp_index = gp_index[keep]
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
    # Get arrival times around contours of -2m, and +1m
    offshore = get_arrival_times_on_contour_level(gp_arrival_times[i], elevation_files[i], contour_level=-2)
    onshore = get_arrival_times_on_contour_level(gp_arrival_times[i], elevation_files[i], contour_level=1)
    return(list(offshore=offshore, onshore=onshore))
}

all_at = mclapply(1:length(gp_arrival_times), parfun, mc.cores=16)

offshore_at = do.call(rbind, lapply(all_at, function(x) x$offshore))
onshore_at = do.call(rbind, lapply(all_at, function(x) x$onshore))

png('offshore_vs_onshore.png', width=6, height=8, units='in', res=300)
# To highlight the difference between offshore and onshore, we ensure the 'offshore' points have elevation < 0,
# while the 'onshore' points have elevation > 0.8. 
k = which(offshore_at$elev < 0)
plot(offshore_at$arrival_time[k]/3600, offshore_at$y[k], pch=19, cex=0.2, xlim=c(2, 5), ylim=c(-34, -28), 
    xlab='Minimum arrival time (hours)', ylab='Latitude', cex.lab=1.3, cex.axis=1.3)
k = which(onshore_at$elev > 0.8)
points(onshore_at$arrival_time[k]/3600, onshore_at$y[k], pch=19, cex=0.2, col='red')

grid(col='orange')
title(main='Minimum arrival times near the shoreline, sunda2 \n (includes mainland, islands & estuaries)', cex.main=1.4)

source('cities.R')
text(rep(2.5, length(cities$lon)), cities$lat, cities$name)
legend('topleft', c('Offshore, elevation < 0 m', 'Onshore, elevation > 0.8 m'), col=c('black', 'red'), pch=19, cex=1.3, bg='white')
dev.off()


