library(geosphere)
# SWALS plot codes
swals = new.env()
#source('~/Code_Experiments/fortran/Structured_shallow_water/plot.R', local=swals)
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R', local=swals)

# Interface for our tide-gauge data
gauges = new.env()
source('../../gauges/nukualofa/get_gauge_data_for_event.R', local=gauges, chdir=TRUE)


historic_event_gauge_plot<-function(md_dir){

    # The event name should be at the start of the multidomain dirname
    event_name = strsplit(basename(dirname(md_dir)), '_')[[1]][1]

    # Check it matches with an event name in the gauge data
    allowed_event_names = names(gauges$gauge_files)
    if(!any(event_name == allowed_event_names)){
        stop(paste0('Could not match event_name: ', event_name))
    }

    event_data = gauges$get_gauge_data(event_name)

    # Add a swals$get_gauge_nearest_point(lonlat, md_dir)
    target_gauge_domains = swals$find_domain_containing_point(
        matrix(event_data$coord, ncol=2), 
        md=NULL, multidomain_dir=md_dir)

    gauges = swals$get_gauges(target_gauge_domains$domain_dir)

    if(event_data$gauge_name == 'nukualofa'){
        # Find a gauge near to the desired site, but also in the ocean. With
        # coarse models we can end up having the nearest gauge on land, so better
        # to move it offshore
        distance_metric = distHaversine(
            cbind(gauges$lon, gauges$lat),
            cbind(gauges$lon*0 + event_data$coord[1], gauges$lat*0+event_data$coord[2]))
        distance_metric[gauges$static_var$elev > 0] = 9e+09 # Dry land is a big penalty
        gi = which.min(distance_metric)
    }

    model_time = event_data$start_time + as.difftime(gauges$time, units='secs')
    model_stage = gauges$time_var$stage[gi,]
    k = which(event_data$obs$time >= model_time[1] &
              event_data$obs$time <= model_time[length(model_time)])
    plot_ylim = range(c(range(model_stage), range(event_data$obs$resid[k])))

    plot(model_time, gauges$time_var$stage[gi,],t='l', col='red', ylim=plot_ylim,
         xlab='Time (UTC)', ylab='Stage (m)')
    points(event_data$obs$time, event_data$obs$resid, t='o', pch=19, cex=0.3)
    title(paste0(event_name, " @ Nuku'alofa \n"), cex.main=1.8)
    title(sub=md_dir, cex.sub=0.7)
    abline(h=seq(-2,2,by=0.1), col='orange', lty='dotted')
    abline(h=0, col='orange', lty='solid')

}
