library(rptha)

MAKE_PDF_PLOT = FALSE # Try to keep file size down.

maximum_elevation_of_gauges = -1 # Avoid 'dry' model points

# Get swals plotting code
swals = new.env()
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R', local=swals)

# Get gauge data interface
gd = new.env()
source('/g/data/w85/tsunami/DATA/TIDES/CODE_INTERFACE/gauge_data_links.R', local=gd, chdir=TRUE)

# Read the event data where it exists
event_data = gd$get_data_for_event(event_id)

# Find the tidal gauges
target_gauges = matrix(unlist(lapply(event_data, f<-function(x) x$coord)), ncol=2, byrow=TRUE)
target_gauge_domains = swals$find_domain_containing_point(target_gauges, md=NULL, multidomain_dir=multidomain_dir)

if(any(target_gauge_domains$number_of_possible_domains > 1)){
    print('Warning: detected gauge with more than one candidate domain, but we do not search all options')
}

# Model start time in days (based on earthquake time)
model_start_julian = julian(model_start_time)

if(MAKE_PDF_PLOT){
    gauges_plot_name = paste0('gauges_plot_', event_id, '_', basename(dirname(multidomain_dir)), '.pdf')
    pdf(paste0(multidomain_dir, '/', gauges_plot_name), width=12, height=11)
}

outlist = vector(mode='list', length=length(event_data))
names(outlist) = names(event_data)

# Plot each gauge
for(i in 1:length(event_data)){

    # Get the gauges on the domain
    gauges = swals$get_gauges(target_gauge_domains$domain_dir[i])

    # Move on if gauge is not on domain
    if(is(gauges, 'try-error')) next

    # Limit the search to gauges that satisfy the elevation criteria
    ki = which(gauges$static_var$elevation0 < maximum_elevation_of_gauges)

    # Find index nearest to the target gauge
    gi = ki[ which.min(distHaversine(
                cbind(gauges$lon[ki], gauges$lat[ki]), 
                matrix(c(target_gauges[i,1:2]), ncol=2, nrow=length(gauges$lon[ki]), byrow=TRUE))
        ) ]

    # Record the distance, to be careful that we are comparing data with a nearby model point.
    distance_from_model_to_gauge = distHaversine(
        cbind(gauges$lon[gi], gauges$lat[gi]), 
        matrix(c(target_gauges[i,1:2]), ncol=2, nrow=1))

    # Store a bunch of outputs
    outlist[[i]] = list()
    outlist[[i]]$event_data = event_data[[i]]
    #outlist[[i]]$gauges = gauges # All gauges, in case we want to check neighbours
    outlist[[i]]$gauges_dir = target_gauge_domains$domain_dir[i]
    outlist[[i]]$gi = gi
    outlist[[i]]$model_time = gauges$time
    outlist[[i]]$model_stage = gauges$time_var$stage[gi,]
    outlist[[i]]$model_gaugeloc = c(gauges$lon[gi], gauges$lat[gi], gauges$gaugeID[gi])
    outlist[[i]]$model_start_julian = model_start_julian
    outlist[[i]]$distance_from_model_to_gauge = distance_from_model_to_gauge


    if(MAKE_PDF_PLOT){

        par(mfrow=c(3,1))
        par(mar=c(4,2,4,1))

        # Set first plot y-limits, with some 'trimming' of the observations in case there are spikes
        YLIM = range(
            c(range(gauges$time_var$stage[gi,])*1.1, 
            quantile(event_data[[i]]$obs$resid, c(0.01, 0.99), na.rm=TRUE)*1.2 + 
                0.15*c(-1,1)*diff(range(gauges$time_var$stage[gi,]))))

        # First plot -- modelled tsunami vs de-tided obs
        plot(model_start_julian + gauges$time/(24*3600), 
            gauges$time_var$stage[gi,], t='l', col='red',
            xlab='Time (days)', ylab='Tsunami (m)',
            main=paste0(event_id, ': ', names(event_data)[i]), 
            cex.main=1.8, cex.lab=1.5,
            ylim=YLIM)
        points(event_data[[i]]$obs$juliant, event_data[[i]]$obs$resid, t='p', pch=19, cex=0.2)
        abline(h=0, col='orange')
        grid(col='orange')
        #title(main=paste0(event_id, ': ', names(event_data)[i]), cex.main=1.8)

        # Second plot -- observations alone -- if we have them
        have_raw_obs = !all(is.na(event_data[[i]]$obs$height))
        if(have_raw_obs){
            plot(event_data[[i]]$obs$juliant, event_data[[i]]$obs$height, t='o', pch=19, cex=0.2,
                xlab='Time (days)', ylab='Obs tide (m)', main = 'Observations only', 
                xlim=range(model_start_julian + gauges$time/(24*3600), na.rm=TRUE) , 
                cex.lab=1.5, cex.main=1.8 )
            abline(h=0, col='orange')
            grid(col='orange')
        }else{
            plot(event_data[[i]]$obs$juliant, event_data[[i]]$obs$resid, t='o', pch=19, cex=0.2,
                xlab='Time (days)', ylab='Obs tide (m)', main = '')
            title('Do not have raw data for this one')
            grid(col='orange')
        }

        # Third plot -- model vs obs for first 6 hours
        k = which(abs(gauges$time_var$stage[gi,]) > 1.0e-03)[1]
        if(is.na(k) | length(k) == 0){
            plot(c(0, 1), c(0, 1), main = 'No stage values > 1mm, skipping', col='white')
        }else{
            t0 = model_start_julian + gauges$time[k]/(3600*24)
            t1 = t0 + 1/24*12
            plot(model_start_julian + gauges$time/(24*3600), 
                gauges$time_var$stage[gi,], t='l', col='red',
                xlab='Time (days)', ylab='Tsunami (m)',
                main=paste0('First 12 hours from model arrival'), 
                cex.main=1.8, cex.lab=1.5,
                ylim=YLIM, xlim=c(t0, t1))
            points(event_data[[i]]$obs$juliant, event_data[[i]]$obs$resid, t='o', pch=19, cex=0.2)
            abline(h=0, col='orange')
            grid(col='orange')
        }
    }

}

if(MAKE_PDF_PLOT) dev.off()

saveRDS(outlist, 
    file=paste0(multidomain_dir, '/', 'gauges_plot_', event_id, '_', 
                basename(dirname(multidomain_dir)), '.RDS'))

#
# Export some additional offshore gauges
#
 
# Find a desired set of tidal gauges
target_gauges = as.matrix(read.csv('../extra_gauges_to_extract.csv', header=TRUE)) 
target_gauge_domains = swals$find_domain_containing_point(target_gauges, md=NULL, multidomain_dir=multidomain_dir)
outlist = vector(mode='list', length=nrow(target_gauges))
for(i in 1:nrow(target_gauges)){

    # Get the gauges on the domain
    gauges = swals$get_gauges(target_gauge_domains$domain_dir[i])

    # Move on if gauge is not on domain
    if(is(gauges, 'try-error')) next

    # Find index nearest to the target gauge
    gi = which.min(distHaversine(
        cbind(gauges$lon, gauges$lat), 
        matrix(c(target_gauges[i,1:2]), ncol=2, nrow=length(gauges$lon), byrow=TRUE))
        )

    outlist[[i]] = list()

    outlist[[i]]$gauges_dir = target_gauge_domains$domain_dir[i]
    outlist[[i]]$gi = gi
    outlist[[i]]$model_time = gauges$time
    outlist[[i]]$model_stage = gauges$time_var$stage[gi,]
    outlist[[i]]$model_gaugeloc = c(gauges$lon[gi], gauges$lat[gi], gauges$gaugeID[gi])
    outlist[[i]]$model_start_julian = model_start_julian
}


saveRDS(outlist, 
    file=paste0(multidomain_dir, '/', 'extra_gauges_to_extract_', event_id, '_', 
            basename(dirname(multidomain_dir)), '.RDS'))
