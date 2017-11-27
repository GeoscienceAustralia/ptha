#
# Single station stage exceedance rate computation
#
# This is actually a slow computational method ('quick' means 'quick to code'), but is useful
# to cross check the other results

lon = 151.41
lat = -34.08

tsunami_files = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_*.nc')

library(rptha)

# Get hazard points -- faster to not use the '_tsunami' file
fid = nc_open(tsunami_files[1], readunlim=FALSE)
n = length(fid$dim$station$vals)
hp = data.frame(
    lon     = rep(NA, n),
    lat     = rep(NA, n),
    elev    = rep(NA, n),
    gaugeID = rep(NA, n)
    )
for(i in 1:n){
    hp$lon[i] = ncvar_get(fid, 'lon', start=i, count=1)
    hp$lat[i] = ncvar_get(fid, 'lat', start=i, count=1)
    hp$elev[i] = ncvar_get(fid, 'elev', start=i, count=1)
    hp$gaugeID[i] = ncvar_get(fid, 'gaugeID', start=i, count=1)
    if(i%%100 == 1) print(i)
}
nc_close(fid)

# Find index of point nearest to lon/lat
ni = lonlat_nearest_neighbours(cbind(lon, lat), cbind(hp$lon, hp$lat))

# Get stage and rates, for each source

stage_rate = vector(mode='list', length=length(tsunami_files))
names(stage_rate) = basename(dirname(dirname(tsunami_files)))
for(i in 1:length(tsunami_files)){
    print(basename(tsunami_files[i]))
    fid = nc_open(tsunami_files[i], readunlim=FALSE)

    event_rate = ncvar_get(fid, 'event_rate_annual')
    event_rate_upper = ncvar_get(fid, 'event_rate_annual_upper_ci')
    event_rate_lower = ncvar_get(fid, 'event_rate_annual_lower_ci')
    peak_stage = ncvar_get(fid, 'max_stage', start=c(1, ni), count=c(-1,1))
    event_Mw = round(ncvar_get(fid, 'event_Mw'), 3) # Deal with floating point imperfections in netcdf
    site = rep(basename(dirname(dirname(tsunami_files[i]))), length=length(event_rate))
    row_index = 1:length(event_rate)

    stage_rate[[i]] = data.frame(
        event_rate = event_rate,
        event_rate_upper = event_rate_upper,
        event_rate_lower = event_rate_lower,
        peak_stage = peak_stage,
        site = site,
        row_index=row_index,
        event_Mw = event_Mw)

    nc_close(fid)
}

# Back-calculate the stage-vs-rate curves
stage_rate_all = do.call(rbind, stage_rate)

odr = rev(order(stage_rate_all$peak_stage))

stg = stage_rate_all$peak_stage[odr]
er = cumsum(stage_rate_all$event_rate[odr])
er_up = cumsum(stage_rate_all$event_rate_upper[odr])
er_lo = cumsum(stage_rate_all$event_rate_lower[odr])

#
# Plot the data
#
plot(stg, er, t='l', log='xy', xlim=c(0.01, max(stg)))
grid()
points(stg, er_up, t='l', col='red')
points(stg, er_lo, t='l', col='red')
title(paste0('Stage-vs-exceedance rate @ (', round(lon,2), ', ', 
    round(lat, 2), 
    ') \n Comparison of file values (points) and separate calculation (lines)'))

# Compare with the values in the file
fid = nc_open('tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc', 
    readunlim=FALSE)
stages = fid$dim$stage$vals
ers = ncvar_get(fid, 'stochastic_slip_rate'            , start=c(1,ni), count=c(-1,1))
ers_up = ncvar_get(fid, 'stochastic_slip_rate_upper_ci', start=c(1,ni), count=c(-1,1))
ers_lo = ncvar_get(fid, 'stochastic_slip_rate_lower_ci', start=c(1,ni), count=c(-1,1))
nc_close(fid)

points(stages, ers, pch=19, cex=1.0, col='brown')
points(stages, ers_up, pch=19, cex=1.0, col='pink')
points(stages, ers_lo, pch=19, cex=1.0, col='pink')


#
# Function to examine the distribution of earthquake magnitudes
#
peak_stage_magnitude_summary<-function(stage_threshold, source_zone){

    is_site = as.character(stage_rate_all$site) == source_zone
    
    k = which( (stage_rate_all$peak_stage > stage_threshold) & is_site)

    if(length(k) == 0){
        output = NA
    }else{
        output = aggregate(stage_rate_all$event_rate[k], 
            by=list(stage_rate_all$event_Mw[k]), 
            f<-function(x) c(sum(x), length(x)))
        output = data.frame(Mw=output[,1], rate_exceeding=output$x[,1], n=output$x[,2])
    
        # Compute overall numbers of Mw events
        fracs = rep(NA, length(output[,1]))
        for(i in 1:length(output[,1])){
            n = sum(is_site & (stage_rate_all$event_Mw == output[i,1]))
            fracs[i] = output[i,3]/n
        }
    
        output = cbind(output, data.frame(fraction_events = fracs))
    }

    return(output)
}

plot_deaggregation_summary<-function(stage_threshold){

    k = which( (stage_rate_all$peak_stage > stage_threshold) & (stage_rate_all$event_rate > 0))

    if(length(k) == 0){
        plot(c(0, 1), c(0, 1), 
            main=paste0('No events exceeding stage_threshold = ', stage_threshold))
    }else{

        par(mfrow=c(2,2))

        rate_by_source = aggregate(stage_rate_all$event_rate[k], 
            by=list(source_zone=as.character(stage_rate_all$site[k])), sum)

        dotchart(rate_by_source$x, labels=rate_by_source[,1], 
            main=paste0('Rate of events from each source-zone with \n stage exceeding ', stage_threshold),
            xlab='Rate (events/year)')

        # Also plot rate-vs-Mw for the 3 largest contributors 
        m1 = order(rate_by_source$x, decreasing=TRUE)
        for(i in 1:min(length(m1), 3)){
            # Name of source-zone
            sz = rate_by_source[m1[i], 1]
            rate_by_Mw = peak_stage_magnitude_summary(stage_threshold, sz) 
            dotchart(rate_by_Mw[,2], 
                labels=paste0(rate_by_Mw[,1], ' (', round(rate_by_Mw$fraction_events, 3), ')'),
                xlab='Rate (events/year)', ylab='')
            mtext(side=2, 'Magnitude (fraction that exceed)', line=2.3)
            title(paste0(sz, ': Exceedance rates vs magnitude'))
        }
    }
}

png(paste0('deagg_summary_', lon, '_', lat, '.png'), width=10, height=7, res=200, units='in')

plot_deaggregation_summary(1.5)

dev.off()

