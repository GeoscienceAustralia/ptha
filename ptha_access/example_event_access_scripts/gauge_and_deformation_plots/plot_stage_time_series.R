# Get the PTHA access codes
source('plot_parameters.R', local=TRUE)

desired_rows = stage_plot_desired_rows 

kt2_events = ptha$get_source_zone_events_data(source_zone=source_zone, 
    slip_type='stochastic', desired_event_rows = desired_rows)


#
# Wave time series
#
wave_ts = vector(mode='list', length=nrow(wave_sites))

for(i in 1:length(wave_ts)){

    target_point = wave_sites[i,]
    wave_ts[[i]] = ptha$get_flow_time_series_at_hazard_point(
        kt2_events, target_points = target_point, event_ID=1:nrow(kt2_events$events))
}

if(!is.null(dart_data_site1)){
    dt1 = read.csv(dart_data_site1, stringsAsFactors=FALSE)

    dt1_tme = difftime(
        strptime(dt1$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'),
        strptime(event_time_GMT, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'), 
        units='secs')
}

if(!is.null(dart_data_site2)){
    dt2 = read.csv(dart_data_site2, stringsAsFactors=FALSE)
    dt2_tme = difftime(
        strptime(dt2$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'),
        strptime(event_time_GMT, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT'), 
        units='secs')
}

#
# Plot the waves
#
pdf('waves_check.pdf', width=9, height=8)
for(i in 1:nrow(kt2_events$events)){

    layout(matrix(c(1, 1, 3, 2, 2, 4), nrow=2, byrow=TRUE))

    par(mar = c(5.1, 4.1, 4.1, 2.1))
    plot(wave_ts[[1]]$time, wave_ts[[1]]$flow[[1]][i,,1], t='l', ylim=c(-1,1)*site1_yrange, xlim=c(0, 15000) + site1_tstart)
    grid()
    max_slip = round(max(as.numeric(strsplit(wave_ts[[1]]$events$event_slip_string[i], '_')[[1]])))
    mw = round(wave_ts[[1]]$events$Mw[i], 1)
    vary_mu_mw = round(wave_ts[[1]]$events$variable_mu_Mw[i], 1)
    is_possible = ifelse(wave_ts[[1]]$events$rate_annual[i] > 0, ' Possible ', ' Impossible! ')
    title_words = paste0('Site 1; ', is_possible, 
        ' ; , Peak slip = ', max_slip, '\n event-row ', desired_rows[i],
        '; Mw=', mw, '; vary_mu_Mw = ', vary_mu_mw)
    title(title_words)

    location_title = paste0(c('lon', 'lat', 'elev', 'ID'), '=', round(wave_ts[[1]]$location, 2), collapse=" ")
    mtext(location_title, side=1, line=-2)
    if(!is.null(dart_data_site1)){
        points(dt1_tme, dt1$resid, t='l', col='blue', lty='dashed')
    }
    
    # Add second wave site
    plot(wave_ts[[2]]$time, wave_ts[[2]]$flow[[1]][i,,1], t='l', ylim=c(-1,1)*site2_yrange, xlim=c(0, 15000) + site2_tstart)
    grid()
    #title('DART in Pacific Ocean east of Tonga')
    title('Site 2')
    location_title = paste0(c('lon', 'lat', 'elev', 'ID'), '=', round(wave_ts[[2]]$location, 2), collapse=" ")
    mtext(location_title, side=1, line=-2)
    if(!is.null(dart_data_site2)){
        points(dt2_tme, dt2$resid, t='l', col='blue', lty='dashed')
    }


    # If the event is possible then add the raster
    my_raster = Sys.glob(paste0(output_dir, '/stochastic_*', desired_rows[i], '.png'))
    if(length(my_raster) == 1){
        rast = brick(my_raster)
        ext0 = extent(rast)
        ext1 = ext0
        ext1@ymax = 600
        ext2 = ext0
        ext2@ymin = 600
        plotRGB(rast, ext=ext1)
        title(basename(my_raster))
        plotRGB(rast, ext=ext2)
        title(basename(my_raster))
    }else{
        plot(c(0,1), c(0,1), col='white', main='Impossible (zero rate!)')
        plot(c(0,1), c(0,1), col='white', main='Impossible (zero rate!)')
    }
}
dev.off()

#
# Save the session to keep the data handy
# 
mytime = as.character(julian(Sys.time(), units='secs'))
save.image(paste0('stage_time_series_plot_', mytime, '.RData'))

