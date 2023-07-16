#
# Plotting routine that we reuse in scripts for both tides and pressure
#
library(geosphere)
source('global_variables.R')
source('spectral_highpass_filter.R')

# Plot gauge data for 14th - 19th January. This works for both sea-level data, and pressure data.
three_panel_plot<-function(gauge_data, site_name, add_highpass_filtered_series=TRUE, filter_threshold=1/(3*3600), site_lon=NA, site_lat=NA){

    plot_start = strptime('2022-01-14 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    plot_end   = strptime('2022-01-19 00:00:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    event_start = explosion_start_time #strptime('2022-01-15 04:15:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')

    # Can plot either stage, or pressure
    if('stage' %in% names(gauge_data)){
        data_s = gauge_data$stage
        plot_y_label = 'Stage (m)'
        final_plot_ylim = c(-1.2, 1.2)
    }else if('pressure' %in% names(gauge_data)){
        data_s = gauge_data$pressure
        plot_y_label = 'MSL Pressure (hPa)'
        final_plot_ylim = c(-4, 4)
    }else{
        stop('No stage or pressure inputs')
    }

    # Find data indices overlapping plot_start/plot_end
    i0 = max(c(which(gauge_data$juliant < julian(plot_start) & !is.na(data_s)), 1))
    i1 = min(c(which(gauge_data$juliant > julian(plot_end  ) & !is.na(data_s)), length(gauge_data$juliant)))
   
    # Indices for 1 day from the event start 
    j0 = max(c(which(gauge_data$juliant < (julian(event_start)  ) & !is.na(data_s)), 1))
    j1 = min(c(which(gauge_data$juliant > (julian(event_start)+1) & !is.na(data_s)), length(gauge_data$juliant)))

    # 2 panel plot, if using high-pass filter
    if(add_highpass_filtered_series) {
        par(mfrow=c(3,1))
        par(mar=c(4,4,2,1))
    }

    if(!is.finite(i0) | !is.finite(i1)){
        # Something has failed. Make a plot anyway, to keep the plot ordering correct.
        plot(c(0, 1), c(0, 1))
        title(paste0('Read failed for ', site_name))
        if(add_highpass_filtered_series){
            plot(c(0, 1), c(0, 1))
            plot(c(0, 1), c(0, 1))
        }
    }else{
        # Typical case
        ii = i0:i1 # Plot indices
        plot(gauge_data$time[ii], data_s[ii], t='l', xlab='Time', ylab=plot_y_label)
        title(site_name, cex.main=1.5)
        grid(col='orange')
        points(c(event_start, event_start), c(-9e+09, 9e+09), col='red', t='l')
        # Hourly vertical lines
        abline(v=(event_start + as.difftime(seq(1, 200), units='hours')), col='green', lty='dotted')

        if(add_highpass_filtered_series){
            data_t = as.numeric(gauge_data$juliant - gauge_data$juliant[1])*3600*24
            filt = try(spectral_highpass_filter(data_t, data_s, interp_dt = 15, cutoff_frequency=filter_threshold))
            if(!is(filt, 'try-error')){
                # Plot de-tided series
                plot(gauge_data$time[ii], filt$highfreq[ii], t='l', xlab='Time', ylab='Residual')
                points(c(event_start, event_start), c(-9e+09, 9e+09), col='red', t='l')
                # Hourly vertical lines
                abline(v=(event_start + as.difftime(seq(1, 200), units='hours')), col='green', lty='dotted')

                # Include the theoretical lamb-wave arrivals
                lamb_wave_speed = LAMB_WAVE_SPEED #316 # m/s -- this is approximate
                earth_radius = EARTH_RADIUS #6378137 # m -- consistent with distHaversine
                dist_tonga = distHaversine(matrix(c(site_lon, site_lat), ncol=2), matrix(hunga_tonga_volcano, ncol=2))
                # Expected arrival times of wave travelling the short path from Tonga to the site
                wave_1 = ( dist_tonga + c(0, 1, 2, 3, 4) * 2 * pi * earth_radius)/lamb_wave_speed
                # Expected arrival times of wave travelling the long path from Tonga to the site
                wave_2 = ( (2*pi*earth_radius - dist_tonga) + c(0, 1, 2, 3, 4) * 2 * pi * earth_radius)/lamb_wave_speed

                # Include rough theoretical lamb-wave arrivals
                abline( v=(explosion_start_time + as.difftime(wave_1,units='secs')), col='blue')
                abline( v=(explosion_start_time + as.difftime(wave_2,units='secs')), col='brown')

                # Plot only first 24 hours
                jj = j0:j1
                plot(gauge_data$time[jj], filt$highfreq[jj], t='l', xlab='Time', ylab='Residual',
                    ylim=final_plot_ylim)
                points(c(event_start, event_start), c(-9e+09, 9e+09), col='red', t='l')
                # Hourly vertical lines
                abline(v=(event_start + as.difftime(seq(1, 200), units='hours')), col='green', lty='dotted')

                # Include rough theoretical lamb-wave arrivals
                abline( v=(explosion_start_time + as.difftime(wave_1,units='secs')), col='blue')
                abline( v=(explosion_start_time + as.difftime(wave_2,units='secs')), col='brown')
                
            }else{
                # If something fails, make a plot anyway, to keep the plot ordering correct.
                plot(c(0, 1), c(0, 1), main='error')
                plot(c(0, 1), c(0, 1), main='error')
            }
        }
    }
}
       

