# Extraction of tides at specific locations

source('/home/gareth/Code_Experiments/TIDAL_PREDICTION/TPXO9/TPXO9_R_interface/predict_tide.R', chdir=TRUE) # Also located on NCI

sites_to_get = list(
    'Steep_Point_Java2006' =  list(
        # Site coordinates as decimal degrees (longitude, latitude)
        site_coordinates = c(113.18897, -26.15979),

        # Timezone MUST be GMT / UTC. Note that -10 means 'ahead 10 hours', so the
        # time-zone 'Etc/GMT-10' corresponds to eastern Australia. Confusingly, the
        # timezone is often called 'GMT+10'
        start_time = strptime('2006-07-17 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'), # WA time
        end_time   = strptime('2006-07-18 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('2006-07-17 08:19:28', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        ),
    'NW_Cape_Java1994' = list(
        site_coordinates = c(114.03, -21.83),
        start_time = strptime('1994-06-03 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        end_time   = strptime('1994-06-04 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('1994-06-02 18:17:36', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        ),
    'Onslow_Java1994' = list(
        site_coordinates = c(115.12, -21.63),
        start_time = strptime('1994-06-03 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        end_time   = strptime('1994-06-04 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('1994-06-02 18:17:36', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        ),
    'Dampier_Java1994' = list(
        site_coordinates = c(116.75, -20.56),
        start_time = strptime('1994-06-03 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        end_time   = strptime('1994-06-04 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('1994-06-02 18:17:36', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        ),
    'Dampier_Sumba1977' = list(
        site_coordinates = c(116.75, -20.56),
        start_time = strptime('1977-08-19 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        end_time   = strptime('1977-08-20 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('1977-08-19 06:08:55', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        ),
    'PointSamson_Sumba1977' = list(
        site_coordinates = c(117.20764, -20.62827),
        start_time = strptime('1977-08-19 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        end_time   = strptime('1977-08-20 00:00:00', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8'),
        eq_time = strptime('1977-08-19 06:08:55', format = '%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8') + as.difftime(8, units='hours')
        )
    )

for(site_name in names(sites_to_get)){        

    site_coordinates = sites_to_get[[site_name]]$site_coordinates
    start_time = sites_to_get[[site_name]]$start_time
    end_time = sites_to_get[[site_name]]$end_time
    eq_time = sites_to_get[[site_name]]$eq_time

    time_interval = '5 min'
    tide_data = get_tidal_prediction(site_name, site_coordinates, 
        start_time, end_time, time_interval)

    png(paste0('Tides_', site_name, '.png'), width=9, height=5, units='in', res=300)
    plot(tide_data$time, tide_data$tide, t='l', main=paste0('TPXO9 @ ', site_name), 
        sub='Earthquake time in red, hours post earthquake as orange dotted lines.', 
        xlab='Time (WAST)', ylab = 'Stage (m)', cex.main=1.7, cex.axis=1.5, cex.lab=1.5, col='blue')
    #grid(col='orange')
    abline(h=seq(-10, 10, by=0.2), col='orange', lty='dotted')
    abline(h=0, col='orange')

    # Add vertical hourly lines
    for(i in seq(1, 24, by=1)){
        tme = eq_time + as.difftime(i, units='hours')
        points(c(tme, tme), c(-10, 10), col='orange', lty='dotted', t='l')
    }

    points(c(eq_time, eq_time), c(-10, 10), t='l', col='red')

    dev.off()
}
