#
# Plot model time-series at Neds beach, Lord Howe Island, for the 2021 New Hebrides tsunami
#

md_dir = commandArgs(trailingOnly=TRUE)[1]

# Get the SWALS routines
file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

# Gauges at Neds beach, Lord Howe Island
target_gauge_locs = matrix(c(
    159.0658, -31.51774,
    159.0650, -31.51791,
    159.0636, -31.51592), byrow=TRUE, ncol=2)

gauges = get_gauges_near_xy(md_dir, target_gauge_locs, lonlat_coords=TRUE, verbose=TRUE)

event_start_time = strptime('2021-02-10 13:19:55', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

# Bureau observer (Amy) observed the waves at Ned's beach from 3.45am LHI time (4.45pm GMT time)
observer_start_time = strptime('2021-02-10 16:45:00', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')

model_time = as.difftime(gauges$time, units='secs') + event_start_time
k = which(gauges$time < 6*3600)

png(paste0('Neds_beach_timeseries_', basename(dirname(md_dir)), '-', basename(md_dir), '.png'),
    width=10, height=5, units='in', res=200)
plot(model_time[k], gauges$time_var$stage[1,k], t='l', 
    ylim=c(-1.2, 1.2),
    main="Modelled tsunami for 3 sites at Ned's beach, Lord Howe Is. \n 2021-02-10 Loyalty Islands tsunami",
    xlab='GMT Time (add 11 hrs to get time at Lord Howe)', 
    ylab='Modelled tsunami water level, no tides (m)', 
    cex.lab=1.5, cex.main=1.5, cex.axis=1.5)
points(model_time[k], gauges$time_var$stage[2,k], t='l', col='green')
points(model_time[k], gauges$time_var$stage[3,k], t='l', col='purple')

abline(v=as.difftime(0, units='secs') + observer_start_time, col='red', lwd=3) # For some reason the 'as.difftime' part is needed
abline(h=seq(-2,2,by=0.25), lty='dotted', col='orange')
abline(v=as.difftime(seq(0, 40)*3600, units='secs') + event_start_time, lty='dotted', col='orange')
abline(h=0.4, lty='dotted')
title(sub='Bureau observer saw 40cm waves (above regular sea level) during a 15min period, starting 16:45 GMT (red line)')
dev.off()

