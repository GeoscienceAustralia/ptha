#
# This code processes some data for Nuku'alofa tide-gauge from late 2014 to 2020.
# It plots a monthly series of "maximum-sea-level - mean-sea-level", and also extracts
# the data around the 2015 Mw 8.3 Chile tsunami event, which is de-tided and written to
# a separate file.
#

parse_BOM_2014_2020<-function(){
    x = read.csv('./tg_1min.csv')
    times = strptime(x[,1], format='%Y%m%dT%H%M%SZ', tz='Etc/UTC')
    juliant = julian(times)
    output = data.frame(time=times, juliant=juliant, stage=x[,2], std_dev=x[,3])
    return(output)
}

BOM_2014_2020 = parse_BOM_2014_2020()

monthly_mean_sl = aggregate(BOM_2014_2020$stage, list(month=format(BOM_2014_2020$time, '%Y-%m')), mean)
monthly_max_sl = aggregate(BOM_2014_2020$stage, list(month=format(BOM_2014_2020$time, '%Y-%m')), max)
stopifnot(all(monthly_mean_sl[,1] == monthly_max_sl[,1])) # Check months line up

# Plot monthly max - monthly mean
png('monthly_sea_level_maxima_above_mean.png', width=8, height=6, units='in', res=300)
plot(monthly_max_sl$x - monthly_mean_sl$x, xlab='Month (2014-2020 dataset)', 
     ylab='Sea level: Monthly "max - mean"', 
     main=paste0('Monthly sea level maxima (above monthly MSL) \n Median = ', 
                 median(monthly_max_sl$x - monthly_mean_sl$x)))
abline(h=0.8, col='red')
dev.off()

output_data = data.frame(month=monthly_mean_sl[,1], mean_sea_level=monthly_mean_sl$x, max_sea_level=monthly_max_sl$x)
write.csv(output_data, 'monthly_sea_levels_summary.csv', row.names=FALSE)


#
# Get a de-tided record around the 2015 Chile tsunami
#
target_dates = strptime(c('2015-09-12 00:00:00', '2015-09-26 00:00:00'), 
                        format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
inds = which(BOM_2014_2020$time >= target_dates[1] &
             BOM_2014_2020$time <= target_dates[2])

nuku_chile2015 = BOM_2014_2020[inds,]
source('../spectral_highpass_filter.R')
nuku_chile2015_filt = spectral_highpass_filter(
    as.numeric(nuku_chile2015$juliant - nuku_chile2015$juliant[1])*3600*24,
    nuku_chile2015$stage)
output = data.frame(time=nuku_chile2015$time, juliant=nuku_chile2015$juliant,
                    stage=nuku_chile2015$stage, std_dev=nuku_chile2015$std_dev,
                    resid=nuku_chile2015_filt$highfreq)
write.csv(output, file='nukualofa_chile2015_detided.csv', row.names=FALSE)

png('nukualofa_chile2015.png', width=14, height=4, units='in', res=300)
target_dates = strptime(c('2015-09-17 00:00:00', '2015-09-20 00:00:00'), 
                        format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
k = which(output$time > target_dates[1] & output$time < target_dates[2])
plot(output$time[k], output$resid[k], t='l', xlab='Time', ylab='Stage residual (m)')
eq_time = strptime('2015/09/16 22:54:32', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC')
points(c(eq_time, eq_time), c(-10, 10), t='l', col='red')
grid()
title("Nuku'alofa Chile 2015 tsunami")
dev.off()


make_plots<-function(){

    #
    # Clear record for the 2015 Mw 8.3 Chile tsunami
    #
    png('Chile2015_earthquake_tsunami.png', width=10, height=4, units='in', res=300)
    inds = 548000:555000
    eq_time = strptime('2015/09/16 22:54:32', format='%Y/%m/%d %H:%M:%S', tz='Etc/UTC')
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$stage[inds], t='l',
         main = 'Tsunami from Mw 8.3 Chile \n Earthquake time: 2015/09/16 22:54:32 UTC')
    points(c(eq_time, eq_time), c(-10, 10), t='l', col='red')
    dev.off()


    #
    # Cyclone Gita -- impact shows up more in the standard-deviation column
    #
    png('Cyclone_Gita_2018.png', width=10, height=7, units='in', res=300)
    start_end_time = strptime(c('2018-02-01 00:00:00', '2018-02-20 00:00:00'), 
                              format='%Y-%m-%d %H:%M:%S', tz='Etc/UTC')
    inds = which(BOM_2014_2020$time > start_end_time[1] & BOM_2014_2020$time < start_end_time[2])
    par(mfrow=c(2,1))
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$stage[inds], t='l',
         main = 'Feb 2018 including Cyclone Gita: Stage ', 
         xlab='Time', ylab='Stage (m)')
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$std_dev[inds], t='l',
         main = 'Feb 2018 including Cyclone Gita: Standard Deviation', 
         xlab='Time', ylab='Standard deviation (m)')
    dev.off()

    #
    # Cyclone Harold 2020!
    #
    png('Cyclone_Harold_2019.png', width=10, height=7, units='in', res=300)
    inds = 2930000:2960000
    par(mfrow=c(2,1))
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$stage[inds], t='l',
         main = 'Cyclone Harold, Early April 2020: Stage ', 
         xlab='Time', ylab='Stage (m)')
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$std_dev[inds], t='l',
         main = 'Cyclone Harold, Early April 2020: Standard Deviation', 
         xlab='Time', ylab='Standard deviation (m)')
    dev.off()

    # This inundation certainly caused damage
    # See https://en.wikipedia.org/wiki/Cyclone_Harold
    wikipedia_text = "Power outages began affecting parts of Tonga due to falling trees caused by the storm on April 9.[154] The center of Harold passed 90–100 km (55–60 mi) south of Tongatapu, lashing Tonga with heavy rains and wind; a peak gust of 80 km/h (50 mph) was registered at 'Eua Airport.[155][156] Damage to food crops and water supplies occurred in 'Eua and Tongatapu.[157] Storm surge, reaching 0.86 m (2 ft 10 in) above king tide,[155] inundated coastal extents of Tongatapu,[158] with their greatest impacts on the island's central and western shores.[157] Three tourist resorts west of Nuku'alofa were destroyed; their beach-side cottages, events complexes, and residences were razed by the surge.[154][159][160] Of the islands, 'Eua was most badly affected, with serious damage wrought to its wharf. Some houses were unroofed and electricity was lost throughout the island.[158] Casualties were reported in the kingdom on April 10, although cut communications by the storm prevented confirmation of them.[161] Farther inland, vegetation and crops were damaged by the storm.[160] On April 23, Tonga's Minister of Finance revealed that the total Damages from Cyclone Harold in Tonga is estimated to in excess of US$111 million.[162]"

    #
    # Full record -- useful to see the tide-range, for example
    #
    png('Full_record.png', width=10, height=7, units='in', res=300)
    inds = 1:nrow(BOM_2014_2020)
    par(mfrow=c(2,1))
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$stage[inds], t='l',
         main = "Nuku'alofa 2014-2020: Stage ", 
         xlab='Time', ylab='Stage (m)')
    grid(col='orange')
    plot(BOM_2014_2020$time[inds], BOM_2014_2020$std_dev[inds], t='l',
         main = "Nuku'alofa 2014-2020: Standard Deviation", 
         xlab='Time', ylab='Standard deviation (m)',
         ylim=c(0, max(BOM_2014_2020$std_dev[inds])) )
    grid(col='orange')
    dev.off()
}



