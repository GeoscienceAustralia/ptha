library(rptha)


# INPUTS
all_nc = Sys.glob('event_time_series/*.nc')
reference_point = NULL #c(153.9422, -27.2858)
# END INPUTS


# Get coastline data for plot
if(file.exists("GADM_2.8_AUS_adm1.rds")){
    aus_coastline = readRDS('GADM_2.8_AUS_adm1.rds')
}else{
    aus_coastline = try(getData('GADM', country='AUS', level=1))
}


pdf('README_summary_plot.pdf', width=12, height=12)

#
# Plot stage vs return period at the reference point
#
fid = nc_open(all_nc[1])

# Read the reference point in the file. 
txt = ncatt_get(fid, varid=0, 'reference_point_coordinates_where_stage_exceedance_rates_were_computed')
reference_point_in_file = as.numeric(scan(text=txt$value, sep=","))
if(is.null(reference_point)[1]) reference_point = reference_point_in_file

rp_stages = ncvar_get(fid, 'rp_curve_stage')
rp_rate = ncvar_get(fid, 'rp_curve_rate')
rp_rate_lower = ncvar_get(fid, 'rp_curve_rate_lower_ci')
rp_rate_upper = ncvar_get(fid, 'rp_curve_rate_upper_ci')
gauges = data.frame(ncvar_get(fid, 'lon'), ncvar_get(fid, 'lat'), ncvar_get(fid, 'elev'), ncvar_get(fid, 'gaugeID'))

nc_close(fid)

par(mfrow=c(2,1))
plot(rp_stages, rp_rate, log='xy', t='o', 
    xlab=paste0(' Peak stage'),
    ylab='Exceedance rate (events/year)',
    ylim=c(1.0e-06, 1))
points(rp_stages, rp_rate_upper, t='l', col='red')
points(rp_stages, rp_rate_lower, t='l', col='red')
title(paste0('Peak tsunami stage vs exceedance rate (due to all earthquake source-zones in this PTHA) \n', 
    round(reference_point_in_file[1],4), ', ', round(reference_point_in_file[2],4)))
abline(h=10**(seq(-6, 0)), col='brown', lty='dotted')
abline(v=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20), 
    col='brown', lty='dotted')
legend('topright', c('Mean rate (over all logic tree branches)', '95% credible intervals'),
    col=c('black', 'red'), lty=c(1,1), pch=c(1,NA), bg='white')

#
# Location of gauges
#
#plot(wrld_simpl, xlim=range(gauges[,1]), ylim=range(gauges[,2]), asp=1, 
#    axes=TRUE, col='grey', xlab='lon', ylab='lat')
plot(gauges[,1], gauges[,2], col='red', asp=1, xlab='lon', ylab='lat')
if(class(aus_coastline) != 'try-error'){
    plot(aus_coastline, add=TRUE, col=rgb(0, 0, 0, alpha=0.2))
}
grid(col='brown')
title(main='Locations of gauge time-series (red) and reference point where time-series are output (green)')
points(reference_point[1], reference_point[2], col='green', pch=19)

#
# Show which sources matter for each exceedance rate
#
sourcename = c()
desired_rate = c()
proportion_rate = c()
for(i in 1:length(all_nc)){
    fid = nc_open(all_nc[i])
    sourcename = c(sourcename, strsplit(basename(all_nc[i]), '_tsunami_')[[1]][1])
    desired_rate = c(desired_rate, ncvar_get(fid, 'desired_exceedance_rate'))
    proportion_rate = c(proportion_rate, ncvar_get(fid, 'proportion_rate_from_this_source'))
    nc_close(fid)
}
udr = unique(desired_rate)
n0 = ceiling(sqrt(length(udr)))
par(mfrow=c(n0, ceiling(length(udr)/n0)))
for(i in 1:length(udr)){
    k = which(desired_rate == udr[i])
    source_contrib = aggregate(proportion_rate[k], by=list(sourcename[k]), FUN=mean)

    total_source_contrib = udr
    other_source_contrib = udr - sum(source_contrib$x)
    
    plot_names = c(source_contrib[,1], 'other')
    plot_weight = c(source_contrib[,2], max(0, other_source_contrib))

    pie(plot_weight, labels=plot_names, border='black')
    title(main = paste0('Source-zone contribution to hazard at rate = 1/', round(1/udr[i], 3)))
}




#
# Plot gauge time-series at the reference point
#
par(mfrow=c(6,1))
par(mar=c(1,2,3,1))
for(i in 1:length(all_nc)){

    fid = nc_open(all_nc[i], readunlim=FALSE)

    gauges = data.frame(ncvar_get(fid, 'lon'), ncvar_get(fid, 'lat'), ncvar_get(fid, 'elev'), ncvar_get(fid, 'gaugeID'))
    gauge_ind = lonlat_nearest_neighbours(cbind(reference_point[1], reference_point[2]), gauges[,1:2])
    distance_km = distHaversine(reference_point, gauges[gauge_ind,1:2])/1000
    print(paste0('Distance between reference point and nearest point (km) = ', distance_km))

    time = ncvar_get(fid, 'time')
    stage = ncvar_get(fid, 'stage')
    expected_peak_stage = ncvar_get(fid, 'event_peak_stage_at_reference_point')
   
    plot(time, stage[,gauge_ind], t='l')
    abline(v=seq(0, 200000, by=3600), col='grey')
    abline(h=expected_peak_stage, col='red')
    title(basename(all_nc[i]))

    nc_close(fid)
}
dev.off()
