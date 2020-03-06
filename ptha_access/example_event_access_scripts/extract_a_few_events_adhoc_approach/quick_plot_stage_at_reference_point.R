library(rptha)
ptha = new.env()
ptha_access_script_path = '../../get_PTHA_results.R'

if(!file.exists(ptha_access_script_path)) stop('Please edit ptha_access_script_path to point to the correct script location')

source(ptha_access_script_path, local=ptha, chdir=TRUE)

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
gauges = data.frame(ncvar_get(fid, 'lon'), ncvar_get(fid, 'lat'), ncvar_get(fid, 'elev'), ncvar_get(fid, 'gaugeID'))
nc_close(fid)
#
par(mfrow=c(2,1))
return_period_info = ptha$get_stage_exceedance_rate_curve_at_hazard_point(
    target_point=reference_point, make_plot=TRUE)

#
# Location of gauges
#
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

    total_source_contrib = udr[i]
    other_source_contrib = 1 - sum(source_contrib$x)
    
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
