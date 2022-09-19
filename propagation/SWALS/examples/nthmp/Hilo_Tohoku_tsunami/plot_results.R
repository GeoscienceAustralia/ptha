model_label = commandArgs(trailingOnly=TRUE)[1]
#
# Get the recent model run
#
source('../../../plot.R')

# Just extract the gauges
md_dir = sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1]
gauges = merge_multidomain_gauges(multidomain_dir = md_dir)

# Multi-panel plot
png(paste0('Gauges_plot_', model_label, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(2,1))

#
# Read data boundary
#
input_bc = read.csv('boundary/se_dat_converted.csv')
model_ind = which(gauges$gaugeID == 1)
# At time zero in the boundary data, it is this many hours since the earthquake
# Note we work with data in seconds, but convert to hours for plotting
# This offset gives a visual match to time-series plotted in Lynett et al. (2017)
boundary_hours_offset = 6.6

plot(input_bc[,1]/3600 + boundary_hours_offset, 
     input_bc[,2], t='l', xlab='Time (hours)', ylab='Stage (m)', 
     main='Idealised boundary forcing at target point', cex.main=1.5, 
     cex.lab=1.5, ylim=c(-2,2))
points(gauges$time/3600 + boundary_hours_offset, 
       gauges$time_var$stage[model_ind,], t='l', col='red')
grid(col='orange')
legend('topright', 
       c('Benchmark (based on a model)', 'Modelled (radiation condition --> inexact agreement)'), 
       lty=c(1,1), pch=c(NA, NA), col=c('black', 'red'), bty='n')

# Make some error metric to ensure the BC is OK
input_bc_fun = approxfun(input_bc[,1], input_bc[,2], rule=2)
err_stat = sum( (input_bc_fun(gauges$time) - gauges$time_var$stage[model_ind,])^2 ) / sum( (input_bc_fun(gauges$time) )^2)
if(err_stat < 0.07){
    print('PASS')
}else{
    print('FAIL')
}

#
# Read gauge at Hilo
#
hilo_detided = read.table('test_data/TG_1617760_detided.txt')
names(hilo_detided) = c('time', 'stage')
# The test problem description indicates that a time shift should be applied, estimated based on first peak arrival.
hilo_time_offset = hilo_detided$time[1] + 50220 - 5700
hilo_detided$time = hilo_detided$time - hilo_time_offset
model_ind = which(gauges$gaugeID == 2)

i_top = read.csv('digitize_Lynett2017/Model_mean_free_surface_envelope_top.csv')
i_bottom = read.csv('digitize_Lynett2017/Model_mean_free_surface_envelope_bottom.csv')

plot(hilo_detided$time/3600 + boundary_hours_offset, 
     hilo_detided$stage, t='l',
     main='Hilo harbour gauge \n [Compare with Fig.7 in Lynett et al. (2017)] ', 
     cex.main=1.5, xlab='Time (hours)', ylab='Stage (m)',
     cex.lab=1.5, 
     xlim=range(input_bc[,1]/3600)+boundary_hours_offset, ylim=c(-2.5,2.5))
points(gauges$time/3600 + boundary_hours_offset, 
       gauges$time_var$stage[model_ind,], t='l', col='red')
grid(col='orange')

points(i_top, t='l', col='grey')
points(i_bottom, t='l', col='grey')

legend('topleft', c('Observed (detided)', 'Modelled'), lty=c(1,1), 
       pch=c(NA, NA), col=c('black', 'red'), bty='n')
legend('topright', c('Model intercomparison: average crest/trough envelope'), lty=c(1), 
       pch=c(NA), col=c('grey'), bty='n')

dev.off()

# Nominal error check on the Hilo gauge -- the model should not agree 'very
# well' because the forcing is idealised.
dat_fun = approxfun(hilo_detided$time, hilo_detided$stage, rule=2)
gauge_fun = approxfun(gauges$time, gauges$time_var$stage[model_ind,], rule=2)
target_times = seq(5000, 12500, by=10)
err_stat = sum( (dat_fun(target_times) - gauge_fun(target_times))^2 ) / sum( (gauge_fun(target_times) )^2)
if(err_stat < 0.5){
    print('PASS')
}else{
    print('FAIL')
}


png(paste0('Speed_comparisons_', model_label, '.png'), width=9, height=7, units='in', res=300)

par(mfrow=c(2, 1))
#
# Velocity time-series, HA125
#
ha125 = read.table('test_data/HAI1125_detided_harmonic.txt')
names(ha125) = c('time', 'u', 'v')
# Offset (like for the boundary stage) and convert to seconds.
ha125$time = ha125$time*3600 - hilo_time_offset #
ha125[,2:3] = ha125[,2:3]/100 # Convert velocity to m/s
model_ind = which(gauges$gaugeID == 3)
ha125$speed = sqrt(ha125$u**2 + ha125$v**2)

# Lynett et al (2017) data 
i_top = read.csv('digitize_Lynett2017/Model_mean_speed_envelope_HA25.csv')

plot(ha125$time/3600 + boundary_hours_offset, 
     ha125$speed, 
     xlim=range(gauges$time/3600)+boundary_hours_offset, 
     ylim=c(0,1.8),
     xlab='Time (hours)', ylab='Speed (m/s)', 
     main='Speed (6min sampling average) @ HA25 \n [Compare with Fig.8 in Lynett et al. 2017]',
     cex.main=1.5, cex.lab=1.5)
model_u = gauges$time_var$uh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_v = gauges$time_var$vh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_speed = sqrt(model_u**2 + model_v**2)
points(gauges$time/3600 + boundary_hours_offset, model_speed, t='l', col='red')
#points(gauges$time/3600, model_speed, t='o', pch=19, cex=0.2, col='red')
grid(col='orange')

# Compute smoothed speed
dt = mean(diff(gauges$time)) # Almost constant
N = round(6*60/dt)
if(N%%2 == 0) N = N+1 # Filter should be odd
smooth_speed = filter(model_speed, filter=rep(1, N)/N)
smooth_speed[is.na(smooth_speed)] = model_speed[is.na(smooth_speed)]
points(gauges$time/3600 + boundary_hours_offset, smooth_speed, t='l', col='purple')

points(i_top, t='l', col='grey')

legend('topright', c('Data', 'Model (no averaging)', 'Model (6min average)', 'Model intercomparison: average 6min speed envelope'), 
       pch=c(1, NA, NA, NA), lty=c(NA, 1, 1, 1),
       col=c('black', 'red', 'purple', 'grey'), bty='n')

#
# Velocity time-series, HA126
#
ha126 = read.table('test_data/HAI1126_detided_harmonic.txt')
names(ha126) = c('time', 'u', 'v')
# Offset (like for the boundary stage) and convert to seconds.
ha126$time = ha126$time*3600 - hilo_time_offset #
ha126[,2:3] = ha126[,2:3]/100 # Convert velocity to m/s
model_ind = which(gauges$gaugeID == 4)
ha126$speed = sqrt(ha126$u**2 + ha126$v**2)

# Lynett et al (2017) data 
i_top = read.csv('digitize_Lynett2017/Model_mean_speed_envelope_HA26.csv')

plot(ha126$time/3600 + boundary_hours_offset, 
     ha126$speed, 
     xlim=range(gauges$time/3600) + boundary_hours_offset, 
     ylim=c(0,1.8), 
     xlab='Time (hours)', ylab='Speed (m/s)', 
     main='Speed (6min sampling average) @ HA26 \n [Compare with Fig.8 in Lynett et al. 2017]',
     cex.main=1.5, cex.lab=1.5)
model_u = gauges$time_var$uh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_v = gauges$time_var$vh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_speed = sqrt(model_u**2 + model_v**2)
points(gauges$time/3600 + boundary_hours_offset, model_speed, t='l', col='red')
#points(gauges$time/3600, model_speed, t='o', pch=19, cex=0.2, col='red')
grid(col='orange')

# Compute smoothed speed
dt = mean(diff(gauges$time)) # Almost constant
N = round(6*60/dt)
if(N%%2 == 0) N = N+1 # Filter should be odd
smooth_speed = filter(model_speed, filter=rep(1, N)/N)
smooth_speed[is.na(smooth_speed)] = model_speed[is.na(smooth_speed)]
points(gauges$time/3600 + boundary_hours_offset, smooth_speed, t='l', col='purple')

points(i_top, t='l', col='grey')

legend('topright', c('Data', 'Model (no averaging)', 'Model (6min average)', 'Model intercomparison: average 6min speed envelope'), 
       pch=c(1, NA, NA, NA), lty=c(NA, 1, 1, 1),
       col=c('black', 'red', 'purple', 'grey'), bty='n')

dev.off()

#
# Plot some fields
#

max_speed = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=1, desired_var='max_speed', return_raster=TRUE)
elevation0 = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=1, desired_var='elevation0', return_raster=TRUE)
max_stage = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=1, desired_var='max_stage', return_raster=TRUE)

max_stage[max_stage < elevation0 + 1.0e-03] = NA

# Plot range following Lynett (2017) to avoid eastern boundary
XLIM = c(204.91, 204.95) 
YLIM = c(19.72, 19.75)

png(paste0('Model_fields_', model_label, '.png'), width=12, height=5.1, units='in', res=200)
par(mfrow=c(1,2))
par(mar=c(5,4,3,4))

# Stage max
stagecol = rev(hcl.colors(255, 'spectral'))
plot(max_stage, col=stagecol, xlim=XLIM, ylim=YLIM, main='Maximum stage (m) near the harbour',
    xlab='Lon', ylab='Lat')

# Speed max
suppressMessages(library(fields))
speedcol = tim.colors(255) #rev(hcl.colors(255, 'spectral'))
plot(min(max_speed, max_speed*0 + 3), col=speedcol, xlim=XLIM, ylim=YLIM, main='Maximum speed (m/s) near the harbour',
    xlab='Lon', ylab='Lat', zlim=c(0, 3)) # Limit to 3m/s like in Lynett (2017)

dev.off()

# Elevation and domain partition
model_gauge_sites = list(
    'Hilo tide station' = rev(c( 19.7308, 204.9447 )),
    'HA25' = rev(c(19.7452, 204.9180)),
    'HA26' = rev(c(19.7417, 204.9300)),
    'Synthetic boundary forcing' = rev(c(19.7576, 204.93)))

png(paste0('Model_elevation_', model_label, '.png'), width=6, height=5.1, units='in', res=200)
elevcol = rev(hcl.colors(255, 'Earth'))
plot(elevation0, col=elevcol, 
    # xlim=XLIM, ylim=YLIM, 
    # main='Elevation (m) near the harbour',
    # xlim=XLIM, ylim=YLIM, 
    main='Elevation (m) and domain partitioning',
    xlab='Lon', ylab='Lat')
contour(elevation0, level=0, add=TRUE, col='blue')

# Add bounding boxes (without merging)
bbox = get_domain_interior_bbox_in_multidomain(md_dir)
for(i in 1:length(bbox$domain_interior_bbox)){
  bb = bbox$domain_interior_bbox[[i]]
  points(rbind(bb, bb[1,]), t='l', col='red', lty='dotted')
}

for(i in 1:length(model_gauge_sites)){
    xy = model_gauge_sites[[i]]
    points(xy[1], xy[2], col='black', pch=19, cex=0.5)
    text(xy[1], xy[2], names(model_gauge_sites)[i], adj=c(0,1), col='black')
}

dev.off()

