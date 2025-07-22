file_home = '/home/gareth/Code_Experiments/fortran/Structured_shallow_water/plot.R'
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
source(ifelse(file.exists(file_home), file_home, file_nci))

md_dir = commandArgs(trailingOnly=TRUE)[1]
if(!file.exists(md_dir)) stop(paste0('Could not find file ', md_dir))

# Get the gauge
bb_gauge_RDS = paste0(md_dir, '/bunbury_beacon_gauge.RDS')
if(!file.exists(bb_gauge_RDS)){
    # If the RDS file doesn't exist, then make it
    bb = get_gauges_near_xy(multidomain_dir=md_dir, 
        xy_sites=matrix(c(115.647, -33.294), ncol=2), 
        lonlat_coords=TRUE, verbose=TRUE)

    saveRDS(bb, bb_gauge_RDS)
}else{
    bb = readRDS(bb_gauge_RDS)
}

model_u = bb$time_var$uh[1,]/(bb$time_var$stage[1,] - bb$static_var$elevation0[1])
model_v = bb$time_var$vh[1,]/(bb$time_var$stage[1,] - bb$static_var$elevation0[1])
model_speed = sqrt(bb$time_var$uh[1,]**2 + bb$time_var$vh[1,]**2)/(bb$time_var$stage[1,] - bb$static_var$elevation0[1])
model_time = strptime('2004-12-26 00:58:50', format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT') + as.difftime(bb$time, units='secs')

# Get the data
bb_data = read.csv('../../../field_obs/2004_tsunami_survey_Brian_Gaull/GD_work/BunburyBeaconCurrents/Bunbury_B3_current_data_meterspersecond_degrees_combined.csv', comment.char='#')
bb_data_time = strptime(bb_data$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT-8')
# In the data, bins b1-b5 are similar most of the time, while bins b6-b7 are clearly spurious often (not too surprising, they are at or near the surface). Take the average of b1-b5
bb_data_speed_b1to5 = (bb_data$speed_b1 + bb_data$speed_b2 + bb_data$speed_b3 + bb_data$speed_b4 + bb_data$speed_b5)/5

angle_speed_to_u_v<-function(angle, speed){
    U = sin(angle/180*pi)*speed
    V = cos(angle/180*pi)*speed
    return(list(U=U, V=V))
}

## Quick plot
#plot(model_time, model_speed, t='l', ylim=c(0, 0.6))
#points(bb_data_time, bb_data_speed_b1to5, t='l', col='red')

output_file = paste0('plots/BunburyBeacon_ADCP_Speed_East_North', basename(md_dir), '.png')
png(output_file, width=10, height=7, units='in', res=300)
par(mfrow=c(2,1))
par(mar=c(4.5,4.5,3.5,1))
k = which(bb_data_time > model_time[1] & bb_data_time < model_time[length(model_time)])
# East velocity
plot(bb_data_time[k], angle_speed_to_u_v(bb_data$dir_b1, bb_data$speed_b1)$U[k], t='l', ylim=c(-0.5, 0.5),
    main='Easterly velocity component @ Bunbury Beacon #3\n Model vs data in bins 1-5 (bottom-top)', 
    cex.main=1.5, cex.axis=1.4, cex.lab=1.3,
    xlab='Time (GMT + 8), begins at Sumatra 2004 earthquake, 26/12/2004', ylab='Speed (m/s)')
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b2, bb_data$speed_b2)$U, t='l', col=2)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b3, bb_data$speed_b3)$U, t='l', col=3)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b4, bb_data$speed_b4)$U, t='l', col=4)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b5, bb_data$speed_b5)$U, t='l', col=5)
points(model_time, model_u, t='l', lwd=2)
grid(col='orange')
legend('topleft', c('Model', 'Bin 1', 'Bin 2'), col=c('black', 'black', 'red'), lwd=c(2,1,1), horiz=TRUE, bty='n')
legend('bottomleft', c('Bin 3', 'Bin 4', 'Bin 5'), col=c(3, 4, 5), lwd=c(1,1,1), horiz=TRUE, bty='n')

# North velocity
plot(bb_data_time[k], angle_speed_to_u_v(bb_data$dir_b1, bb_data$speed_b1)$V[k], t='l', ylim=c(-0.5, 0.5),
    main='Northerly velocity component @ Bunbury Beacon #3\n Model vs data in bins 1-5 (bottom-top)', 
    cex.main=1.5, cex.axis=1.4, cex.lab=1.3,
    xlab='Time (GMT + 8), begins at Sumatra 2004 earthquake, 26/12/2004', ylab='Speed (m/s)')
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b2, bb_data$speed_b2)$V, t='l', col=2)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b3, bb_data$speed_b3)$V, t='l', col=3)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b4, bb_data$speed_b4)$V, t='l', col=4)
points(bb_data_time, angle_speed_to_u_v(bb_data$dir_b5, bb_data$speed_b5)$V, t='l', col=5)
points(model_time, model_v, t='l', lwd=2)
grid(col='orange')
legend('topleft', c('Model', 'Bin 1', 'Bin 2'), col=c('black', 'black', 'red'), lwd=c(2,1,1), horiz=TRUE, bty='n')
legend('bottomleft', c('Bin 3', 'Bin 4', 'Bin 5'), col=c(3, 4, 5), lwd=c(1,1,1), horiz=TRUE, bty='n')

dev.off()

# Bunbury Beacon tides
bb_tides = read.csv('../../../field_obs/2004_tsunami_survey_Brian_Gaull/GD_work/BunburyTides/Tide_data_from_Gaull_report_BunburyBeacon3.csv')
bb_tides_time = strptime(bb_tides$time, format='%Y-%m-%d %H:%M:%S', tz='Etc/GMT')
model_stage = bb$time_var$stage[1,]
png('plots/BunburyBeacon_tides_comparison.png', width=9, height=4, units='in', res=300)
plot(model_time + as.difftime(8, units='hours'), model_stage, t='l', col='red', ylim=c(-1,1),
    xlab='Time (GMT + 8)', ylab='Stage residual (m)', main='Bunbury Beacon #3 tide',
    cex.main=1.5, cex.lab=1.3, cex.axis=1.4)
points(bb_tides_time + as.difftime(8, units='hours'), bb_tides$resid, t='l', col='black')
legend('topleft', c('Obs', 'Model'), lty=c(1,1), col=c('black', 'red'), pch=c(NA, NA))
grid(col='orange')
dev.off()
