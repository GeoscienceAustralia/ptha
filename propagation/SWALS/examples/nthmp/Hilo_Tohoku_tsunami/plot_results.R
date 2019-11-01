model_label = commandArgs(trailingOnly=TRUE)[1]
#
# Get the recent model run
#
source('../../../plot.R')
#target_dir = sort(Sys.glob('OUTPUTS/RUN*/RUN*'), decreasing=TRUE)[1]
#x = get_all_recent_results(target_dir)
gauges = merge_multidomain_gauges(multidomain_dir = sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])

# Multi-panel plot
png(paste0('Gauges_plot_', model_label, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(2,1))

#
# Read data boundary
#
input_bc = read.csv('boundary/se_dat_converted.csv')
model_ind = which(gauges$gaugeID == 1)

plot(input_bc, t='l', xlab='Time (s)', ylab='Stage (m)', 
     main='Idealised boundary forcing at target point', cex.main=1.5, 
     cex.lab=1.5, ylim=c(-2,2))
points(gauges$time, gauges$time_var$stage[model_ind,], t='l', col='red')
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
# Time shift (estimated based on first peak arrival)
hilo_time_offset = hilo_detided$time[1] + 50220 - 5700
hilo_detided$time = hilo_detided$time - hilo_time_offset
model_ind = which(gauges$gaugeID == 2)

plot(hilo_detided$time, hilo_detided$stage, t='l',
     main='Hilo harbour gauge [Compare with Fig.7 in Lynett et al. (2017)] ', 
     cex.main=1.5, xlab='Time (s)', ylab='Stage (m)',
     cex.lab=1.5, xlim=range(input_bc[,1]), ylim=c(-2,2))
points(gauges$time, gauges$time_var$stage[model_ind,], t='l', col='red')
grid(col='orange')
legend('topleft', c('Observed (detided)', 'Modelled'), lty=c(1,1), 
       pch=c(NA, NA), col=c('black', 'red'), bty='n')

dev.off()

# Nominal error check on the Hilo gauge -- the model should not agree 'very
# well' because the forcing is idealised.
dat_fun = approxfun(hilo_detided[,1], hilo_detided[,2], rule=2)
gauge_fun = approxfun(gauges$time, gauges$time_var$stage[model_ind,], rule=2)
target_times = seq(5000, 12500, by=10)
err_stat = sum( (dat_fun(target_times) - gauge_fun(target_times))^2 ) / sum( (gauge_fun(target_times) )^2)
if(err_stat < 0.45){
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

plot(ha125$time, ha125$speed, xlim=range(gauges$time), ylim=c(0,1.5),
     xlab='Time (s)', ylab='Speed (m/s)', main='Speed @ HA25 [Compare with Fig.8 in Lynett et al. 2017]',
     cex.main=1.5, cex.lab=1.5)
model_u = gauges$time_var$uh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_v = gauges$time_var$vh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_speed = sqrt(model_u**2 + model_v**2)
points(gauges$time, model_speed, t='l', col='red')
grid(col='orange')

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

plot(ha126$time, ha126$speed, xlim=range(gauges$time), ylim=c(0,1.5), 
     xlab='Time (s)', ylab='Speed (m/s)', main='Speed @ HA26[Compare with Fig.8 in Lynett et al. 2017]',
     cex.main=1.5, cex.lab=1.5)
model_u = gauges$time_var$uh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_v = gauges$time_var$vh[model_ind,]/(gauges$time_var$stage[model_ind,] - gauges$static_var$elevation0[model_ind])
model_speed = sqrt(model_u**2 + model_v**2)
points(gauges$time, model_speed, t='l', col='red')
grid(col='orange')

dev.off()



