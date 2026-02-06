source('../../../plot.R')
MSL = 0.78
# Time offset, needed because of our initial condition forcing
TM = 0.0 #-0.6

md_dir = get_recent_output_folder()

gauges = merge_multidomain_gauges(multidomain_dir=md_dir)

gauge_files = paste0('../../nthmp/Conical_shelf_lab/obs_timeseries/WG', gauges$gaugeID, '.txt')

#
# Plot stage
#

png('Stage_gauges.png', width=8, height=5, units='in', res=200)
par(mfrow = c(3,3))
par(mar=c(2,2,1,1))
for(i in 1:9){
    id = gauges$gaugeID[i]
    obs_data = read.table(gauge_files[i], header=FALSE)
    # Convert cm to m
    obs_data[,2] = obs_data[,2] / 100

    model_stage = gauges$time_var$stage[i,] - MSL
    model_time = gauges$time - TM

    plot(obs_data, t='l', xlab='Time (s)', ylab='Stage (m above MSL)',
        xlim=c(0, max(model_time)), ylim=c(-0.1, 0.5))
    abline(h=0, col='orange')
    points(model_time, model_stage, t='l', col='red')
    grid(col='orange')
    title(main=paste0(gsub('.txt', '', basename(gauge_files[i])), ' @ x=', gauges$lon[i], ', y=', gauges$lat[i]), 
          line=-1.5, cex.main=1.7)
    if(i == 1) legend('center', c('Data', 'Model'), col=1:2, lty=c(1,1), pch=c(NA, NA), bty='n', cex=1.6)

    # Rough error checks
    m1 = max(model_stage[model_time < 30])
    o1 = max(obs_data[obs_data[,1] < 30,2])
    #print(c(i, m1, o1, abs(m1-o1), abs(m1-o1)/o1)) 
    if(abs(m1 - o1) < 0.11){
        print('PASS')
    }else{
        print(c('FAIL stage ', i, m1, o1, abs(m1-o1)))
    }
}
dev.off()


#
# Plot speeds
#

png('velocity_gauges.png', width=9, height=7, units='in', res=200)
vel_inds = c(10, 11, 12)
vel_gauge = c('A', 'B', 'C')
par(mfrow=c(3,2))
par(mar=c(4,4,2,1))
for(i in 1:length(vel_inds)){
    vi = vel_inds[i]
    li = vel_gauge[i]
    model_time = gauges$time - TM

    model_u = gauges$time_var$uh[vi,]/(gauges$time_var$stage[vi,]-gauges$static_var$elevation0[vi])
    plot(model_time, model_u, 
         t='l', col='red', xlab='Time (s)', ylab='U-velocity component (m/s)', 
         main=paste0('U-velocity @ site ', li, ', x=', gauges$lon[vi], ', y=', gauges$lat[vi]),
         ylim=c(-2,2), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    u = read.table(paste0('../../nthmp/Conical_shelf_lab/obs_timeseries/U_Velocity_Average', li, '.txt'))
    points(u, t='l', col='black')
    abline(h=0, col='orange')
    grid(col='orange')

    min_time = 0
    if(li == 'B'){
        min_time = 7.5
    }

    # Basic error check
    m1 = max(model_u[model_time < 30 & model_time > min_time])
    o1 = max(u[u[,1] < 30 & u[,1] > min_time, 2], na.rm=TRUE)
    if(abs(m1 - o1) < 0.25*o1){
        print('PASS')
    }else{
        print(c('FAIL vel ', i, m1, o1, abs(m1-o1)/o1))
    }

    model_v = gauges$time_var$vh[vi,]/(gauges$time_var$stage[vi,]-gauges$static_var$elevation0[vi])
    plot(model_time, model_v, 
         t='l', col='red', xlab='Time (s)', ylab='V-velocity component (m/s)', 
         main=paste0('V-velocity @ site ', li, ', x=', gauges$lon[vi], ', y=', gauges$lat[vi]),
         ylim=c(-2,2), cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
    v = read.table(paste0('../../nthmp/Conical_shelf_lab/obs_timeseries/V_Velocity_Average', li, '.txt'))
    points(v, t='l', col='black')
    abline(h=0, col='orange')
    grid(col='orange')
    legend('topright', c('Data', 'Model'), col=1:2, lty=c(1,1), pch=c(NA, NA), bty='n', cex=1.6)
    #print(c(i, 'V', max(model_v[model_time < 30]), max(v[v[,1] < 30, 2], na.rm=TRUE)))

    # Don't do error check on V velocities -- literature suggests no-one can model them too well
    # At y=0 sites we might not expect such large v-velocities (due to the symmetry of the domain),
    # so not surprising that the observations are tricky to model
}
dev.off()

#
# Background plot
#
png('domain_setup.png', width=8.5, height=5, units='in', res=200)
multidomain_image(md_dir, variable='elevation0', time_index=NA, xlim=c(-9, 44), ylim=c(-13,13), zlim=c(-0.1, 1.2),
                  cols=rev(hcl.colors('Spectral', n=200)), use_fields=TRUE)
title(main='Elevation and gauge locations', cex.main=1.5)

points(gauges$lon, gauges$lat, pch=19)
text(gauges$lon[1:9], gauges$lat[1:9], paste0('WG', 1:9), pos=3, cex=1.2)
text(gauges$lon[10:12], gauges$lat[10:12], vel_gauge, pos=1, cex=1.2)

dev.off()
