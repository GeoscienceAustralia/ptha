# Plot the submerged island experiment
source('../../../plot.R')

x = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])

fix_time_order<-function(input_df){
    # The observations sometimes have times with order reversal. 
    # Make sure the times are ordered
    time_order = order(input_df[,1])
    return(input_df[time_order,])
}

#
# Plot the last timestep to show the model layout
#
suppressMessages(library(fields))
nt = dim(x[[1]]$ud)[3]
png('Model_elevation_and_speed.png', width=12, height=5.4, units='in', res=200)
par(mfrow=c(3,1))
par(mar = c(2, 3.3, 2, 2))
image.plot(x[[1]]$xs, x[[1]]$ys, x[[1]]$elev0, asp=1, zlim=c(-0.06, 0.01), las=1, 
           main='Flume elevation with conical island. The flow is from left to right.', 
           legend.args=list('Elevation (m)', side=2), cex.main=2.,
           xlab="", ylab="")
mtext('x (m)', side=1, line=1)
mtext('y (m)', side=2, line=2.3)

speed = sqrt(x[[1]]$ud[,,nt]**2 + x[[1]]$vd[,,nt]**2)/pmax(x[[1]]$stage[,,nt] - x[[1]]$elev0, 1e-06)
image.plot(x[[1]]$xs, x[[1]]$ys, speed, asp=1, las=1, 
           main='Flow speed at the final model timestep. Eddies develop downstream of the island.', 
           legend.args=list('Speed (m/s)', side=2), cex.main=2.,
           xlab="", ylab="")
mtext('x (m)', side=1, line=1)
mtext('y (m)', side=2, line=2.3)

image.plot(x[[1]]$xs, x[[1]]$ys, x[[1]]$stage[,,nt], asp=1, las=1, zlim=c(-1,1)*0.003, 
           main='Flow stage at the final model timestep. Eddies develop downstream of the island.', 
           legend.args=list('Stage (m)', side=2), cex.main=2.,
           xlab="", ylab="")
mtext('x (m)', side=1, line=1)
mtext('y (m)', side=2, line=2.3)

dev.off()


# Plot gauges
sites = vector(mode='list', length=2)
sites[[1]] = list()
sites[[1]]$u_vel = fix_time_order(read.table('obs/SL_S1_U.DAT'))
sites[[1]]$v_vel = fix_time_order(read.table('obs/SL_S1_V.DAT'))
sites[[2]] = list()
sites[[2]]$u_vel = fix_time_order(read.table('obs/SL_S2_U.DAT'))
sites[[2]]$v_vel = fix_time_order(read.table('obs/SL_S2_V.DAT'))

png('Velocities_at_2_sites.png', width=4, height=6, units='in', res=300)
par(mfrow=c(4,1))
par(mar=c(3,3,2,1))
for(gi in 1:2){
    model_time = x[[1]]$gauges$time
    model_elev = x[[1]]$gauges$static_var$elevation0[gi]
    model_stage = x[[1]]$gauges$time_var$stage[gi,]
    model_uh = x[[1]]$gauges$time_var$uh[gi,]
    model_vh = x[[1]]$gauges$time_var$vh[gi,]

    model_depth = model_stage - model_elev

    # Only plot the last 100s of model time (for comparison with data in
    # "statistically-steady-state-turbulent" regime)
    toffset = max(model_time) - 100

    #
    # Plot u-velocity
    #

    plot(sites[[gi]]$u_vel, t='l', ylim=c(0, 0.2), xlab='Time (s)', ylab='Speed (m/s)')
    grid(col='orange')
    abline(h=0, col='orange', lty='dashed', lwd=2)
    points(model_time - toffset, model_uh/model_depth, t='l', col='red')
    title(main=paste0('Site ', gi, ': u-velocity'))

    # Compare model and data -- mean
    k = which(model_time > toffset)
    model_mean_u = mean(model_uh[k]/model_depth[k])
    obs_mean_u = weighted.mean(sites[[gi]]$u_vel[,2], 
        w=0.5*(c(0, diff(sites[[gi]]$u_vel[,1])) + c(diff(sites[[gi]]$u_vel[,1]), 0)))
    if(abs(model_mean_u - obs_mean_u) < max(0.006, 0.2*abs(obs_mean_u))){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Compare model and data -- sd
    model_sd_u = sd(model_uh[k]/model_depth[k])
    obs_sd_u = sqrt(weighted.mean((sites[[gi]]$u_vel[,2] - obs_mean_u)**2, 
        w=0.5*(c(0, diff(sites[[gi]]$u_vel[,1])) + c(diff(sites[[gi]]$u_vel[,1]), 0))))
    if(abs(model_sd_u - obs_sd_u) < max(0.006, 0.2*obs_sd_u)){
        print('PASS')
    }else{
        print('FAIL')
    }


    #
    # Plot v-velocity
    #

    plot(sites[[gi]]$v_vel, t='l', ylim=c(-0.1, 0.1))
    grid(col='orange')
    abline(h=0, col='orange', lty='dashed', lwd=2)
    points(model_time - toffset, model_vh/model_depth, t='l', col='red')
    title(main=paste0('Site ', gi, ': v-velocity'))

    # Compare model and data -- mean
    k = which(model_time > toffset)
    model_mean_v = mean(model_vh[k]/model_depth[k])
    obs_mean_v = weighted.mean(sites[[gi]]$v_vel[,2], 
        w=0.5*(c(0, diff(sites[[gi]]$v_vel[,1])) + c(diff(sites[[gi]]$v_vel[,1]), 0)))
    if(abs(model_mean_v - obs_mean_v) < max(0.006, 0.2*abs(obs_mean_v))){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Compare model and data -- sd
    model_sd_v = sd(model_vh[k]/model_depth[k])
    obs_sd_v = sqrt(weighted.mean((sites[[gi]]$v_vel[,2] - obs_mean_v)**2, 
        w=0.5*(c(0, diff(sites[[gi]]$v_vel[,1])) + c(diff(sites[[gi]]$v_vel[,1]), 0))))
    if(abs(model_sd_v - obs_sd_v) < max(0.008, 0.2*abs(obs_sd_v))){
        print('PASS')
    }else{
        print('FAIL')
    }

}
dev.off()
