# Append to plot filenames
model_name = commandArgs(trailingOnly=TRUE)[1]

# Gage stage time-series -- gauges 1-6
gauge_h = read.table(file='tjhr_a_9521830_sup_0001/building_gauges_h.txt', skip=2)
obs_h = list(time = gauge_h[,1], stage = gauge_h[,2:7])

# Gauge uv -- only gauges 1-5
gauge_uv = read.delim(file='tjhr_a_9521830_sup_0001/building_gauges_uv.txt', skip=2)
obs_uv = list(time=gauge_uv[,1], u=gauge_uv[,seq(2, 10, by=2)], v=gauge_uv[,seq(3, 11, by=2)])

source('../../plot.R')
gauges = merge_multidomain_gauges(multidomain_dir=sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])

model_metadata = lapply(Sys.glob(paste0(gauges$multidomain_dir, '/RUN*')), 
                        f<-function(x) get_all_recent_results(x,read_grids=FALSE, read_gauges=FALSE))

# Stage at 6 gauges
png(paste0('Gauges_stage_plot_', model_name, '.png'), width=9, height=6, units='in', res=300)
par(mfrow=c(3,2))
par(mar=c(4, 4, 2, 1))
for(i in 1:6){
    if(i != 6){
        YLIM = c(0, 0.15)
    }else{
        YLIM = c(0, 0.4)
    }
    plot(obs_h$time, obs_h$stage[,i], t='l', main=paste0('Stage @ Gauge ', i), cex.main=1.5,
         xlab='Time (s)', ylab='Stage (m)', ylim=YLIM)
    points(gauges$time, gauges$time_var$stage[i,],t='l',col='red')
    grid(col='orange')

    # SOME CRUDE TESTS. Because observations and models for this problem are quite variable,
    # we make some ad-hoc decisions about what to test. 
    if(i == 1){
        # Mean stage between 15-25s ~ 0.08m
        inds = which(gauges$time > 15 & gauges$time < 25)
        if(mean(gauges$time_var$stage[i,inds] - 0.08) < 0.01){
            print('PASS')
        }else{
            print('FAIL')
        }
    }else if(i == 2){
        # Mean stage between 20-30s ~ 0.1m
        inds = which(gauges$time > 20 & gauges$time < 30)
        if(mean(gauges$time_var$stage[i,inds] - 0.1) < 0.01){
            print('PASS')
        }else{
            print('FAIL')
        }

    }else if(i == 4){
        # Mean stage between 15-20s ~ 0.08m
        inds = which(gauges$time > 15 & gauges$time < 20)
        if(mean(gauges$time_var$stage[i,inds] - 0.08) < 0.01){
            print('PASS')
        }else{
            print('FAIL')
        }
    }else if(i == 6){
        # Mean stage between 18-23s ~ 0.2m
        inds = which(gauges$time > 18 & gauges$time < 23)
        if(mean(gauges$time_var$stage[i,inds] - 0.2) < 0.01){
            print('PASS')
        }else{
            print('FAIL')
        }
    }

}
dev.off()

# U-velocity. This is the main flow direction. No velocity for G6
png(paste0('Gauges_Uvel_plot_', model_name, '.png'), width=9, height=6, units='in', res=300)
par(mfrow=c(3,2))
par(mar=c(4, 4, 2, 1))
for(i in 1:5){
    YLIM = c(-0.2, 1)*2.5
    plot(obs_uv$time, obs_uv$u[,i], t='l', main=paste0('U-vel @ Gauge ', i), cex.main=1.5,
         xlab='Time (s)', ylab='U-vel (m)', ylim=YLIM)
    model_u = gauges$time_var$uh[i,]/(gauges$time_var$stage[i,]-gauges$static_var$elevation0[i])
    points(gauges$time, model_u,t='l',col='red')
    grid(col='orange')
}
dev.off()

# V-velocity. This is the main flow direction. No velocity for G6
png(paste0('Gauges_Vvel_plot_', model_name, '.png'), width=9, height=6, units='in', res=300)
par(mfrow=c(3,2))
par(mar=c(4, 4, 2, 1))
for(i in 1:5){
    YLIM = c(-1, 1)*1.0
    plot(obs_uv$time, obs_uv$v[,i], t='l', main=paste0('V-vel @ Gauge ', i), cex.main=1.5,
         xlab='Time (s)', ylab='V-vel (m)', ylim=YLIM)
    model_v = gauges$time_var$vh[i,]/(gauges$time_var$stage[i,]-gauges$static_var$elevation0[i])
    points(gauges$time, model_v,t='l',col='red')
    grid(col='orange')
}
dev.off()


#
# Velocity snapshots
#
vel_times = c(1, 3, 5, 10, 15)
vel_times_char = c('01', '03', '05', '10', '15')
vel_files = paste0('tjhr_a_9521830_sup_0001/building_vel_t', vel_times_char, '.txt')
zero_vel_threshold = 0.0 # Ignore smaller velocities
PLOT_XLIM = c(8, 14)
all_polys = lapply(Sys.glob('poly/*.csv'), read.csv)
for(i in 1:length(vel_times)){

    png(paste0('velocity_field_t', vel_times_char[i],'_', model_name, '.png'),
        width=6, height=8, res=300, units='in')
    par(mfrow=c(2, 1))
    par(mar=c(3, 3, 3, 1))
    vel_scale = 20

    # The measurements have x/y origin being 6.75, 1.8 relative to our model
    vel_measured = read.table(vel_files[i], skip=2)
    vel_measured[,1] = vel_measured[,1] + 6.75 + 0.8
    vel_measured[,2] = vel_measured[,2] + 1.8
    zero_vels = ((vel_measured[,3]**2 + vel_measured[,4]**2) > zero_vel_threshold**2)

    plot(vel_measured[,1:2], asp=1, col=0, xlim=PLOT_XLIM, 
         main=paste0('VELOCITY MEASURED @ t', vel_times_char[i]),
         cex.main=1.5, xlab='x', ylab='y')
    arrows(vel_measured[,1], vel_measured[,2], 
           vel_measured[,1] + vel_measured[,3]/vel_scale * zero_vels, 
           vel_measured[,2] + vel_measured[,4]/vel_scale * zero_vels,
           length=0, col=zero_vels)
    for(j in 1:length(all_polys)) polygon(all_polys[[j]][,1], all_polys[[j]][,2], fill=NA, border='red')

    # Compare with model -- get UH, VH, and depth(=stage because elevation=0)
    model_time_ind = which.min(abs(model_metadata[[1]]$time - vel_times[i]))
    model_uh = merge_domains_nc_grids(multidomain_dir=gauges$multidomain_dir, 
        domain_index=1, desired_var='uh', desired_time_index=model_time_ind)
    model_vh = merge_domains_nc_grids(multidomain_dir=gauges$multidomain_dir, 
        domain_index=1, desired_var='vh', desired_time_index=model_time_ind)
    model_h = merge_domains_nc_grids(multidomain_dir=gauges$multidomain_dir, 
        domain_index=1, desired_var='stage', desired_time_index=model_time_ind)

    model_u = model_uh$uh/(model_h$stage + 1.0e-10)
    model_v = model_vh$vh/(model_h$stage + 1.0e-10)

    nearest_x_ind = sapply(vel_measured[,1], f<-function(x) which.min(abs(model_h$xs - x)))
    nearest_y_ind = sapply(vel_measured[,2], f<-function(x) which.min(abs(model_h$ys - x)))

    plot(model_h$xs[nearest_x_ind], model_h$ys[nearest_y_ind], col=0, asp=1, xlim=PLOT_XLIM, 
         main=paste0('VELOCITY MODELLED @ t', vel_times_char[i]),
         cex.main=1.5, xlab='x', ylab='y')

    mx = model_h$xs[nearest_x_ind]
    my = model_h$ys[nearest_y_ind]
    mu = mapply(f<-function(ix, iy) model_u[ix, iy], nearest_x_ind, nearest_y_ind)
    mv = mapply(f<-function(ix, iy) model_v[ix, iy], nearest_x_ind, nearest_y_ind)
    zero_vels = ((mu*mu + mv*mv) > zero_vel_threshold**2)
    #arrows(mx, my, mx + mu*zero_vels/vel_scale, my + mv*zero_vels/vel_scale, length=0)
    arrows(vel_measured[,1], vel_measured[,2], 
           vel_measured[,1] + mu*zero_vels/vel_scale, 
           vel_measured[,2] + mv*zero_vels/vel_scale, length=0, col=zero_vels)
    for(j in 1:length(all_polys)) polygon(all_polys[[j]][,1], all_polys[[j]][,2], fill=NA, border='red')
    dev.off()
}

