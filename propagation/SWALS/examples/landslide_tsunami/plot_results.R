# Plot code
source('../../plot.R')

# Get the most recent 002 domain (nested domain)
x = get_all_recent_results(sort(Sys.glob('OUTPUTS/RUN_*/RUN_*0002_*'), decreasing=TRUE)[1])

# Analytical shoreline solution
shore = read.csv('DATA/Shoreline.csv')

ERR_TOL = 0.05

# Plot instants in time
ts = c(160, 175, 220)
yi = ceiling(length(x$ys)/2) # Use this y-index for model result
vel_ylims = list(c(-10, 15), c(-16, 10), c(-4, 1.5)) # Plot limits
stage_ylims = list(NULL, NULL, c(0, 20))
for(i in 1:length(ts)){

    png_filename = paste0('Snapshot_near_shoreline_time_', ts[i], '.png')
    png(png_filename, width=8, height=8, units='in', res=300)

    # Get the analytical solution from a file
    m1 = read.csv(paste0('DATA/t', round(ts[i]), '.csv'))
    par(mfrow=c(2,1))
    # This "analytical solution" can go below the bed in some regions -- find them
    beach_slope = -0.1
    nullpts = which(m1[,2] <  beach_slope * m1[,1])
    #if(length(nullpts) > 0){
    #    # Fix stage
    #    m1[nullpts,2] = beach_slope * m1[nullpts,1] 
    #    # Fix speed
    #    m1[nullpts,3] = NA
    #}
    #analytical_shoreline = m1[max(which(m1[,2] < (beach_slope * m1[,1]))),1]
    analytical_shoreline = approx(shore[,1], shore[,2], xout=ts[i])$y

    # Plot free-surface
    plot(m1[,1], m1[,2], t='p', 
        #main=paste0('Stage at time=', round(ts[i]), 's', 
        #    ' \n Note the analytical solution can go below bed (suggesting it is not exact)'), 
        main=paste0('Stage at time=', round(ts[i]), 's'),
        xlab='x', ylab='Stage (m)',
        ylim=stage_ylims[[i]])
    ind = which.min(abs(ts[i] - x$time))
    points(x$xs, x$stage[,yi,ind],t='o', col='red', pch=19, cex=0.2)
    points(x$xs, x$elev0[,yi],t='l', col='blue')
    abline(v=analytical_shoreline, col='brown', lty='dashed')

    legend('right', c('Analytical', 'Numerical', 'Bed', 'Analytical shoreline'), 
           col=c('black', 'red', 'blue', 'brown'),
           pch=c(1, NA, NA, NA), lty=c(NA, 'solid', 'solid', 'dashed'), 
           bty='n', cex=1.3)


    # Plot velocity
    plot(m1[,1], m1[,3], t='p', 
         main=paste0('Velocity at time=', round(ts[i]), 's'), 
         xlab='x', ylab='Velocity (m/s)', ylim=vel_ylims[[i]])
    points(x$xs, x$ud[,yi,ind]/(x$stage[,yi,ind]-x$elev0[,yi] + 1e-100), t='o', 
           col='red', pch=19, cex=0.2)

    abline(v=analytical_shoreline, col='brown', lty='dashed')

    # Test 1 -- stage (PASS/FAIL)
    model_at_analytical_x = approx(x$xs, x$stage[,yi, ind], xout=m1[,1])
    err_stat = mean(abs(m1[,2] - model_at_analytical_x$y), na.rm=TRUE)/diff(range(m1[,2]))
    if(err_stat < ERR_TOL){
        print(c('PASS'))
    }else{
        print(c('FAIL', err_stat))
    }

    # Test 2 -- velocity (PASS/FAIL)
    model_at_analytical_x = approx(x$xs, x$ud[,yi, ind]/(x$stage[,yi,ind] -x$elev0[,yi] + 1e-100), xout=m1[,1])
    err_stat = mean(abs(m1[,3] - model_at_analytical_x$y), na.rm=TRUE)/diff(range(m1[,2]))
    if(err_stat < ERR_TOL){
        print(c('PASS'))
    }else{
        print(c('FAIL', err_stat))
    }
    dev.off()
}

# Plot "shoreline", defined as the site with smallest x-coordinate that has
# depth exceeding "shore_threshold"
shore = read.csv('DATA/Shoreline.csv')
shore_threshold = 0.1
shore_x = x$time * 0
shore_vel = shore_x
for(i in 1:length(shore_x)){
    ind = min(which(x$stage[,yi,i] - x$elev0[,yi] > shore_threshold))
    shore_x[i] = x$xs[ind]
    shore_vel[i] = x$ud[ind,yi,i]/(x$stage[ind,yi,i] - x$elev0[ind,yi])
}
png_filename = 'Shoreline_timeseries.png'
png(png_filename, width=8, height=8, units='in', res=300)
par(mfrow=c(2,1))
plot(shore[,1], shore[,2], t='l', 
     main=paste0('Shoreline \n (Numerical shoreline defined by depth > ', shore_threshold, 'm)'), 
     xlab='Time', ylab='Position', lwd=3)
points(x$time, shore_x, t='o', col='red', pch=19, cex=0.2)
grid()
legend('topleft', c('Analytical', 'Numerical'), col=c('black', 'red'), lty=c(1,1), pch=c(NA, NA), 
       bty='n', cex=1.3)
plot(shore[,1], shore[,3], t='l', main='Shoreline velocity', xlab='Time', ylab = 'Velocity (m/s)',
    ylim=c(-20, 15), lwd=3)
points(x$time, shore_vel, t='o', col='red', pch=19, cex=0.2)
grid()
dev.off()
