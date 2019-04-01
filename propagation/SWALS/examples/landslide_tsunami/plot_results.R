# Plot code
source('../../plot.R')
# Get the most recent 002 domain
x = get_all_recent_results(sort(Sys.glob('OUTPUTS/RUN_*/RUN_*0002_*'), decreasing=TRUE)[1])

ERR_TOL = 0.05

pdf('benchmark_plots.pdf', width=8, height=8)

# Plot instants in time
ts = c(160, 175, 220)
yi = ceiling(length(x$ys)/2) # Use this y-index for model result
vel_ylims = list(c(0, 15), c(-16, 10), c(-4, 1.5))
stage_ylims = list(NULL, NULL, c(0, 20))
for(i in 1:length(ts)){
    m1 = read.csv(paste0('DATA/t', round(ts[i]), '.csv'))
    #m1 = matrix(scan(paste0('DATA/t', round(ts[i]), '.csv'), sep=",", skip=1, quiet=TRUE), ncol=3, byrow=TRUE)
    par(mfrow=c(2,1))
    plot(m1[,1], m1[,2], t='l', 
        main=paste0('Stage ', round(ts[i]), 
            ' \n Note analytical solution can go below bed, suggesting it is not exact'), 
        xlab='x', ylab='Stage (m)',
        ylim=stage_ylims[[i]], lwd=3)
    ind = which.min(abs(ts[i] - x$time))
    points(x$xs, x$stage[,yi,ind],t='o', col='red', pch=19, cex=0.2)
    points(x$xs, x$elev0[,yi],t='l', col='blue')

    plot(m1[,1], m1[,3], t='l', main=paste0('Velocity ', round(ts[i])), xlab='x', ylab='Velocity (m/s)',
        ylim=vel_ylims[[i]], lwd=3)
    points(x$xs, x$ud[,yi,ind]/(x$stage[,yi,ind]-x$elev0[,yi]), t='o', col='red', pch=19, cex=0.2)

    model_at_analytical_x = approx(x$xs, x$stage[,yi, ind], xout=m1[,1])

    # Test -- velocity
    err_stat = mean(abs(m1[,2] - model_at_analytical_x$y))/diff(range(m1[,2]))
    if(err_stat < ERR_TOL){
        print(c('PASS'))
    }else{
        print(c('FAIL', err_stat))
    }

    # Test 2 -- velocity
    model_at_analytical_x = approx(x$xs, x$ud[,yi, ind]/(x$stage[,yi,ind] -x$elev0[,yi]), xout=m1[,1])
    err_stat = mean(abs(m1[,3] - model_at_analytical_x$y))/diff(range(m1[,2]))
    if(err_stat < ERR_TOL){
        print(c('PASS'))
    }else{
        print(c('FAIL', err_stat))
    }
}

# Plot shoreline
shore = read.csv('DATA/Shoreline.csv')
shore_threshold = 0.15
shore_x = x$time * 0
shore_vel = shore_x
for(i in 1:length(shore_x)){
    ind = min(which(x$stage[,yi,i] - x$elev0[,yi] > shore_threshold))
    shore_x[i] = x$xs[ind]
    shore_vel[i] = x$ud[ind,yi,i]/(x$stage[ind,yi,i] - x$elev0[ind,yi])
}
plot(shore[,1], shore[,2], t='l', main='Shoreline \n Beware sensitivity to shore_threshold', xlab='Time', ylab='Position', lwd=3)
points(x$time, shore_x, t='o', col='red', pch=19, cex=0.2)
grid()
plot(shore[,1], shore[,3], t='l', main='Shoreline velocity \n See also Dellis et al 2008, who uses a 0.5m mesh', xlab='Time', ylab = 'Velocity (m/s)',
    ylim=c(-20, 15), lwd=3)
points(x$time, shore_vel, t='o', col='red', pch=19, cex=0.2)
grid()
dev.off()
