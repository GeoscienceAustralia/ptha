source('../../../plot.R')
library(raster)

# Tolerance for agreement in model/obs time-series at a few gauges.
# This restricts the " mean(abs(error)) / observed_range "
GAUGE_TS_ERR_TOL = 0.1

# Note this will read the 'final' domain, which is the high res one.
x = get_all_recent_results(sort(Sys.glob('OUTPUTS/RUN_*/RUN_*'), decreasing=TRUE)[1], quiet=TRUE)
m1 = read.csv('../test_repository/BP07-DmitryN-Monai_valley_beach/output_ch5-7-9.csv')

#pdf('gauges_plot.pdf', width=18, height=12)
png('gauges_plot.png', width=12, height=8, units='in', res=200)
par(mfrow=c(3,1))
par(mar=c(5.1,5.1, 4.1, 2.1))
plot(x$gauges$time, x$gauges$time_var$stage[1,], t='o', ylim=c(-0.01, 0.045), xlim=c(0,30),
    xlab='Time (s)', ylab='Stage (m)', cex.lab=2, cex.axis=2)
points(m1[,1], m1[,2]/100, t='l', col='red', lwd=2)
grid(); title('Gauge 5', cex.main=2.5)
legend('topleft', c('Model', 'Observed'), lty=rep('solid',2), pch=c(1, NA), col=c('black', 'red'), 
       bty='n', cex=2)
plot(x$gauges$time, x$gauges$time_var$stage[2,], t='o', ylim=c(-0.01, 0.045), xlim=c(0,30),
    xlab='Time (s)', ylab='Stage (m)', cex.lab=2, cex.axis=2)
points(m1[,1], m1[,3]/100, t='l', col='red', lwd=2)
grid(); title('Gauge 7', cex.main=2.5)
plot(x$gauges$time, x$gauges$time_var$stage[3,], t='o', ylim=c(-0.01, 0.045), xlim=c(0,30),
    xlab='Time (s)', ylab='Stage (m)', cex.lab=2, cex.axis=2)
points(m1[,1], m1[,4]/100, t='l', col='red', lwd=2)
grid(); title('Gauge 9', cex.main=2.5)
dev.off()

#
# Pass/Fail tests that we reasonably match the observed time-series at the
# gauges plotted above.
#
# Restrict the test to the model time, and AFTER wave arrives (noting the
# initial stages in the experiment are irregular)
k = which(m1[,1] < 33 & m1[,1] > 8)
# Gauges 5, 7, 9 are in columns 2:4 of m1, in units of cm
for(i in 2:4){

    model_at_gauge_times = approx(x$gauges$time, x$gauges$time_var$stage[i-1,], xout=m1[k,1])
    err_stat = mean(abs(model_at_gauge_times$y*100 - m1[k, i])/diff(range(m1[k,i])))
    #print(err_stat)

    if(err_stat < GAUGE_TS_ERR_TOL){
        print('PASS')
    }else{
        print('FAIL')
    }
}

im_xlim = c(4.6, 5.2)
im_ylim = c(1.5, 2.2)

# The ts is flexible (photo times are not exact, but spaced by 0.5s), note the GEOCLAW validation
ts = c(15.3, 15.8, 16.3, 16.8, 17.3) - 0.2
rasts = paste0('../test_repository/BP07-DmitryN-Monai_valley_beach/Frame_', c(10, 25, 40, 55, 70),'_line.png')
for(i in 1:5){


    png(paste0('snapshot_time_', round(ts[i], 3), '.png'), width=18, height=12, units='in', res=200)
    par(mfrow=c(1,2))
    image(x$xs, x$ys, x$elev[,,1], col=grey(seq(0,1,len=255)), xlim=im_xlim, ylim=im_ylim, asp=1,
        cex.lab=2, cex.axis=2, xlab='x(m)', ylab='y(m)')
    contour(x$xs, x$ys, x$elev[,,1], add=TRUE, levels=seq(-.1, 0.1, by=0.01))

    ind = which.min(abs(x$time - ts[i]))
    wet = x$stage[,,ind] > (x$elev[,,1] + 2.0e-03)
    
    contour(x$xs, x$ys, wet, add=TRUE, levels=c(0.5), col='red', lwd=3) 
    title(paste0('Time ', ts[i], ' s'), cex.main=2.5)

    points(x$gauges$lon, x$gauges$lat, col='green', pch=19)
    text(x$gauges$lon, x$gauges$lat, col='green', c('G5', 'G7', 'G9'), pos=3, cex=1.5)

    snapshot1 = brick(rasts[i])
    plotRGB(flip(t(snapshot1), direction='x'))

    dev.off()
}

#dev.off()


# Peak runup at points, from 6 runs of the experiment [included in NTHMP test problem]
## 
##                                             Case NO.
##                      109_105  109_106  109_107  210_101  210_102  210_103  
##    x(m)      y(m)                          runup(m)
## 
##   5.1575   1.8800    0.0875   0.09     0.08     0.09     0.1      0.09
##   5.0300   2.2062    0.0525   0.0575   0.055    0.065    0.0675   0.065
##   4.9975   2.3200    0.0525   0.055    0.055    0.0575   0.0575   0.0575

# Indices at each point
inds_list = list(
    inds1 = c(which.min(abs(x$xs - 5.1575)), which.min(abs(x$ys - 1.8800))),
    inds2 = c(which.min(abs(x$xs - 5.0300)), which.min(abs(x$ys - 2.2062))),
    inds3 = c(which.min(abs(x$xs - 4.9975)), which.min(abs(x$ys - 2.3200)))
)
# Observed range in repetitions of the experiment
range_list = list(
    site1 = c(0.08, 0.1),
    site2 = c(0.0525, 0.0675),
    site3 = c(0.0525, 0.0575)
)

#print('Checking peak stage at 3 reported gauges')
modelled_stage_maxima = rep(NA, 3)
experimental_maxima = rep(NA, 3)
experimental_minima = rep(NA, 3)
for(i in 1:3){

    # Check that the gauge level is between the observed levels in repeated experiments
    stagei = x$maxQ[inds_list[[i]][1], inds_list[[i]][2]]
    rangei = range_list[[i]]
    if(stagei < rangei[2] & stagei > rangei[1]){
        #print(c('PASS', stagei))
        print('PASS')
    }else{
        #print(c('FAIL', stagei))
        print('FAIL')
    }

    modelled_stage_maxima[i] = stagei
    experimental_maxima[i] = rangei[2]
    experimental_minima[i] = rangei[1]
}

# Store the runup values
outputs = data.frame(
    x = c(5.1575, 5.0300, 4.9975),
    y = c(1.8800, 2.2062, 2.3200),
    model = modelled_stage_maxima,
    experiment_lower_limit = experimental_minima,
    experiment_upper_limit = experimental_maxima)

write.csv(round(outputs, 4), 'model_vs_experiment_test_result.csv', row.names=FALSE)
