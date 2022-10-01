sl = 0.97 # The ambient sea-level of the model run
source('../../../plot.R')

recent_output_dir = sort(Sys.glob('OUTPUTS/RUN_*'), decreasing=TRUE)[1]

# 2 multidomains
x = get_multidomain(recent_output_dir)

# Get gauges and labels for them. Note we do not have gauges WG1, WG2, and the
# one at the wavemaker, because our domain starts at x=5m (where it uses a boundary forcing)
gauges = merge_multidomain_gauges(x)
gauge_labind = floor(gauges$gaugeID/10)
gauge_label = paste0(c('wg', 'A', 'B', 'C', 'D')[gauge_labind+1], gauges$gaugeID - gauge_labind*10)

#
# Plot the boundary forcing, to ensure we got it right
#
png('boundary_check.png', width=8, height=5, units='in', res=300)
if(round(x[[1]]$lower_left_corner[1]) == 5){
    ts5 = read.table('problem_data/ts_5m.txt')
}else if(round(x[[1]]$lower_left_corner[1]) == 2){
    ts5 = read.table('problem_data/ts_from_observations_at_wg1.txt', header=FALSE)
}
wg_all = read.table('problem_data/Wavegage.txt', header=TRUE)

plot(x[[1]]$time, x[[1]]$stage[1,54,]-sl, t='l', col='purple',
    main='Gauge observations (@ x=2m), provided forcing (@ x=5m), \n and left-most-model-time-series (@ x=5m).',
    xlab='Time', ylab='Stage (m)', cex.axis=1.3, cex.lab=1.3)
points(ts5[,1:2], t='l', col='brown')
points(wg_all[,1], wg_all[, 5], t='l', col='orange')
points(wg_all[,1], wg_all[, 4], t='l', col='orange')
legend('topleft', c('Model @ x=5m', 'Provided forcing @ x=5m', 'Observed @ x=2m'), 
       col=c('purple', 'brown', 'orange'),
       lty=c(1,1,1), pch=c(NA, NA, NA), bty='n')
dev.off()

suppressMessages(library(fields))
png('Model_elevation_and_gauges.png', width=8.4, height=5, units='in', res=200)
image.plot(x[[1]]$xs, x[[1]]$ys, x[[1]]$elev0, zlim=c(0, 1.6), add=FALSE, 
           asp=1, main='Multidomain elevation and gauges',
           xlab='x (m)', ylab='y (m)', cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
image.plot(x[[2]]$xs, x[[2]]$ys, x[[2]]$elev0, zlim=c(0, 1.6), add=TRUE)
domain_bbox = get_domain_interior_bbox_in_multidomain(recent_output_dir)
polygon(domain_bbox$domain_interior_bbox[[1]], border='red')
polygon(domain_bbox$domain_interior_bbox[[2]], border='red')
all_gauges = read.csv('problem_data/point_gauge_locations.csv', header=TRUE)
points(all_gauges[,1:2], pch=19, col='black')
text(all_gauges[4:5, 1], all_gauges[4:5, 2], c('WG3', 'WG4'), pos=1)
dev.off()

png('Model_elevation_and_gauges_zoom.png', width=6, height=6.8, units='in', res=200)
image.plot(x[[2]]$xs, x[[2]]$ys, x[[2]]$elev0, zlim=c(0, 1.6), add=FALSE, 
           asp=1, main='Interior domain elevation and gauges',
           xlab='x (m)', ylab='y (m)', cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
domain_bbox = get_domain_interior_bbox_in_multidomain(recent_output_dir)
polygon(domain_bbox$domain_interior_bbox[[1]], border='red')
polygon(domain_bbox$domain_interior_bbox[[2]], border='red')
city_gauges = read.csv('problem_data/Gauge_locations_near_city.csv', header=TRUE)
points(city_gauges[,3:4], pch=19, col='black')
text(city_gauges[,3], city_gauges[,4], paste0(city_gauges[,2], city_gauges[,1]), 
     pos=1, #pos=rep(c(1,3), nrow(city_gauges)/2), 
     cex=0.8)
dev.off()

#
# Plot of WG3 and WG4. Note we do not have gauges WG1, WG2
# because our domain either starts at x=5m (where it uses a boundary forcing).
# The benchmark problem description suggested the ideal forcing should be OK, 
# but we see a clear phase-lag at these gauges (also reported by the HySea presentation).
# Probably better to use another forcing.
time_nudge = 0.7
data_wg = read.table('./problem_data/Wavegage.txt', header=TRUE)
k = which(gauge_label %in% c('wg3', 'wg4'))
png('gauges_wg3_wg4.png', width=8, height=5, units='in', res=300)
plot(data_wg$Time, data_wg$wg3, t='l', col='blue', xlim=c(0, 40),
   #main=paste0('Data @ WG3 & WG4, vs model. The phase-lag is noted in other \n',
   #            'papers and suggestings problems with the ideal boundary forcing. \n',
   #            'See Macais et al. (2020), doi: 10.1016/j.coastaleng.2020.103667'))
   main=paste0('Data vs model @ gauges WG3 & WG4. Model time shifted left by ', round(time_nudge, 2), ' sec \n',
               'to account for an apparent time-offset in the ideal boundary forcing, \n',
               'following Macais et al. (2020), doi: 10.1016/j.coastaleng.2020.103667'))

points(data_wg$Time, data_wg$wg4, t='l', col='blue')
points(gauges$time - time_nudge, gauges$time_var$stage[k[1],] - sl, t='l', col='orange')
points(gauges$time - time_nudge, gauges$time_var$stage[k[2],] - sl, t='l', col='orange')
legend('topleft', 
       c('Observed', 'Modelled with ideal forcing'), 
       lty=c(1,1), col=c('blue', 'orange'))
dev.off()

assert<-function(istrue){
    if(istrue){
        print('PASS')
    }else{
        print('FAIL')
    }
}
#
# Plot stage model-vs-data at gauges in the urban area.
# Note that the HySea presentation reports under-estimation at Gauge B1, and
# the TUFLOW paper [quad-tree + sub-grid] reports under-estimation at Gauge A1. 
#
gauge_files = paste0('./problem_data/other/Location_', gauge_label, '.txt')
for(gauge_group in c('A', 'B', 'C', 'D')){
    png(paste0('urban_gauge_group_', gauge_group, '_depth.png'), width=6, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    #par(mfrow=c(5,2))
    par(oma=c(1,1,0,0))
    par(mfrow=c(length(k), 1))
    par(mar=c(1, 2, 1, 1))
    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        plot(data[,1], data[,2], xlim=c(25, 40), ylim=c(0, 0.25), t='l')
        depth = gauges$time_var$stage[j,] - gauges$static_var$elevation0[j]
        points(gauges$time, depth, col='red', t='l')
        title(main=paste0(gauge_label[j], ': depth'), line=-2, cex.main=1.5)
        grid(col='orange')
        
        # Some tests
        err1 = ( max(depth*(gauges$time > 25), na.rm=TRUE)/max(data[,2] * (data[,1] > 25), na.rm=TRUE) - 1 )
        #print(c(gauge_label[j], round(err1, 2)))

        # Here I 'make up' some test criteria, as there are no objective criteria. This
        # will catch accidental big changes, but it may be very reasonable to change the criteria.
        if(gauge_label[j] == 'B1'){
            # At B1 (and C1, A1), presentations + papers from other groups
            # generally report waves smaller than observed (unless you move the
            # point slightly)
            assert(err1 > -0.35 & err1 < -0.1)
        }else if(gauge_label[j] == 'B4'){
            assert(abs(err1) < 0.25)
        }else if(gauge_label[j] == 'B6'){
            assert(abs(err1) < 0.25)
        }else if(gauge_label[j] == 'B9'){
            # This one is very small -- we underestimate, with similar results to Macias et al. (2020).
            assert(err1 > -0.7 & err1 < 0)
        }
    }
    legend('topright', c('Model', 'Data'), lty=c(1,1), col=c('red', 'black'), pch=c(NA, NA), bty='n', cex=1.5)
    dev.off()
}

# Plot of (speed^2 * depth) -- they call it 'flux', note it is the 'convective flux'
gauge_files = paste0('./problem_data/other/Location_', gauge_label, '.txt')
for(gauge_group in c('A', 'B', 'C', 'D')){
    png(paste0('urban_gauge_group_', gauge_group, '_convective_flux_hv2.png'), width=6, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    #par(mfrow=c(5,2))
    #par(mar=c(2, 2, 1, 1))
    par(oma=c(1,1,0,0))
    par(mfrow=c(length(k), 1))
    par(mar=c(1, 2, 1, 1))

    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        depth = pmax(gauges$time_var$stage[j,] - gauges$static_var$elevation0[j], 1.0e-04)
        plot(data[,1], data[,4], xlim=c(25, 40), ylim=c(0, 0.75), t='l')
        #points(gauges$time, (gauges$time_var$uh[j,]^2 + gauges$time_var$vh[j,]^2)/depth, 
        points(gauges$time, (gauges$time_var$uh[j,]^2)/depth,  # "Cross-shore"
               col='red', t='l')
        title(main=paste0(gauge_label[j], ': uuh'), line=-2, cex.main=1.5)
        grid(col='orange')
    }
    legend('topright', c('Model', 'Data'), lty=c(1,1), col=c('red', 'black'), pch=c(NA, NA), bty='n', cex=1.5)
    dev.off()

}

# Speed
gauge_files = paste0('./problem_data/other/Location_', gauge_label, '.txt')
for(gauge_group in c('A', 'B', 'C', 'D')){
    png(paste0('urban_gauge_group_', gauge_group, '_speed.png'), width=6, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    #par(mfrow=c(5,2))
    #par(mar=c(2, 2, 1, 1))
    par(oma=c(1,1,0,0))
    par(mfrow=c(length(k), 1))
    par(mar=c(1, 2, 1, 1))
    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        plot(data[,1], data[,3], xlim=c(25, 40), ylim=c(0, 3.0), t='l')
        depth = pmax(gauges$time_var$stage[j,] - gauges$static_var$elevation0[j], 1.0e-04)
        #speed = (sqrt(gauges$time_var$uh[j,]^2 + gauges$time_var$vh[j,]^2))/depth
        speed = (sqrt(gauges$time_var$uh[j,]^2))/depth # "Cross-shore"
        points(gauges$time, speed, col='red', t='l')
        title(main=paste0(gauge_label[j], ': x-velocity'), line=-2, cex.main=1.5)
        grid(col='orange')

        # Minimum time at which observations start
        minT = min(ifelse(!is.na(data[,3]), data[,1], 1e+06))

        # Some tests
        err1 = ( max(speed*(gauges$time > minT), na.rm=TRUE)/max(data[,3] * (data[,1] > minT), na.rm=TRUE) - 1 )
        #print(c(gauge_label[j], round(err1, 2)))

        # Here I 'make up' some test criteria, as there are no objective criteria. This
        # will catch accidental big changes, but it may be very reasonable to change the criteria.
        if(gauge_label[j] == 'B1'){
            assert(abs(err1) < 0.25)
        }else if(gauge_label[j] == 'B4'){
            assert(abs(err1) < 0.25)
        }else if(gauge_label[j] == 'B6'){
            assert(abs(err1) < 0.25)
        }else if(gauge_label[j] == 'B9'){
            assert(abs(err1) < 0.4)
        }
    }
    legend('topright', c('Model', 'Data'), lty=c(1,1), col=c('red', 'black'), pch=c(NA, NA), bty='n', cex=1.5)
    dev.off()
}
