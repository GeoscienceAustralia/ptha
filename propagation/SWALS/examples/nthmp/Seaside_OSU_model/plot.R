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
plot(x[[1]]$time, x[[1]]$stage[1,54,]-sl, t='l', col='purple',
    main='Overplot of boundary-forcing, and left-most-model-time-series. \n Should agree well.')
points(ts5[,1:2], t='l', col='brown')
dev.off()

#
# Plot of WG3 and WG4. Note we do not have gauges WG1, WG2
# because our domain either starts at x=5m (where it uses a boundary forcing).
# The benchmark problem description suggested the ideal forcing should be OK, 
# but we see a clear phase-lag at these gauges (also reported by the HySea presentation).
# Probably better to use another forcing.
data_wg = read.table('./problem_data/Wavegage.txt', header=TRUE)
k = which(gauge_label %in% c('wg3', 'wg4'))
png('gauges_wg3_wg4.png', width=8, height=5, units='in', res=300)
plot(data_wg$Time, data_wg$wg3, t='l', col='blue', xlim=c(0, 40),
   main=paste0('Data @ WG3 & WG4, vs model. The phase-lag is noted in other \n',
               'papers and suggestings problems with the ideal boundary forcing. \n',
               'See Macais et al. (2020), doi: 10.1016/j.coastaleng.2020.103667'))
points(data_wg$Time, data_wg$wg4, t='l', col='blue')
points(gauges$time, gauges$time_var$stage[k[1],] - sl, t='l', col='orange')
points(gauges$time, gauges$time_var$stage[k[2],] - sl, t='l', col='orange')
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
    png(paste0('urban_gauge_group_', gauge_group, '_stage.png'), width=9, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    par(mfrow=c(5,2))
    par(mar=c(2, 2, 1, 1))
    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        plot(data[,1], data[,2], xlim=c(25, 40), ylim=c(0, 0.25), t='o', pch=1)
        depth = gauges$time_var$stage[j,] - gauges$static_var$elevation0[j]
        points(gauges$time, depth, col='red', t='o', pch=1)
        title(main=gauge_label[j], line=-2, cex.main=1.5)
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
    dev.off()
}

# Plot of (speed^2 * depth) -- they call it 'flux', note it is the 'convective flux'
gauge_files = paste0('./problem_data/other/Location_', gauge_label, '.txt')
for(gauge_group in c('A', 'B', 'C', 'D')){
    png(paste0('urban_gauge_group_', gauge_group, '_convective_flux_hv2.png'), width=9, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    par(mfrow=c(5,2))
    par(mar=c(2, 2, 1, 1))
    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        depth = pmax(gauges$time_var$stage[j,] - gauges$static_var$elevation0[j], 1.0e-04)
        plot(data[,1], data[,4], xlim=c(25, 40), ylim=c(0, 0.75), t='o', pch=1)
        points(gauges$time, (gauges$time_var$uh[j,]^2 + gauges$time_var$vh[j,]^2)/depth, 
               col='red', t='o', pch=1)
        title(main=gauge_label[j], line=-2, cex.main=1.5)
        grid(col='orange')
    }
    dev.off()
}

# Speed
gauge_files = paste0('./problem_data/other/Location_', gauge_label, '.txt')
for(gauge_group in c('A', 'B', 'C', 'D')){
    png(paste0('urban_gauge_group_', gauge_group, '_speed.png'), width=9, 
        height=9, units='in', res=300)
    k = which(startsWith(gauge_label, gauge_group))
    par(mfrow=c(5,2))
    par(mar=c(2, 2, 1, 1))
    for(j in k){
        data = read.table(gauge_files[j], comment="%")        
        plot(data[,1], data[,3], xlim=c(25, 40), ylim=c(0, 3.0), t='o', pch=1)
        depth = pmax(gauges$time_var$stage[j,] - gauges$static_var$elevation0[j], 1.0e-04)
        speed = (sqrt(gauges$time_var$uh[j,]^2 + gauges$time_var$vh[j,]^2))/depth
        points(gauges$time, speed, col='red', t='o', pch=1)
        title(main=gauge_label[j], line=-2, cex.main=1.5)
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
    dev.off()
}
