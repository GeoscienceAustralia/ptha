library(ncdf4)
# Work on most recent output file
nc_file = sort(Sys.glob('OUTPUTS/RUN*/RUN*/Grid*.nc'), decreasing=TRUE)[1]

fid = nc_open(nc_file)
x = ncvar_get(fid, 'x')
elev = ncvar_get(fid, 'elevation0')[,3] # Only get middle entry
stage = ncvar_get(fid, 'stage')[,3,] # Only get middle entry
time = ncvar_get(fid, 'time')

# The code assumes particular parameter values (otherwise we have to change the
# gauge locations, see Grilli et al 1994)
H0 = 0.44
h = 0.2*H0
initial_h = max(stage[,1]*(elev<0))
min_elev = min(elev)
if((abs(initial_h - h) >1e-04) | (abs(min_elev + H0) > 1e-05)){
    print('FAIL: Test code assumes H/d = 0.2 and d=0.44')
}else{
    print('PASS')
}

# Location of gauges g0, g1, ... for this case
gauge_name = paste0('g', 0:9)
# Gauge locations from Table 1 except g0 (-5), and flipped around because my
# model domain has the forcing on the right, and my slope does not begin at x=0
# gauge_x = c(-5, 20.96, 21.41, 22.55, 23.23, 23.68, 24.14, 24.68, 25.34, 26.91)*H0
## NOTE: The gauge-9 value looks like a typo in the original paper and seems to be
## treated that way by Filippina et al. (2015)
gauge_x = c(-5, 20.96, 21.41, 22.55, 23.23, 23.68, 24.14, 24.68, 25.34, 25.91)*H0
gauge_ind = sapply(gauge_x, function(y) which.min(abs(x-y)))

g0 = gauge_ind[which(gauge_name == 'g0')]
#plot(time, stage[g0,], t='l')
model_peak_time_g0 = time[which.max(stage[g0,])]
#abline(v=13.56*sqrt(H0/9.8), col='red')

tNorm = sqrt(H0/9.8) # Used to normalise the time coordinate in the paper

# Get digitised observations at g0
obs_g0 = read.csv('g0_fig.csv')
#obs_peak_time_g0 = 13.56 * tNorm # Peak time used in Grilli et al.
obs_peak_time_g0 = 13.75*tNorm # Better alignment
phase_offset = -obs_peak_time_g0 + model_peak_time_g0

png('Solitary_shoaling_Grilli94.png', width=8, height=10, units='in', res=200)
par(mfrow=c(2,1))
plot((time - phase_offset)/tNorm, stage[g0, ]/H0, t='l', xlim=c(5, 45), 
    main='Observations and model at gauge 0 (aligned with manual phase offset)',
    xlab='Time/sqrt(h0/g)', ylab='Stage/h0', cex.axis=1.4, cex.lab=1.4)
points(obs_g0[,1], obs_g0[,2], col='red')
grid(col='orange')

# Get digitized observations at other gauges
obs_g1_g9 = read.csv('g1_to_g9.csv')
gauge_obs = c('g1', 'g3', 'g5', 'g7', 'g9')
# The digiting software begins all curves at the same time, but this isn't right,
# so remove the initial times depending on the curve
gauge_obs_start_t = c(37, 37.5, 39.35, 41.333, 41.7)
gauge_obs_t = obs_g1_g9[,1]

for(i in 1:length(gauge_obs)){

    ii = gauge_ind[which(gauge_obs[i] == gauge_name)]

    if(gauge_obs[i] == 'g1'){
        plot((time - phase_offset)/tNorm, stage[ii,]/H0, t='l', xlim=c(37, 47), ylim=c(-0.1, 0.4), col=i,
            xlab='Time/sqrt(h0/g)', ylab='Stage/h0', cex.axis=1.4, cex.lab=1.4)
        points(gauge_obs_t, obs_g1_g9[,i+1]*(gauge_obs_t > gauge_obs_start_t[i]), col=i, pch=i)
    }else{
        points((time - phase_offset)/tNorm, stage[ii,]/H0, t='l', xlim=c(37, 47), ylim=c(-0.1, 0.4), col=i)
        points(gauge_obs_t, obs_g1_g9[,i+1]*(gauge_obs_t > gauge_obs_start_t[i]), col=i, pch=i)
    }

    obs_peak = max(obs_g1_g9[,i+1])
    model_peak = max(stage[ii,])/H0
    reltol = 0.03
    if(abs(model_peak - obs_peak) < reltol*obs_peak){
        print(c('PASS'))
    }else{
        print(c('FAIL', obs_peak, model_peak, (obs_peak-model_peak)/obs_peak))
    }
}
legend('bottom', gauge_obs, col=1:length(gauge_obs), lty=1, pch=i, horiz=TRUE, bty='n', cex=1.2)
title(main='Observations (points) vs model (lines) at gauges 1,3,5,7,9')
grid(col='orange')
dev.off()

# Plot the initial condition
png('Initial_condition.png', width=9, height=5, units='in', res=200)
plot(x, elev, t='l', ylim=c(-0.5, 0.5), xlab='x/h0 ', ylab='Stage/h0 ', cex.axis=1.4, cex.lab=1.4,
    main='Initial condition', cex.main=1.7)
points(x, stage[,1], t='l', col='red')
grid(col='orange')
ii = c(1,2,10)
text(gauge_x[ii], 0.4, gauge_name[ii], cex=1.3)
abline(v=gauge_x, lty='dotted')
legend('bottomright', c('Elevation', 'Stage', 'Gauges'), lty=c('solid','solid', 'dotted'), col=c('black', 'red', 'black'))
dev.off()
