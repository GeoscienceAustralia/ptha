#
# Try to create a stage/velocity time-series which will 
# lead to observations at gauges 1-4 being naturally reproduced
#

obs_files = c('../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2a.txt',
              '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2b.txt',
              '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/ts2cnew1.txt')

wavemaker_file = '../test_repository/BP06-FrankG-Solitary_wave_on_a_conical_island/fdbk2abc.txt'
wavemaker_forcing = read.table(wavemaker_file, skip=1)
wavemaker_tstart = 22 # Wavemaker time at model-time = 0

gauge_x_distance_from_wavemaker = list(5.76, 6.82, 7.56)
depth_MSL = 0.32
gauge_tstart = 22 # Gauge time at model-time 0
g = 9.8


gauges_data = lapply(obs_files, f<-function(x) read.table(x, header=FALSE, skip=8))


#
# Take the wave at gauge2, and estimate the wave that would have been at the wavemaker,
# accounting for depth-nonlinearity in the propagation speed, but assuming the wave
# propagates with a linear stage-vs-velocity relation.
#
# This is heuristic. But it's an interesting alternative to the wavemaker boundary
# (which has its own challenges).
#
# Case C currently leads to shock formation. The treatment below is hacky. 
# But given the rough assumptions underlying the approach, it may not be worthwhile to
# try to improve.
#
gauge_time = lapply(gauges_data, f<-function(x) x[,1] - gauge_tstart)
## Get gauge 2 (in column 3), while subtracting the mean of the first 3s of gauge
#gauge_2 = lapply(gauges_data, f<-function(x) x[,3] - mean(x[ which(x[,1] < gauge_tstart + 3.0), 3]))
# Get gauge 2 (in column 3
gauge_2 = lapply(gauges_data, f<-function(x) x[,3])
gauge_time_adjusted = vector(mode='list', length=3)
gauge_fun = vector(mode='list', length=3)
wavemaker_fun = gauge_fun
for(i in 1:length(gauge_time_adjusted)){
    
    # Estimate wave velocity assuming a plane wave
    vel =  sqrt(g/(depth_MSL + gauge_2[[i]]))*(gauge_2[[i]]) 
    
    #
    # wave-speed = gravity wave speed + local velocity 
    # (velocity is in the same direction for this example during the leading
    # part of the wave)
    wave_speed = sqrt(g * (depth_MSL + gauge_2[[i]])) + vel * (vel > 0)
    travel_time = gauge_x_distance_from_wavemaker[[i]] / (wave_speed)
    gauge_time_adjusted[[i]] = gauge_time[[i]] - travel_time
    # Check for time reversals
    k = which(diff(gauge_time_adjusted[[i]]) < 0)
    # Should correct this using conservation-law arguments
    # The following case-specific treatment might fail in general
    if(length(k) > 0){
        #browser()
        gauge_time_adjusted[[i]] = gauge_time_adjusted[[i]][-( min(k-2):max(k+2) )]
        gauge_2[[i]] = gauge_2[[i]][-(min(k-2):max(k+2))]
    }
    gauge_fun[[i]] = approxfun(gauge_time_adjusted[[i]], gauge_2[[i]])
    wavemaker_fun[[i]] = approxfun(wavemaker_forcing[,1] - wavemaker_tstart, wavemaker_forcing[,i+1])
}


tseq = seq(0.0, 10, len=1001)
wavemaker_implied_vel = lapply(wavemaker_fun, f<-function(wf) (wf(tseq+0.1) - wf(pmax(tseq-0.1, 0)))/(0.2) * 1/100) 
gauge_implied_vel = lapply(gauge_fun, f<-function(gf) sqrt(g/depth_MSL) * gf(tseq))

par(mfrow=c(3,1))
for(i in 1:3){
    plot(tseq, wavemaker_implied_vel[[i]], t='l')
    points(tseq, gauge_implied_vel[[i]], t='l', col='red')
}

for(i in 1:3){
    out_file = paste0('gauge_forcing_nonlinear_case', i, '.csv')
    out_data = cbind(tseq, gauge_implied_vel[[i]])
    write.table(out_data, file=out_file, sep=',', row.names=FALSE, col.names=c('time (s)', 'x-vel (m)'))
}
