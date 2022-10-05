#
# Plot the test problem results
#

image_name_flag = commandArgs(trailingOnly=TRUE)
if(length(image_name_flag) != 1) stop('Should only pass one commandline argument')

# Tolerance for error at gauges in test
# The "motu" gauge has lots of high-frequency waves that are harder to model.
ERR_TOL_GAUGES = c(0.02, 0.026, 0.026, 0.07)


source('../../../plot.R')

#
# Get the model results
#
recent_dir = rev(Sys.glob('OUTPUTS/RUN_*'))[1]
model_gauges = merge_multidomain_gauges(multidomain_dir = recent_dir)

#
# Image plot of elevation, gauge locations, and domains
#
png('Model_elevation_and_gauge_locations.png', width=9, height=4.35, units='in', res=200)
par(mar=c(2,2,3,2))
# Elevation
elev_col = c(hcl.colors('Mako', n=40), hcl.colors('Lajolla', n=60)[1:40])
#elev_col = c(hcl.colors('Mako', n=40), hcl.colors('Green-Orange', n=40))
multidomain_image(recent_dir, variable='elevation0', time_index=NA, 
    var_transform_function=function(x) pmin(x, 40),  zlim=c(-40, 40), 
    xlim=c(0, 41000), ylim=c(0, 18600), col=elev_col, use_fields=TRUE)
title(main='Elevation, multidomain structure and gauge locations', cex.main=1.7)
# Domain boxes
bbox = get_domain_interior_bbox_in_multidomain(recent_dir)
for(i in 1:length(bbox$merged_domain_interior_bbox)){
    polygon(bbox$merged_domain_interior_bbox[[i]], border='red', lty='dashed')
}

# Gauge locations
pg = read.csv('point_gauge_locations.csv')
points(pg[,1], pg[,2], pch=4, col='black', cex=1.5)
text(pg[,1], pg[,2], c('A Beacon', 'Tug Berth', 'Sulphur Point', 'Moturiki', 'ADCP', ''), 
     col='black', pos=c(1, 4, 4, 4, 2, 1), cex=1.2)

dev.off()
#
# Extract gauge data. Interpretation of data based on matlab script provided in
# the test problem
#
obs = read.table('test_data/port_data.txt')
tobs = list()
tobs$abeacon = data.frame(time = obs[,1]*3600, stage=obs[,3])
tobs$tug = data.frame(time = obs[,4] * 3600, stage=obs[,6])
#tobs$taut = data.frame(time = obs[,1]*3600, stage=obs[,3]) ## Very near to TUG.
tobs$sulfur = data.frame(time = obs[,7]*3600, stage=obs[,9])
obs = read.table('test_data/tide_gauge.txt')
tobs$motu = data.frame(time = obs[,4]*3600, stage=obs[,6])

#
# Identify gauges associated with each entry in tobs
#
model_gauge_inds = match(c(1, 2, 3, 4), model_gauges$gaugeID)

#
# Make gauges plot
#
png(paste0('gauges_plot_', image_name_flag, '.png'), width=8, height=8, units='in', res=300)
par(mfrow=c(4, 1))
par(mar = c(2, 2, 2, 2))
par(oma = c(0, 2, 0, 0))
for(i in 1:length(tobs)){

    plot(tobs[[i]]$time/3600, tobs[[i]]$stage, t='l', xlim=c(0, 48),
        xlab='Hours', ylab='Stage (m)', cex.axis=1.5, cex.lab=1.7)
    title(paste0('Stage @ ', names(tobs)[i]), cex.main=1.8)
    grid()
    mgi = model_gauge_inds[i]
    points(model_gauges$time/3600, model_gauges$time_var$stage[mgi,], t='l', col='red')

    # SIMPLE TEST
    model_at_obs_times = approx(model_gauges$time/3600, model_gauges$time_var$stage[mgi,], xout=tobs[[i]]$time/3600)
    err_stat = mean(abs(tobs[[i]]$stage - model_at_obs_times$y), na.rm=TRUE)/diff(range(tobs[[i]]$stage, na.rm=TRUE))
    if(err_stat < ERR_TOL_GAUGES[i]){
        print(c('PASS', err_stat))
    }else{
        print(c('FAIL', err_stat))
    }
}
dev.off()

#
# Make ADCP plot
#

adcp = read.table('test_data/currents.txt')
obs_c = data.frame(time = adcp[,1]*3600, speed=sqrt(rowSums(adcp[,6:7]^2)))

#
ind = match(5, model_gauges$gaugeID)
model_time = model_gauges$time
model_speed = sqrt(model_gauges$time_var$uh[ind,]^2 + model_gauges$time_var$vh[ind,]^2)/
    (model_gauges$time_var$stage[ind,] - model_gauges$static_var$elevation0[ind])

# This seems closer to the 'right' location'.
ind = match(6, model_gauges$gaugeID)
model_speed_extra = sqrt(model_gauges$time_var$uh[ind,]^2 + model_gauges$time_var$vh[ind,]^2)/
    (model_gauges$time_var$stage[ind,] - model_gauges$static_var$elevation0[ind])

png(paste0('currents_', image_name_flag, '.png'), width=10, height=6, units='in', res=200)
plot(obs_c$time/3600, obs_c$speed, t='l', xlim=c(0, 48), xlab='Time (hours)', ylab='Speed (m/s)')
points(model_time/3600, model_speed, t='l', col='red')
points(model_time/3600, model_speed_extra, t='l', col='blue', lty='dotted')
title('Currents at ADCP (Improves with higher resolution). \n Underestimates reported by several codes, adcp position uncertainty?')
legend('topleft', c('Data', 'Model @ reported coord', 'Model near reported coord'),
    lty=c('solid','solid','dotted'), col=c('black', 'red', 'blue'))
# In the NTHMP report, several groups discuss issues with the locations:
# See https://nws.weather.gov/nthmp/documents/NTHMP_Currents_Workshop_Report.pdf
#
# For CLIFFS solver, Tolkova notes that (page 85),
#    "The simulated current at the prescribed ADCP position was lower than  the
#    measurements  (Figure  5).  However,  the  current  rapidly  varied
#    across  narrow  harbor entrance where the ADCP was positioned, being fast
#    on the channel centerline, and low next to the banks. The flow transition
#    from low by the shore to high in the  channel occupied a slightly wider
#    zone on the numerical grid. Indeed, simulated current 127 m (3 nodes
#    diagonally) from the prescribed position toward the channel centerline
#    provided nearly perfect fit to the observations (Figure 5, red)"
#
# For the MOST solver, Arcas notes that (p 71) 
#   "Computed values of wave elevation were also in extremely good agreement with observations in
#    the case of tsunami+tide, however, current speed results for this second case tend to underestimate
#    observed values for the tidal component of the signal as seen in Figure 5. This deficiency is most
#    probably due to the lack of knowledge about velocity initial and boundary conditions along grid
#    boundaries when tidal effects are included. A quick analysis of the accuracy of the results suggests
#    that accurate calculation of the velocity solution, particularly when large tidal currents are present
#    will require some a priori knowledge of tidal velocity boundary conditions."
#
# In the presentations from that workshop, the simulations from Kirby (FUNWAVE) have a similar issue.
dev.off()

