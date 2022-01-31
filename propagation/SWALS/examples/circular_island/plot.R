# Assuming the model has been run, we do
source('../../plot.R')

#x = get_all_recent_results(quiet=TRUE)

md = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])
x = md[[length(md)]] # Gives the domain of interest for BOTH single domain model and nested model

# Allowed mean error fraction for 'PASS'
REL_ERR = 0.05

# Approximate the peak wave height with the max from the latter part of the
# simulation, when transients have died down
l = length(x$gauges$time)
ls = floor(2/3*l):l
peak_wave_height = apply(x$gauges$time_var$stage[,ls], 1, max)


source('analytical_solution_zhang.R')
peak_wave_analytical_function = compute_zhang_solution()
peak_wave_analytical = peak_wave_analytical_function(x$gauges$lon, x$gauges$lat)

# Get the 'exact' gauge location (i.e. mid-cell of model)
#gx_ind = sapply(x$gauges$lon, f<-function(y) which.min(abs(y-x$xs)))
#gy_ind = sapply(x$gauges$lat, f<-function(y) which.min(abs(y-x$ys)))
#peak_wave_analytical = peak_wave_analytical_function(x$xs[gx_ind], x$ys[gy_ind])

#pdf('model_data_comparison.pdf', width=12, height=10)
png('model_data_comparison_1.png', width=12, height=10, res=300, units='in')

options(scipen=5)
suppressPackageStartupMessages(library(fields))
par(mfrow=c(1,1))
par(oma=c(1,1,1,1))
image.plot(x$xs, x$ys, x$maxQ * (x$elev0 < 0), col=rainbow(255)[1:200], 
    asp=1, main='Maximum-stage around conical island, with gauge locations', 
    xlim=c(-1,1)*3e+05, ylim=c(-1,1)*3e+05, xlab='x (m)', ylab='y(m)',
    cex.main=2, cex.lab=1.5, cex.axis=1.5, legend.lab='Max-stage (m)',
    legend.cex=2.0, axis.args=list(cex.axis=1.5))
points(x$gauges$lon, x$gauges$lat, pch=19)
dev.off()

png('model_data_comparison_2.png', width=12, height=10, res=300, units='in')
# There are 3 sets of waves
par(mfrow=c(3,1))
par(mar=c(5,5,4,1))
l = length(peak_wave_analytical)/3
plot(peak_wave_height[1:l], t='o', 
    ylab='Maximum-stage (m)', xlab='Gauge number (anti-clockwise from east)', 
    cex.lab=2, cex.axis=1.8)
title(main='Near Island (inner ring of points)', line=0.4, cex.main=3)
points(peak_wave_analytical[1:l], t='h', col='red', lwd=2)
grid()
legend('top', c('Model', 'Analytical'), col=c('black', 'red'), lwd=c(1,1), bty='n',
    cex=2)

# Simple test R
err_val = mean(abs(peak_wave_height[1:l] - peak_wave_analytical[1:l]))/mean(abs(peak_wave_analytical[1:l]))
if( err_val < REL_ERR){
    print('PASS')
}else{
    print(c('FAIL', err_val))
}

plot(peak_wave_height[1:l + l], t='o', 
    ylab='Maximum-stage (m)', xlab='Gauge number (anti-clockwise from east)', 
    cex.lab=2, cex.axis=1.8)
title(main='Halfway along shelf (middle ring of points)', line=0.4, cex.main=3)
points(peak_wave_analytical[1:l + l], t='h', col='red', lwd=2)
grid()

# Simple test
err_val = mean(abs(peak_wave_height[1:l + l] - peak_wave_analytical[1:l + l]))/mean(abs(peak_wave_analytical[1:l + l]))
if(err_val < REL_ERR){
    print('PASS')
}else{
    print(c('FAIL', err_val))
}

plot(peak_wave_height[1:l + 2*l], t='o', 
    ylab='Maximum-stage (m)', xlab='Gauge number (anti-clockwise from east)', 
    cex.lab=2, cex.axis=1.8)
title(main='Along shelf edge (outer ring of points)', line=0.4, cex.main=3)
points(peak_wave_analytical[1:l + 2*l], t='h', col='red', lwd=2)
grid()

# Simple test
err_val = mean(abs(peak_wave_height[1:l + 2*l] - peak_wave_analytical[1:l + 2*l])) / mean(abs(peak_wave_analytical[1:l + 2*l]))
if(err_val  < REL_ERR ){
    print('PASS')
}else{
    print(c('FAIL', err_val))
}

dev.off()

