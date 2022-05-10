# Assuming the model has been run, we do
source('../../plot.R')

md_dir = sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1]
md = get_multidomain(md_dir, read_grids=FALSE, always_read_max_grids=TRUE)
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

## Get the 'exact' gauge location (i.e. mid-cell of model)
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

#
# Check mass conservation
#
md_log_files = Sys.glob(paste0(md_dir, '/*.log'))
log_data = readLines(md_log_files)
k = grep('unexplained change:', log_data)
last_k = k[length(k)]

unexplained_mass_change = as.numeric(strsplit(log_data[last_k], split=':')[[1]][2])
# The volume should be a few lines earier
md_volume_line = log_data[last_k - 3]
# Double check it's the correct line
if(!grepl('Multidomain volume', md_volume_line)){
    print('FAIL: (due to change in log file format?)')
}else{
    # Test mass conservation
    md_volume = as.numeric(strsplit(md_volume_line, split=':')[[1]][2])
     
    if(abs(unexplained_mass_change) < (1.0e-06 * md_volume)){
        print('PASS')
    }else{
        print('FAIL')
    }
}
