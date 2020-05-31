timestepping_method = commandArgs(trailingOnly=TRUE)[1]

library(rptha)
source('../../plot.R')
md_dir = rev(Sys.glob('OUTPUTS/RUN*'))[1]
x = get_multidomain(md_dir)

# Logfile
md_file_lines = readLines(Sys.glob(paste0(md_dir, '/multidom*.log'))[1])

# Get dx/dy in cartesian-like coordinates
dx = distHaversine(cbind(x[[1]]$xs[1], x[[1]]$ys[1]), cbind(x[[1]]$xs[2], x[[1]]$ys[1]))
dy = distHaversine(cbind(x[[1]]$xs[1], x[[1]]$ys[1]), cbind(x[[1]]$xs[1], x[[1]]$ys[2]))

# Get coordinates x/y, in cartesian 
xmid = length(x[[1]]$xs)/2
ymid = length(x[[1]]$ys)/2
xc = ((1:length(x[[1]]$xs)) - xmid)*dx
yc = ((1:length(x[[1]]$ys)) - ymid)*dy

# Analytical solution -- see 
#    Sampson et al. (2006) Moving boundary shallow water flow above
#    parabolic bottom topography, ANZIAM J. 47 (EMAC2005) pp.C373–C387
a = 3000
h0 = 10
G = 10
grav = 9.8
f = 7.2726742371731926E-005
ohm = sqrt(f^2 + 2 * grav * h0 / a**2)

stage<-function(x, y, t){ 
    -G**2 * h0/a**2 * cos(ohm*t)**2 + 2*G*h0/a**2*cos(ohm*t)*x
}
uvel<-function(x, y, t) -G*ohm*sin(ohm*t)
vvel<-function(x, y, t) -G*f*cos(ohm*t)

# x-locations at which to check the solution.
# Beware the analytical stage can sometimes go below the bed elevation.
ks = c(which.min(abs(xc - 2500)), which.min(abs(xc+500)), which.min(abs(xc+1000)))

pdf(paste0('model_vs_analytical_', timestepping_method, '.pdf'), width=9, height=7)
# Stage time-series
par(mfrow=c(3,1))
for(k in ks){
    model_stage = x[[1]]$stage[k, floor(ymid),]
    plot(x[[1]]$time, model_stage, t='l', col='red', main=paste0('Stage timeseries @ x = ', round(xc[k], 2)))
    exact_stage = stage(xc[k], 0, x[[1]]$time)
    points(x[[1]]$time, exact_stage, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'))

    err = sum((model_stage - exact_stage)**2 / sum(exact_stage**2))
    if(err < 0.006) {
        print('PASS')
    }else{
        print(paste0('FAIL', err))
    }
}

# U timeseries
par(mfrow=c(3,1))
for(k in ks){
    depth = x[[1]]$stage[k,floor(ymid),] - x[[1]]$elev0[k,floor(ymid)]
    model_u = x[[1]]$ud[k, floor(ymid),]/depth
    plot(x[[1]]$time, model_u, t='l', col='red', 
         main=paste0('U-velocity timeseries @ x = ', round(xc[k], 2)))
    exact_u = uvel(xc[k], 0, x[[1]]$time)
    points(x[[1]]$time, exact_u, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'))

    err = sum((model_u - exact_u)**2 / sum(exact_u**2))
    if(err < 0.006) {
        print('PASS')
    }else{
        print(paste0('FAIL', err))
    }
}

# V timeseries
par(mfrow=c(3,1))
for(k in ks){
    depth = x[[1]]$stage[k,floor(ymid),] - x[[1]]$elev0[k,floor(ymid)]
    model_v =  x[[1]]$vd[k, floor(ymid),]/depth
    plot(x[[1]]$time, model_v, t='l', col='red', 
         main=paste0('V-velocity timeseries @ x = ', round(xc[k], 2)))
    exact_v = vvel(xc[k], 0, x[[1]]$time)
    points(x[[1]]$time, exact_v, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'))

    err = sum((model_v - exact_v)**2 / sum(exact_v**2))
    if(err < 0.006) {
        print('PASS')
    }else{
        print(paste0('FAIL', err))
    }
}
dev.off()

##
## Check some statistics from the log file
##

expected_initial_energy_total =  2.170550448602E+04
k = grep("Global energy-total", md_file_lines, fixed=TRUE)
energy_totals = as.numeric(md_file_lines[k+1])
# Check the starting energy is as expected
if(abs(energy_totals[1] - expected_initial_energy_total) < 0.001 * expected_initial_energy_total){
    print('PASS')
}else{
    print('FAIL')
}
n = length(energy_totals)
# Check the final energy is not too different -- noting some solvers are more dissipative than others
if(abs(energy_totals[n] - expected_initial_energy_total) < 0.05 * expected_initial_energy_total){
    print('PASS')
}else{
    print('FAIL')
}


expected_initial_volume_total = 1.997761893443E+07
k = grep("Multidomain volume    ", md_file_lines, fixed=TRUE)
volume_totals = sapply(md_file_lines[k], f<-function(x) as.numeric(strsplit(x, ':')[[1]][2]), USE.NAMES=FALSE)
# Check the starting volume is as expected
if(abs(volume_totals[1] - expected_initial_volume_total) < 0.001 * expected_initial_volume_total){
    print('PASS')
}else{
    print('FAIL')
}
n = length(volume_totals)
# Check the final volume is not too different -- noting some solvers are more dissipative than others
if(abs(volume_totals[n] - expected_initial_volume_total) < 0.001 * expected_initial_volume_total){
    print('PASS')
}else{
    print('FAIL')
}


