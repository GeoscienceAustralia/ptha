timestepping_method = commandArgs(trailingOnly=TRUE)[1]

source('../../plot.R')
md_dir = rev(Sys.glob('OUTPUTS/RUN*'))[1]
x = get_multidomain(md_dir)

# Logfile
md_file_lines = readLines(Sys.glob(paste0(md_dir, '/multidom*.log'))[1])

# Get dx/dy in cartesian-like coordinates
dx = as.numeric(strsplit(md_file_lines[grep('_bottom_edge(1)', md_file_lines, fixed=TRUE)], '   ')[[1]][2])
dy = as.numeric(strsplit(md_file_lines[grep('_left_edge(1)'  , md_file_lines, fixed=TRUE)], '   ')[[1]][2])

# Get coordinates x/y, in cartesian 
xmid = length(x[[1]]$xs)/2 + 0.5
ymid = length(x[[1]]$ys)/2 + 0.5
xc = ((1:length(x[[1]]$xs)) - xmid)*dx
yc = ((1:length(x[[1]]$ys)) - ymid)*dy

# Analytical solution -- see 
#    Sampson et al. (2006) Moving boundary shallow water flow above
#    parabolic bottom topography, ANZIAM J. 47 (EMAC2005) pp.C373–C387
a = 3000
h0 = 10
G = 200 # 10
grav = 9.8
f = 7.2925440502248539E-005
ohm = sqrt(f^2 + 2 * grav * h0 / a**2)

stage<-function(x, y, t){ 
    -G**2 * h0/a**2 * cos(ohm*t)**2 + 2*G*h0/a**2*cos(ohm*t)*x
}
uvel<-function(x, y, t) -G*ohm*sin(ohm*t)
vvel<-function(x, y, t) -G*f*cos(ohm*t)
elev<-function(x, y){ - h0 * (1.0 - (x/a)**2) }

# x-locations at which to check the solution.
# Beware the analytical stage can sometimes go below the bed elevation.
ks = rev( c(which.min(abs(xc - 3000)), which.min(abs(xc+500)), which.min(abs(xc+1000))) )

# Stage time-series
png(paste0('model_vs_analytical_stage_', timestepping_method, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(2.7,3,2,1))
for(k in ks){
    model_stage = x[[1]]$stage[k, floor(ymid),]
    plot(x[[1]]$time, model_stage, t='o', col='red', 
        main=paste0('Stage timeseries @ x = ', round(xc[k], 2)),
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
    exact_stage = pmax( stage(xc[k], 0, x[[1]]$time), elev(xc[k], 0) )
    points(x[[1]]$time, exact_stage, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'),  
           bg=rgb(1,1,1,alpha=0.7), cex=1.5)

    err = sum((model_stage - exact_stage)**2 / sum(exact_stage**2))
    if(err < 0.006) {
        print('PASS')
    }else{
        print(paste0('FAIL STG', err))
    }
}
dev.off()

# U timeseries
png(paste0('model_vs_analytical_U_', timestepping_method, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(2.7,3,2,1))
for(k in ks){
    depth = x[[1]]$stage[k,floor(ymid),] - x[[1]]$elev0[k,floor(ymid)]
    model_u = x[[1]]$ud[k, floor(ymid),]/depth #* (depth > 0.001)
    exact_u = uvel(xc[k], 0, x[[1]]$time) * ( stage(xc[k], 0, x[[1]]$time) > elev(xc[k], 0))

    is_wetdry_zone = any(stage(xc[k], 0, x[[1]]$time) < elev(xc[k], 0))

    YLIM = range(c(model_u, exact_u), na.rm=TRUE)

    plot(x[[1]]$time, model_u, t='o', col='red', 
        main=paste0('U-velocity timeseries @ x = ', round(xc[k], 2)),
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=YLIM)
    points(x[[1]]$time, exact_u, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'),
           bg=rgb(1,1,1,alpha=0.7), cex=1.5)

    if(!is_wetdry_zone){
        # Avoid threshold checks on speed in the wet/dry zone, as it can be sensitive to speeds
        # when the flow is nearly dry. Instead we check the flux (below)
        err = sum((model_u - exact_u)**2 / sum(exact_u**2))
        if(err < 0.01) {
            print('PASS')
        }else{
            print(paste0('FAIL U', err))
        }
    }
}
dev.off()


# UH timeseries
png(paste0('model_vs_analytical_UH_', timestepping_method, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(2.7,3,2,1))
for(k in ks){
    #depth = x[[1]]$stage[k,floor(ymid),] - x[[1]]$elev0[k,floor(ymid)]
    model_u = x[[1]]$ud[k, floor(ymid),]# /depth #* (depth > 0.001)
    exact_u = uvel(xc[k], 0, x[[1]]$time) * pmax( stage(xc[k], 0, x[[1]]$time) - elev(xc[k], 0), 0)

    YLIM = range(c(model_u, exact_u), na.rm=TRUE)

    plot(x[[1]]$time, model_u, t='o', col='red', 
        main=paste0('UH-flux timeseries @ x = ', round(xc[k], 2)),
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=YLIM)
    points(x[[1]]$time, exact_u, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'),
           bg=rgb(1,1,1,alpha=0.7), cex=1.5)

    err = sum((model_u - exact_u)**2 / sum(exact_u**2))
    if(err < 0.01) {
        print('PASS')
    }else{
        print(paste0('FAIL UH', err))
    }
}
dev.off()

# V timeseries
png(paste0('model_vs_analytical_V_', timestepping_method, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(2.7,3,2,1))
for(k in ks){
    depth = x[[1]]$stage[k,floor(ymid),] - x[[1]]$elev0[k,floor(ymid)]
    model_v =  x[[1]]$vd[k, floor(ymid),]/depth #* (depth > 0.001)
    exact_v = vvel(xc[k], 0, x[[1]]$time) * ( stage(xc[k], 0, x[[1]]$time) > elev(xc[k], 0))

    is_wetdry_zone = any(stage(xc[k], 0, x[[1]]$time) < elev(xc[k], 0))

    YLIM = range(c(model_v, exact_v), na.rm=TRUE)
    plot(x[[1]]$time, model_v, t='o', col='red', 
        main=paste0('V-velocity timeseries @ x = ', round(xc[k], 2)),
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=YLIM)
    points(x[[1]]$time, exact_v, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'),
           bg=rgb(1,1,1,alpha=0.7), cex=1.5)

    if(!is_wetdry_zone){
        err = sum((model_v - exact_v)**2 / sum(exact_v**2))
        if(err < 0.01) {
            print('PASS')
        }else{
            print(paste0('FAIL V', err))
        }
    }
}
dev.off()

# VH timeseries
png(paste0('model_vs_analytical_VH_', timestepping_method, '.png'), width=9, height=7, units='in', res=300)
par(mfrow=c(3,1))
par(mar=c(2.7,3,2,1))
for(k in ks){
    model_v =  x[[1]]$vd[k, floor(ymid),]
    exact_v = vvel(xc[k], 0, x[[1]]$time) * pmax( stage(xc[k], 0, x[[1]]$time) - elev(xc[k], 0), 0 )
    YLIM = range(c(model_v, exact_v), na.rm=TRUE)
    plot(x[[1]]$time, model_v, t='o', col='red', 
        main=paste0('VH-flux timeseries @ x = ', round(xc[k], 2)),
        cex.main=1.5, cex.axis=1.5, cex.lab=1.5, ylim=YLIM)
    points(x[[1]]$time, exact_v, t='l', col='black')
    legend('topright', c('Analytical', 'Model'), lty=c(1, 1), pch=c(NA,NA), col=c('black', 'red'),
           bg=rgb(1,1,1,alpha=0.7), cex=1.5)

    err = sum((model_v - exact_v)**2 / sum(exact_v**2))
    if(err < 0.01) {
        print('PASS')
    }else{
        print(paste0('FAIL VH', err))
    }
}
dev.off()

##
## Check some statistics from the log file
##

expected_initial_energy_total = 8.703570853927E+006 # 2.170550448602E+04
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
# Check the final volume is not too different
if(abs(volume_totals[n] - expected_initial_volume_total) < 0.001 * expected_initial_volume_total){
    print('PASS')
}else{
    print('FAIL')
}


