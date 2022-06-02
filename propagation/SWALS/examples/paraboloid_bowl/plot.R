source('../../plot.R')

x = get_all_recent_results(rev(Sys.glob('OUTPUTS/*/RUN_*'))[1])

ERR_TOL = 0.02

analytical_sol<-function(t, x, y){

    # run parameters
    gravity = 9.8
    r0 = 2000.0
    L = 2500.0
    D = 1000.0
    A = (L**4 - r0**4)/(r0**4 + L**4)
    omeg = 2.0 * sqrt(2.0 * gravity * D) / L

    r = sqrt(x^2 + y^2)
    z = D * (r^2/L^2 - 1)

    stage = D*((sqrt(1-A^2))/(1-A*cos(omeg*t)) -1 -r^2/L^2*(
        (1 - A^2)/((1 - A*cos(omeg*t))**2) - 1))
    vel_x = 0.5*omeg*r*A*sin(omeg*t) / (1-A*cos(omeg*t)) * sign(x)


    k = which(stage < z)
    stage[k] = z[k]
    vel_x[k] = 0

    output = list(stage=stage, vel_x = vel_x, t = t, elev=z)
}

#pdf('model_vs_data.pdf', width=10, height=12)

png('model_vs_data_stage_at_centre.png', width=10, height=12, units='in', res=300)

# Time-series in centre
analytical_centre = analytical_sol(x$time, 0, 0)
par(mfrow=c(2,1))
plot(analytical_centre$t, analytical_centre$stage, t='l', col=1,
    xlab='Time', ylab='Stage', main='Stage at the centre of the domain',
    cex.main=1.7)

xind = which.min(abs(x$xs))
yind = which.min(abs(x$ys))
points(x$time, x$stage[xind, yind,],t='p', col='red')

# Test the centre stage
model_at_analytical_t = approx(x$time, x$stage[xind, yind,], xout=analytical_centre$t)
err_stat = mean(abs(model_at_analytical_t$y - analytical_centre$stage))/diff(range(analytical_centre$stage))
if(err_stat < ERR_TOL){
    print('PASS')
}else{
    print('FAIL')
}

model_vel_centre = x$ud[xind,yind,]/(x$stage[xind,yind,]-x$elev[xind,yind,])
plot(x$time, model_vel_centre, col='red',
    xlab='Time', ylab='Velocity', 
    main='X-velocity at centre (should be zero)',
    ylim=c(-1,1)*0.01, cex.main=1.7)
points(analytical_centre$t, analytical_centre$vel_x, t='l', col=1)
legend('topright', c('Analytical', 'SWALS'), col=c('black', 'red'), lty=c(1, NA), pch=c(NA, 1),
       cex = 1.7, bty='n')

dev.off()

# Test the centre velocity
if(all(abs(model_vel_centre) < 1.0e-02)){
    print('PASS')
}else{
    print('FAIL')
}

# Check that the time ran long enough (i.e. did not blow up!)
if(max(x$time) < 120){
    print('FAIL')
}else{
    print('PASS')
}

# A few snapshots in time -- stage
png('model_vs_data_stage_over_time.png', width=10, height=12, units='in', res=300)
par(mfrow=c(4,3))
par(mar=c(3,3,3,1))
for(ts in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)){
    model_tind = which.min(abs(x$time - ts))
    analytical_centre = analytical_sol(x$time[model_tind], x$xs, 0)
    plot(x$xs, analytical_centre$stage, t='l', xlab='x', ylab='Stage',
        main=paste0('Stage, t=', round(x$time[model_tind], 3)), lwd=4, ylim=c(-500, 500),
        cex.main = 2, cex.axis=1.5)
    points(x$xs, x$stage[,yind, model_tind], col='red', t='o', cex=0.2, pch=19) 
    points(x$xs, x$elev[,yind, model_tind], col='blue', t='l', lty='dotted')
}
legend('bottom', c('Analytical', 'SWALS', 'Elevation'), lty=c('solid', 'solid', 'dotted'), 
       col=c('black', 'red', 'blue'), cex=2, bty='n')
dev.off()

# A few snapshots in time -- flux-x
png('model_vs_data_flux_x_over_time.png', width=10, height=12, units='in', res=300)
par(mfrow=c(4,3))
par(mar=c(3,3,3,1))
for(ts in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)){
    model_tind = which.min(abs(x$time - ts))
    analytical_centre = analytical_sol(x$time[model_tind], x$xs, 0)
    analytical_flux = analytical_centre$vel_x*(analytical_centre$stage - analytical_centre$elev)
    plot(x$xs, analytical_flux, 
         t='l', xlab='x', ylab='Vel_x',
        main=paste0('Flux_x, t=', round(x$time[model_tind], 3)), lwd=4,
        ylim=c(-1,1)*3e+04, cex.main = 2, cex.axis=1.5)
    points(x$xs, x$ud[,yind, model_tind], 
        col='red', t='o', cex=0.2, pch=19) 
}
legend('bottomright', c('Analytical', 'SWALS'), lty=c('solid', 'solid'), col=c('black', 'red'), cex=2, bty='n')
dev.off()

# A few snapshots in time -- vel-x
png('model_vs_data_vel_x_over_time.png', width=10, height=12, units='in', res=300)
par(mfrow=c(4,3))
par(mar=c(3,3,3,1))
for(ts in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120)){
    model_tind = which.min(abs(x$time - ts))
    analytical_centre = analytical_sol(x$time[model_tind], x$xs, 0)
    plot(x$xs, analytical_centre$vel_x, t='l', xlab='x', ylab='Vel_x',
        main=paste0('Vel_x, t=', round(x$time[model_tind],3)), lwd=4,
        ylim=c(-60,60), cex.main=2, cex.axis=1.5)
    points(x$xs, x$ud[,yind, model_tind]/(x$stage[,yind,model_tind] - x$elev[,yind,model_tind]), 
        col='red', t='o', cex=0.2, pch=19) 
}
legend('bottom', c('Analytical', 'SWALS'), lty=c('solid', 'solid'), col=c('black', 'red'), cex=2, bty='n')
dev.off()


# Animation showing the problematic 'slow drying'
# Note the most extreme wet points (awhile into the simulation) which don't seem to dry,
# looks problematic. In those regions, the change in elevation between cells is large
# compared with the flow depth (~10m vs <1m)
run_animation = FALSE
if(run_animation){
    for(tind in 1:300){
        par(mfrow=c(2,1))
        plot(x$stage[,201,tind] - x$elev[,201,tind], t='o', log='y', ylim=c(1.0e-03, 500))
        exact = analytical_sol(x$time[tind], x$xs, y=0)
        points(exact$stage - exact$elev, t='l', col='red')
        #points(x3$stage[,201,tind] - x3$elev[,201,tind], t='l', col='green')
        grid()
        plot(x$ud[,201,tind]/(x$stage[,201,tind] - x$elev[,201,tind]), t='o')
        exact = analytical_sol(x$time[tind], x$xs, y=0)
        points(exact$vel_x, t='l', col='red')
        #points(x3$ud[,201,tind]/(x3$stage[,201,tind] - x3$elev[,201,tind]), t='l', col='green')
        grid()
        Sys.sleep(0.2)
     }
}
