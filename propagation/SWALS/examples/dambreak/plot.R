#
# Simulation plot
#

# Dam-break problem parameters
H1 = as.numeric(commandArgs(trailingOnly=TRUE)[1]) #1
H0 = as.numeric(commandArgs(trailingOnly=TRUE)[2]) #0.0001
g = 9.8

near_t = 15 # Plot near this time

source('../../plot.R')
x = get_all_recent_results()

yind = round(dim(x$stage)[2]/2)

tind = which.min(abs(x$time - near_t))

# Get the analytical solution. Note it is reversed orientation, so we flip x and u
source('dam_break_analytical.R')
sol = wet_dam_break(t=x$time[tind], H0=H0, H1=H1)
sol0 = wet_dam_break(t=0, H0=H0, H1=H1)

energy_analytical = sol$energy
d_model = (x$stage[,yind,tind]-x$elev[,yind,tind])
u_model = x$ud[,yind,tind]/(d_model + 1.0e-10)
energy_numerical = (x$xs[2] - x$xs[1])*sum(u_model^2 * d_model + d_model^2*g) * 0.5

energy_numerical_0 = (x$xs[2] - x$xs[1])*sum(g*(x$stage[,yind,1]-x$elev[,yind,1])^2) * 0.5

output_file = paste0('dam_break_numerical_vs_analytical_H0_', H0, '_H1_', H1, '.png')
png(output_file, width=10, height=8, units='in', res=300)
par(mfrow=c(2,1))
plot(x$xs, x$stage[,yind,tind], t='o', ylim=c(0, max(sol$h)), 
     main=paste0('Stage, time=', x$time[tind], '\n Energy Numerical: ', 
                 round(energy_numerical, 2), '/', round(energy_numerical_0, 2), 
                 ', analytical: ', round(energy_analytical, 2), '/', round(sol0$energy, 2)), 
     cex=0.3, xlab='x (m)', ylab='Stage (m)')
points(-sol$x, sol$h, t='l', col='red')
legend('topleft', c('Model', 'Analytical'), lty=c(1,1), col=c('black', 'red'))

plot(x$xs, u_model, 
     t='o', ylim=c(min(-sol$u), 0), main='Velocity', cex=0.3, xlab='x (m)', ylab=' Velocity (m/s)')
points(-sol$x, -sol$u, t='l', col='red')
dev.off()

#
# Test
#
model_stage_fun = approxfun(x$xs, x$stage[,yind, tind], rule=2)
model_at_analytical_x = model_stage_fun(-sol$x)
test_stat = sum(abs(sol$h - model_at_analytical_x))/sum(abs(sol$h))
if(test_stat < 0.01){
    print('PASS')
}else{
    print('FAIL')
}
