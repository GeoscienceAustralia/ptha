# Assuming the model has been run, we do
source('../../plot.R')
x = get_all_recent_results()

# Approximate the peak wave height with the max from the latter part of the
# simulation, when transients have died down
l = length(x$gauges$time)
ls = floor(2/3*l):l
peak_wave_height = apply(x$gauges$time_var$stage[,ls], 1, max)

pdf('model_data_comparison.pdf', width=10, height=10)
source('analytical_solution_zhang.R')
peak_wave_analytical_function = compute_zhang_solution()
peak_wave_analytical = peak_wave_analytical_function(x$gauges$lon, x$gauges$lat)

par(mfrow=c(1,1))
image(x$xs, x$ys, x$maxQ[,,1] * (x$maxQ[,,2] < 0), col=rainbow(255)[1:200], 
    asp=1, xlab='x', ylab='y', main='Peak stage around circular island, with gauge locations', 
    xlim=c(-1,1)*3e+05, ylim=c(-1,1)*3e+05)
points(x$gauges$lon, x$gauges$lat, pch=19)

# There are 3 sets of waves
par(mfrow=c(3,1))
l = length(peak_wave_analytical)/3
plot(peak_wave_height[1:l], t='o', main='Near Island (inner ring of points)', 
    ylab='Peak wave height', xlab='Gauge')
points(peak_wave_analytical[1:l], t='h', col='red')
grid()
legend('top', c('Model', 'Analytical'), col=c('black', 'red'), lwd=c(1,1), bty='n')

plot(peak_wave_height[1:l + l], t='o', main='Halfway along shelf (middle ring of points)', 
    ylab='Peak wave height', xlab='Gauge')
points(peak_wave_analytical[1:l + l], t='h', col='red')
grid()
plot(peak_wave_height[1:l + 2*l], t='o', main='Along shelf edge (outer ring of points)', 
    ylab='Peak wave height', xlab='Gauge')
points(peak_wave_analytical[1:l + 2*l], t='h', col='red')
grid()

dev.off()
