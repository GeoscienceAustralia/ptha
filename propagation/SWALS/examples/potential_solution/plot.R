suppressMessages(library(raster))
source('../../plot.R')

source('potential_solution.R')

md_dir = rev(sort(Sys.glob('OUTPUTS/RUN*')))[1]
x = get_multidomain(md_dir, read_grids=FALSE, read_gauges=FALSE)
ti = length(x[[1]]$time)

# Useful to have the timestepping method
fid = nc_open(Sys.glob(paste0(x[[1]]$output_folder, '/Grid*.nc')))
ts_method = ncatt_get(fid, 0, 'timestepping_method')$value
nc_close(fid)

# Create the solution at the final time
if(ts_method %in% c('leapfrog_nonlinear', 'linear', 'leapfrog_linear_plus_nonlinear_friction')){
    # For the staggered grid models we should subtract the timestep/2 to get the time-staggering correct.
    solend = make_final_solution(x[[1]]$time[ti] - 1.2*401/201*0.5)
}else{
    solend = make_final_solution(x[[1]]$time[ti])
}

STAGE = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=1, desired_var='stage', desired_time_index=ti)

png('numerical_vs_potential_solution.png', width=7, height=6, units='in', res=200)
plot(xs, solend[,floor(length(ys)/2)], t='l', xlab='Distance (m)', ylab='Stage (m) @ time 100',
    ylim=c(-1.4, 1.0), lwd=5, lty='dashed', main=paste0('Solution at final time = ', round(x[[1]]$time[ti],3) ))
# Numerical solution along y==0
points(STAGE$xs, STAGE$stage[,floor(length(STAGE$ys)/2)], t='l', col='red')
# Numerical solution along x==0 -- note we have an uneven grid size
points(STAGE$ys, STAGE$stage[floor(length(STAGE$xs)/2), ], t='l', col='orange')
dx = mean(diff(STAGE$xs))
dy = mean(diff(STAGE$ys))
legend('bottomleft', 
       c('Potential wave theory', 
         paste0('SWALS with dispersion, rectangular cells (y==0, dx=', signif(dx,2), ')'), 
         paste0('SWALS with dispersion, rectangular cells (x==0, dy=', signif(dy,2), ')')), 
       lty=c('dashed','solid','solid'), lwd=c(5,1,1), pch=c(NA, NA, NA), 
       col=c('black', 'red', 'orange'), bg=NA, bty='n', cex=1.2)
dev.off()
