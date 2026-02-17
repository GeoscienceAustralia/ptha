suppressMessages(library(raster))
source('../../../plot.R')
source('potential_solution.R')

plot_tag = commandArgs(trailingOnly=TRUE)[1] # Include in the plot filenames

# Get the simulations
md_dir = rev(sort(Sys.glob('OUTPUTS/RUN*')))[1]
x = get_multidomain(md_dir, read_grids=FALSE, read_gauges=FALSE)
ti = length(x[[1]]$time)
di = get_domain_indices_in_multidomain(md_dir)
has_two_domains = (2 %in% di)

# Get the timestepping method, because for staggered schemes we will consider
# the time staggering.
timestepping_methods = rep("", length(x))
for(i in 1:length(x)){
    fid = nc_open(Sys.glob(paste0(x[[i]]$output_folder, '/Grid*.nc')))
    timestepping_methods[i] = ncatt_get(fid, 0, 'timestepping_method')$value
    nc_close(fid)
}

# Read the final time stage
STAGE = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=1, desired_var='stage', desired_time_index=ti)
if(has_two_domains) STAGE2 = merge_domains_nc_grids(multidomain_dir=md_dir, domain_index=2, desired_var='stage', desired_time_index=ti)

# Create the solution at the final time
solend = flip_like_raster(make_final_solution(x[[1]]$time[ti]))
staggered_schemes = c('leapfrog_nonlinear', 'linear', 'leapfrog_linear_plus_nonlinear_friction')
using_staggered_schemes = any(timestepping_methods %in% staggered_schemes)
if(using_staggered_schemes){
    # For the staggered grid models there is some ambiguity in the 'final'
    # timestep (since stage is offset by dt/2 compared to UH/VH). To reflect this 
    # we can also overlap the solution timestep/2 behind the recorded final time, using
    # the outer-grid timestep (which anyway influences the inner grid too)
    solend_lag = flip_like_raster( make_final_solution(x[[1]]$time[ti] - 1.2*401/dim(STAGE$stage)[1]*0.5) )
}else{
    solend_lag = solend
}

#
# Plot transects of the numerical and analytical solution at the final time
#

png(paste0('numerical_vs_potential_solution_', plot_tag, '.png'), width=7, height=6, units='in', res=200)
potential_y = ys
potential_stage = solend[floor(length(xs)/2),]
plot(potential_y, potential_stage, t='l', xlab='Distance (m)', ylab='Stage (m)',
    ylim=c(-1, 1.0)*2, lwd=5, lty='solid', main=paste0('Solution at final time = ', round(x[[1]]$time[ti],3) ))
if(using_staggered_schemes) points(ys, solend_lag[floor(length(xs)/2),], t='l', col='green', lwd=2)

points_and_pass_fail_check<-function(model_x, model_stage, col){
    # Convenience function to add a line to the plot, and also do a PASS/FAIL check on the results

    points(model_x, model_stage, t='l', col=col)

    # PASS/FAIL check
    reference_result = approx(potential_y, potential_stage, xout=model_x)$y
    k = which(!(is.na(reference_result) | is.na(model_stage)))
    err_stat = mean(abs(reference_result - model_stage)[k])
    err_tol = 0.015
    if(err_stat < err_tol){
        print('PASS')
    }else{
        print(c('FAIL', err_stat))
    }

}

# Numerical solution along x==0 -- note we have an uneven grid size
points_and_pass_fail_check(STAGE$ys, STAGE$stage[floor(length(STAGE$xs)/2), ], col='red')
if(has_two_domains){
    # Nested grid
    k = which.min(abs(STAGE2$xs)) 
    points_and_pass_fail_check(STAGE2$ys, STAGE2$stage[k, ], col='red')
}

dx = mean(diff(STAGE$xs))
dy = mean(diff(STAGE$ys))
if(using_staggered_schemes){
    legend('bottomleft', 
           c('Potential wave theory', 
             'Potential wave theory (1/2 dt offset staggered schemes)',
             paste0('SWALS with dispersion, (x==0, coarse dy=', signif(dy,2), ')')), 
           lty=c('solid','solid', 'solid'), lwd=c(3,2,1), pch=c(NA, NA, NA), 
           col=c('black','green', 'red'), bg=NA, bty='n', cex=1.2)
}else{
    legend('bottomleft', 
           c('Potential wave theory', 
             paste0('SWALS with dispersion, (x==0, coarse dy=', signif(dy,2), ')')), 
           lty=c('solid','solid'), lwd=c(3,1), pch=c(NA, NA), 
           col=c('black', 'red'), bg=NA, bty='n', cex=1.2)
}
title(sub=paste0(unique(timestepping_methods), collapse=" & "))
grid()
dev.off()

#
# Plot the final water surface, and the analytical solution
#
tmp = get_domain_interior_bbox_in_multidomain(md_dir)
if(has_two_domains){
    interior_box = tmp$merged_domain_interior_bbox[[2]]
    interior_box = rbind(interior_box, interior_box[1,])
}
png(paste0('numerical_free_surface_', plot_tag, '.png'), width=14, height=7, units='in', res=200)
par('mfrow'=c(1,2))
col_pal = colorRampPalette(c('darkblue', 'blue', 'lightblue', 'green', 'lightgreen', 'yellow', 'beige', 'orange', 'red', 'purple'))(200)
colz = col_pal #c(rainbow(125)[1:100], hcl.colors('Blue-Red 3', n=100))
col_breaks = 2.02 * asin(seq(-1,1,len=length(colz)+1))/(pi/2)
image(STAGE$xs, STAGE$ys, STAGE$stage, zlim=range(col_breaks), col=colz, breaks=col_breaks, xlab='x', ylab='y', main='Numerical')
if(has_two_domains){
    image(STAGE2$xs, STAGE2$ys, STAGE2$stage, zlim=range(col_breaks), col=colz, breaks=col_breaks, add=TRUE)
}
title(sub=paste0(unique(timestepping_methods), collapse=" & "))
if(has_two_domains) points(interior_box[,1], interior_box[,2], lty='dotted', t='l', col='grey')

image(xs, ys, solend, zlim=range(col_breaks), col=colz, breaks=col_breaks, main='Analytical', xlab="x", ylab='y')
dev.off()

