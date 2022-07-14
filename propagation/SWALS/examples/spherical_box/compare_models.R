# Compare the 3 model types

source('../../plot.R')

md_runs = sort(Sys.glob('OUTPUTS/RUN*'))

# We should have 3 runs -- leapfrog_nonlinear, cliffs, rk2
stopifnot(length(md_runs) == 3)

# Read the multidomains
all_md = lapply(md_runs, get_multidomain)

#
# Check the timestepping methods are as expected, as we use
# this logic below
# 
all_nc_files = lapply(all_md, function(x) Sys.glob(paste0(x[[1]]$output_folder, '/Grid*.nc')))
all_timestepping_methods = unlist(lapply(all_nc_files, function(x){
    fid = nc_open(x, readunlim=FALSE)
    ts_method = ncatt_get(fid, var=0)$timestepping_method
    nc_close(fid)
    return(ts_method)
}) )

expected_ts_methods = c('leapfrog_nonlinear', 'cliffs', 'rk2')
if(all(all_timestepping_methods == expected_ts_methods)){
    print('PASS')
}else{
    print('FAIL')
}

# Get the stage at the last timestep for comparison
ind = length(all_md[[1]][[1]]$time)
all_last_stage = lapply(all_md, function(x) x[[1]]$stage[,,ind])
names(all_last_stage) = expected_ts_methods

xs = all_md[[1]][[1]]$xs
ys = all_md[[1]][[1]]$ys

# Check that differences in the final-time stage are small
# There are some isolated larger differences in shallower waters where length-scales
# become smaller, but in general the differences are not important, and images of 
# the three models look pretty much the same
ERR_TOL = 0.01
for(model_ind in 2:3){
    if(mean(abs(all_last_stage[[1]] - all_last_stage[[model_ind]]) ) < ERR_TOL){
        print('PASS')
    }else{
        print('FAIL')
    }
}

library(fields)
png('initial_conditions.png', width=10, height=4, units='in', res=200)
par(mfrow=c(1,2))
par(oma = c(0, 0, 0, 2))
image.plot(xs, ys, all_md[[1]][[1]]$stage[,,1], zlim=c(-4,4), main='Initial stage (all)', cex.main=1.5, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
image.plot(xs, ys, all_md[[1]][[1]]$elev[,,1], main='Initial elevation (all)', cex.main=1.5, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
dev.off()

png('three_models_at_final_time.png', width=12, height=4, units='in', res=200)
par(mfrow=c(1,3))
image.plot(xs, ys, all_last_stage[[1]], zlim=c(-4,4), main='leapfrog_nonlinear \n stage @ 1h', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
image.plot(xs, ys, all_last_stage[[2]], zlim=c(-4,4), main='cliffs \n stage @ 1h', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
image.plot(xs, ys, all_last_stage[[3]], zlim=c(-4,4), main='rk2 \n stage @ 1h', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
dev.off()


png('model_differences_at_final_time.png', width=12, height=4, units='in', res=200)
par(mfrow=c(1,3))
image.plot(xs, ys, all_last_stage[[1]] - all_last_stage[[2]], zlim=c(-4,4), main='stage difference @1h \n (leapfrog_nonlinear - cliffs)', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
image.plot(xs, ys, all_last_stage[[1]] - all_last_stage[[3]], zlim=c(-4,4), main='stage difference @1h \n (leapfrog_nonlinear - rk2)', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
image.plot(xs, ys, all_last_stage[[2]] - all_last_stage[[3]], zlim=c(-4,4), main='stage difference @1h \n (cliffs - rk2)', cex.main=2, xlab='Lon', ylab='Lat', cex.axis=1.3, cex.lab=1.3)
dev.off()
