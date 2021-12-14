#
# Compare the test cases run with OPENMP / COARRAYS
#
# Ideally we'd like exactly the same simulation result in either case.
# This is nontrivial though, because of floating point round-off and
# non-associativity. 
#
# By default, if we don't pass a load-balance-partition.txt file to the model,
# then SWALS will partition the domains in the COARRAY case, but not
# the OPENMP case. This leads to tiny differences (e.g. 1e-10) in the 
# initial condition (I assume because of round-off scale differences in
# coordinates and other domain properties). Those small difference will slowly
# grow over time, and eventually be noticable (albeit not important).
#
# The same 'slowly growing difference' occurs if we compare 2 OPENMP models,
# one having initial conditions perturbed by 1e-10. Given that, it is not surprising
# that differences emerge when the domain is partitioned differently.
#
# A similar thing seems to happen in the ROMS model: https://www.myroms.org/forum/viewtopic.php?f=30&t=3741
#
# However, if we desire "exactly the same" answer with OPENMP and COARRAY, we can
# give the code a 'load-balance-partition.txt' file which specifies the same partition
# in both cases. Then, I seem to be getting 'near perfect' reproducibility in
# the flow variables at the end of the simulation. Interestingly the 'mass
# conservation round-off' still shows some differences -- perhaps because of the
# double-precision CO_SUM operations to get the full domain volume in the COARRAY case?
#


multidomain_log_openmp = readLines(rev(Sys.glob('OUTPUTS/openmp_results/multidomain_log.log')))
multidomain_log_coarray = readLines(rev(Sys.glob('OUTPUTS/coarray_results/multidomain_*.log')[1]))

# Get max/min stage, and max speed, globally over the domains
k = grep('Global stage range', multidomain_log_openmp)
max_stage_openmp = as.numeric(multidomain_log_openmp[k+1])
min_stage_openmp = as.numeric(multidomain_log_openmp[k+2])
k = grep('Global speed range', multidomain_log_openmp)
max_speed_openmp = as.numeric(multidomain_log_openmp[k+1])
k = grep('Time:', multidomain_log_openmp, fixed=TRUE)
times_openmp = sort(unique(as.numeric(multidomain_log_openmp[k+1])))

k = grep('Global stage range', multidomain_log_coarray)
max_stage_coarray = as.numeric(multidomain_log_coarray[k+1])
min_stage_coarray = as.numeric(multidomain_log_coarray[k+2])
k = grep('Global speed range', multidomain_log_coarray)
max_speed_coarray = as.numeric(multidomain_log_coarray[k+1])
k = grep('Time:', multidomain_log_coarray, fixed=TRUE)
times_coarray = sort(unique(as.numeric(multidomain_log_coarray[k+1])))

if(all( abs(times_coarray - times_openmp) < 1.0e-06)){
    print('PASS')
}else{
    print('FAIL')
}

#
# Compare max stage. These tolerances were made with an earlier version of the code
#
k1 = which(times_coarray < 200) # Well behaved
k2 = which(times_coarray > 200 & times_coarray < 300) # Higher peak, which later reduces
k3 = which(times_coarray > 300 & times_coarray < 1500) # Quite good
k4 = which(times_coarray > 1500) # Offset issue occurs at the end
err_stat = abs(max_stage_openmp - max_stage_coarray)/diff(range(max_stage_coarray))
if(all(err_stat[k1] < 1.0e-8) & all(err_stat[k2] <  3e-05) & all(err_stat[k3] < 1e-06) & all(err_stat[k4] < 5e-02)){
    print('PASS')
}else{
    print('FAIL')
}
# Compare min stage 
err_stat = abs(min_stage_openmp - min_stage_coarray)/diff(range(min_stage_coarray))
if(all(err_stat[k1] < 1.0e-8) & all(err_stat[k2] <  3e-04) & all(err_stat[k3] < 3e-01) & all(err_stat[k4] < 3e-01)){
    print('PASS')
}else{
    print('FAIL')
}

# Compare max speed
err_stat = abs(max_speed_openmp - max_speed_coarray)/diff(range(max_speed_coarray))
#if(all(err_stat[k] < 1.0e-3) & all(err_stat[-k] < 5e-02)){
if(all(err_stat[k1] < 1.0e-11) & all(err_stat[k2] <  2e-05) & all(err_stat[k3] < 2e-02) & all(err_stat[k4] < 3e-02)){
    print('PASS')
}else{
    print('FAIL')
}

#
# Time-series of max stage and speed from the logfiles
#

png('Compare_openmp_coarray.png', width=7, height=9, units='in', res=300)
par(mfrow=c(3,1))
plot(times_openmp, max_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum stage (m)',
     main='Maximum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_stage_coarray, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')

plot(times_openmp, min_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Minimum stage (m)',
     main='Minimum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, min_stage_coarray, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')

plot(times_openmp, max_speed_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum speed (m/s)',
     main='Maximum speed over time (note the linear domain may \n produce artificially high transient values, even if stable)', 
     col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_speed_coarray, t='l', col='black', lty='dotted')
legend('topright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')
dev.off()


#
# Compare coarray and openmp at a particular time-slice
#

source('../../../plot.R')
suppressMessages(library(fields))

md_omp = 'OUTPUTS/openmp_results'
md_ca = 'OUTPUTS/coarray_results'

if(md_omp == md_ca){
    print('FAIL: Cannot find separate openmp/coarray directories')
}

# Plot difference between coarray/openmp at a particular time
time_inds = c(40, 400)
desired_var = 'uh' # 
domain_inds = 1:6
err_stats = rep(NA, length(domain_inds))
for(time_ind in time_inds){
    png(paste0('Compare_omp_coarray_time_index_', time_ind, '.png'), width=18, height=10, res=300, units='in')
    par(mfrow = c(3,4))
    par(oma=c(0, 1, 0, 1))
    for(domain_ind in domain_inds){
        xOMP = merge_domains_nc_grids(multidomain_dir=md_omp, domain_index=domain_ind, desired_var=desired_var, desired_time_index=time_ind)
        xCA  = merge_domains_nc_grids(multidomain_dir=md_ca,  domain_index=domain_ind, desired_var=desired_var, desired_time_index=time_ind)
        image.plot(xOMP$xs, xOMP$ys, xOMP[[desired_var]] - xCA[[desired_var]], asp=1, main=paste0(domain_ind, ' difference'))
        err_stats[domain_ind] = max(abs(xOMP[[desired_var]] - xCA[[desired_var]]), na.rm=TRUE)/diff(range(xOMP[[desired_var]], na.rm=TRUE))
        image.plot(xOMP$xs, xOMP$ys, xOMP[[desired_var]], asp=1, main=paste0(domain_ind, ' omp'))
    }
    dev.off()

    # Use different 'allowed errors' at short and long times
    if(time_ind == time_inds[1]){
        # Small shorter-time difference
        if(all(err_stats < 1.0e-4)){
            print('PASS')
        }else{
            print('FAIL')
        }
    }else if(time_ind == time_inds[2]){
        # By now the models differ more significantly on some domains.
        # These were about 10x larger than the err_stats at the time of writing.
        target_errs = c(6e-04, 8e-04, 6e-03, 5e-01, 1e+00, 4e-04)
        if(all(err_stats < target_errs)){
            print('PASS')
        }else{
            print('FAIL')
        }

    }else{
        stop('NEED TO PROVIDE AN ERROR THRESHOLD FOR THIS CASE')
    }
}
