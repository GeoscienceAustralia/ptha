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

multidomain_log_openmp = readLines(rev(Sys.glob('multidomain_log*openmp.log')))
multidomain_log_coarray = readLines(rev(Sys.glob('multidomain_log*coarray.log')))


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

# Compare max stage. 
err_stat = abs(max_stage_openmp - max_stage_coarray)/diff(range(max_stage_coarray))
if(all(err_stat < 1.0e-6)){
    print('PASS')
}else{
    print('FAIL')
}

# Compare min stage 
err_stat = abs(min_stage_openmp - min_stage_coarray)/diff(range(min_stage_coarray))
if(all(err_stat < 1.0e-6)){
    print('PASS')
}else{
    print('FAIL')
}

# Compare max speed
err_stat = abs(max_speed_openmp - max_speed_coarray)/diff(range(max_speed_coarray))
if(all(err_stat < 1.0e-6)){
    print('PASS')
}else{
    print('FAIL')
}

#
# Time-series plot of differences
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
     main='Maximum speed over time (note linear domains may \n produce artificially high transient values, even if stable)', 
     col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_speed_coarray, t='l', col='black', lty='dotted')
legend('topleft', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')
dev.off()

#
# Compare coarray and openmp at a particular time-slice
#

source('../../../plot.R')
suppressMessages(library(fields))

md_omp = Sys.glob('OUTPUTS/RUN*')[1]
md_ca = Sys.glob('OUTPUTS/RUN*')[2]

if(md_omp == md_ca){
    print('FAIL: Cannot find separate openmp/coarray directories')
}

# Plot difference between coarray/openmp at a particular time
time_ind =  700
desired_var = 'vh' # 'stage' does not show any differences, perhaps because the range is high (due to topography)?
domain_inds = 1
err_stats = rep(NA, length(domain_inds))
png(paste0('Compare_omp_coarray_time_index_', time_ind, '.png'), width=6, height=5, res=300, units='in')
par(mfrow = c(1,1))
for(domain_ind in domain_inds){
    xOMP = merge_domains_nc_grids(multidomain_dir=md_omp, domain_index=domain_ind, desired_var=desired_var, desired_time_index=time_ind)
    xCA  = merge_domains_nc_grids(multidomain_dir=md_ca,  domain_index=domain_ind, desired_var=desired_var, desired_time_index=time_ind)
    image.plot(xOMP$xs, xOMP$ys, xOMP[[desired_var]] - xCA[[desired_var]], asp=1)
    err_stats[domain_ind] = max(abs(xOMP[[desired_var]] - xCA[[desired_var]]), na.rm=TRUE)
}
dev.off()

if(all(err_stats < 1.0e-8)){
    print('PASS')
}else{
    print('FAIL')
}
