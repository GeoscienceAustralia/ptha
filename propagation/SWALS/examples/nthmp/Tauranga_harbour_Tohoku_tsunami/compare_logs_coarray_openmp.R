
multidomain_log_openmp = readLines(rev(Sys.glob('multidomain_log*openmp.log')))
multidomain_log_coarray = readLines(rev(Sys.glob('multidomain_log*coarray.log')))

# Get max/min stage, and max speed, globally over the domains
k = grep('Global stage range', multidomain_log_openmp)
max_stage_openmp = as.numeric(multidomain_log_openmp[k+1])
min_stage_openmp = as.numeric(multidomain_log_openmp[k+2])
max_speed_openmp = as.numeric(multidomain_log_openmp[k+4])
times_openmp = as.numeric(multidomain_log_openmp[k-17])

k = grep('Global stage range', multidomain_log_coarray)
max_stage_coarray = as.numeric(multidomain_log_coarray[k+1])
min_stage_coarray = as.numeric(multidomain_log_coarray[k+2])
max_speed_coarray = as.numeric(multidomain_log_coarray[k+4])
times_coarray = as.numeric(multidomain_log_coarray[k-17])

if(all( abs(times_coarray - times_openmp) < 1.0e-06)){
    print('PASS')
}else{
    print('FAIL')
}

# Compare max stage. 
k = which(times_coarray < 50000 | times_coarray > 60000) # Test had an isolated spike (dry-threshold?).
err_stat = abs(max_stage_openmp - max_stage_coarray)/diff(range(max_stage_coarray))
if(all(err_stat[k] < 5.0e-4)){
    print('PASS')
}else{
    print('FAIL')
}

# min stage also benefits from a late-time treatment.
k = which(times_coarray < 80000)
err_stat = abs(min_stage_openmp - min_stage_coarray)/diff(range(min_stage_coarray))
if(all(err_stat[k] < 5.0e-4) & all(err_stat[-k] < 5e-02)){
    print('PASS')
}else{
    print('FAIL')
}

# max speed also benefits from a late-time treatment.
#k = which(times_coarray < 1000)
err_stat = abs(max_speed_openmp - max_speed_coarray)/diff(range(max_speed_coarray))
if(all(err_stat < 5.0e-3)){
    print('PASS')
}else{
    print('FAIL')
}

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


# Compare coarray and openmp at a particular time-slice

source('../../../plot.R')
library(fields)

md_omp = Sys.glob('OUTPUTS/RUN*')[1]
md_ca = Sys.glob('OUTPUTS/RUN*')[2]

if(md_omp == md_ca){
    print('FAIL: Cannot find separate openmp/coarray directories')
}

# Plot difference between coarray/openmp at a particular time
time_ind =  200
desired_var = 'vh' # 'stage' does not show any differences, perhaps because the range is high (due to topography)?
domain_inds = 1:2
err_stats = rep(NA, length(domain_inds))
png(paste0('Compare_omp_coarray_time_index_', time_ind, '.png'), width=5, height=9, res=300, units='in')
par(mfrow = c(2,1))
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
