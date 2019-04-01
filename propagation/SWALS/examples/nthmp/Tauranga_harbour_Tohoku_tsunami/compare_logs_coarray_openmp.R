
multidomain_log_openmp = readLines(rev(Sys.glob('multidomain_log*openmp.log')))
multidomain_log_coarray = readLines(rev(Sys.glob('multidomain_log*coarray.log')))

# Get max/min stage, and max speed, globally over the domains
k = grep('Global stage range', multidomain_log_openmp)
max_stage_openmp = as.numeric(multidomain_log_openmp[k+1])
min_stage_openmp = as.numeric(multidomain_log_openmp[k+2])
max_speed_openmp = as.numeric(multidomain_log_openmp[k+4])
times_openmp = as.numeric(multidomain_log_openmp[k-15])

k = grep('Global stage range', multidomain_log_coarray)
max_stage_coarray = as.numeric(multidomain_log_coarray[k+1])
min_stage_coarray = as.numeric(multidomain_log_coarray[k+2])
max_speed_coarray = as.numeric(multidomain_log_coarray[k+4])
times_coarray = as.numeric(multidomain_log_coarray[k-15])

if(all( abs(times_coarray - times_openmp) < 1.0e-06)){
    print('PASS')
}else{
    print('FAIL')
}

# Compare max stage. Currently the test shows a differece at late times.
# In particular the cell with highest inundation is treated 'dry' in a 
# threshold manner, and slight roundoff significantly affects the result.
k = which(times_coarray < 125000)
err_stat = abs(max_stage_openmp - max_stage_coarray)/diff(range(max_stage_coarray))
if(all(err_stat[k] < 5.0e-4) & all(err_stat[-k] < 5e-02)){
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
#if(all(err_stat[k] < 1.0e-3) & all(err_stat[-k] < 5e-02)){
if(all(err_stat < 1.0e-3)){
    print('PASS')
}else{
    print('FAIL')
}

png('Compare_openmp_coarray.png', width=7, height=9, units='in', res=300)
par(mfrow=c(3,1))
plot(times_openmp, max_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum stage (m)',
     main='Maximum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_stage_coarray, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2)

plot(times_openmp, min_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Minimum stage (m)',
     main='Minimum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, min_stage_coarray, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2)

plot(times_openmp, max_speed_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum speed (m/s)',
     main='Maximum speed over time (note linear domains may \n produce artificially high transient values, even if stable)', 
     col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_speed_coarray, t='l', col='black', lty='dotted')
legend('topleft', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2)
dev.off()

