#
# Compare the test cases run with OPENMP with and without local time-stepping
#
# Local timestepping is not expected to lead to exact reproducibily. But here
# we check for similarity of the results (so that future updates leading to big changes will be detected)
#


multidomain_log_openmp = readLines(rev(Sys.glob('OUTPUTS/openmp_results/multidomain_log.log')))
multidomain_log_ompLocalTS = readLines(rev(Sys.glob('OUTPUTS/openmp_results_localtimestep/multidomain_log.log')))

# Get max/min stage, and max speed, globally over the domains
k = grep('Global stage range', multidomain_log_openmp)
max_stage_openmp = as.numeric(multidomain_log_openmp[k+1])
min_stage_openmp = as.numeric(multidomain_log_openmp[k+2])
k = grep('Global speed range', multidomain_log_openmp)
max_speed_openmp = as.numeric(multidomain_log_openmp[k+1])
k = grep('Time:', multidomain_log_openmp, fixed=TRUE)
times_openmp = sort(unique(as.numeric(multidomain_log_openmp[k+1])))

k = grep('Global stage range', multidomain_log_ompLocalTS)
max_stage_ompLocalTS = as.numeric(multidomain_log_ompLocalTS[k+1])
min_stage_ompLocalTS = as.numeric(multidomain_log_ompLocalTS[k+2])
k = grep('Global speed range', multidomain_log_ompLocalTS)
max_speed_ompLocalTS = as.numeric(multidomain_log_ompLocalTS[k+1])
k = grep('Time:', multidomain_log_ompLocalTS, fixed=TRUE)
times_ompLocalTS = sort(unique(as.numeric(multidomain_log_ompLocalTS[k+1])))

if(all( abs(times_ompLocalTS - times_openmp) < 1.0e-06)){
    print('PASS')
}else{
    print('FAIL')
}

#
# Compare max stage. 
#
k1 = which(times_ompLocalTS < 200) # Well behaved
k2 = which(times_ompLocalTS > 200 & times_ompLocalTS < 300) # Higher peak, which later reduces
k3 = which(times_ompLocalTS > 300 & times_ompLocalTS < 1500) # Quite good
k4 = which(times_ompLocalTS > 1500) # Offset issue occurs at the end
err_stat = abs(max_stage_openmp - max_stage_ompLocalTS)/diff(range(max_stage_ompLocalTS))
if(quantile(err_stat, 0.95) < 1.0e-05){
    print('PASS')
}else{
    print('FAIL')
}
# Compare min stage 
err_stat = abs(min_stage_openmp - min_stage_ompLocalTS)/diff(range(min_stage_ompLocalTS))
if(quantile(err_stat, 0.95) < 1.0e-03){
    print('PASS')
}else{
    print('FAIL')
}

# Compare max speed
err_stat = abs(max_speed_openmp - max_speed_ompLocalTS)/diff(range(max_speed_ompLocalTS))
if(quantile(err_stat, 0.95) < 1.0e-02 & quantile(err_stat, 0.5) < 1.0e-03){
    print('PASS')
}else{
    print('FAIL')
}

#
# Time-series of max stage and speed from the logfiles
#

png('Compare_openmp_ompLocalTS.png', width=7, height=9, units='in', res=300)
par(mfrow=c(3,1))
plot(times_openmp, max_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum stage (m)',
     main='Maximum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_ompLocalTS, max_stage_ompLocalTS, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'OpenmpLocalTS'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')

plot(times_openmp, min_stage_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Minimum stage (m)',
     main='Minimum wet stage over time', col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_ompLocalTS, min_stage_ompLocalTS, t='l', col='black', lty='dotted')
legend('bottomright', c('Openmp', 'OpenmpLocalTS'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')

plot(times_openmp, max_speed_openmp, t='l', lty='dashed', xlab='Time (s)', ylab='Maximum speed (m/s)',
     main='Maximum speed over time (note the linear domain may \n produce artificially high transient values, even if stable)', 
     col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_ompLocalTS, max_speed_ompLocalTS, t='l', col='black', lty='dotted')
legend('topright', c('Openmp', 'OpenmpLocalTS'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')
dev.off()


#
# Compare ompLocalTS and openmp at a particular time-slice
#

source('../../../plot.R')
suppressMessages(library(fields))

md_omp = 'OUTPUTS/openmp_results'
md_ca = 'OUTPUTS/openmp_results_localtimestep'

if(md_omp == md_ca){
    print('FAIL: Cannot find separate openmp/ompLocalTS directories')
}

# Plot difference between ompLocalTS/openmp at a particular time
time_inds = c(40, 400)
desired_var = 'uh' # 
domain_inds = 1:6
err_stats = rep(NA, length(domain_inds))
for(time_ind in time_inds){
    png(paste0('Compare_omp_ompLocalTS_time_index_', time_ind, '.png'), width=18, height=10, res=300, units='in')
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
        if(max(err_stats) < 1.0e-01 & min(err_stats) < 1e-05 & median(err_stats) < 0.02){
            print('PASS')
        }else{
            print('FAIL')
        }
    }else if(time_ind == time_inds[2]){
        # By now the models differ more significantly on some domains.
        target_errs = c(1e-03, 1e-03, 6e-03, 0.4, 1e+00, 0.02)
        if(all(err_stats < target_errs)){
            print('PASS')
        }else{
            print('FAIL')
        }

    }else{
        stop('NEED TO PROVIDE AN ERROR THRESHOLD FOR THIS CASE')
    }
}
