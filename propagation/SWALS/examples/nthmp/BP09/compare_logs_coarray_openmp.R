
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

# Compare max stage. Currently the test shows a differece at late times.
# In particular the cell with highest inundation is treated 'dry' in a 
# threshold manner, and slight roundoff significantly affects the result.
#
# Also, like with other tests, we see transient times of higher difference, followed
# by 'going back to similarity'. This tests considers that
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
# min stage also benefits from a late-time treatment.
err_stat = abs(min_stage_openmp - min_stage_coarray)/diff(range(min_stage_coarray))
if(all(err_stat[k1] < 1.0e-8) & all(err_stat[k2] <  3e-04) & all(err_stat[k3] < 3e-01) & all(err_stat[k4] < 3e-01)){
    print('PASS')
}else{
    print('FAIL')
}

# max speed also benefits from a late-time treatment.
#k = which(times_coarray < 1000)
#k2 = which(times_coarray > 2600)
err_stat = abs(max_speed_openmp - max_speed_coarray)/diff(range(max_speed_coarray))
#if(all(err_stat[k] < 1.0e-3) & all(err_stat[-k] < 5e-02)){
if(all(err_stat[k1] < 1.0e-11) & all(err_stat[k2] <  2e-05) & all(err_stat[k3] < 2e-02) & all(err_stat[k4] < 3e-02)){
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
     main='Maximum speed over time (note the linear domain may \n produce artificially high transient values, even if stable)', 
     col='red', cex.main=2, cex.axis=1.5, cex.lab=1.6)
points(times_coarray, max_speed_coarray, t='l', col='black', lty='dotted')
legend('topright', c('Openmp', 'Coarray'), lty=c('dashed', 'dotted'), col=2:1, cex=2, bty='n')
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
# Note: Ultimately the models do differ significantly
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

    # Eventually the models diverge (apparently related to vortices, which
    # may be quite sensitive to even minor floating point differences)
    # So test the models in 2 ways
    if(time_ind == time_inds[1]){
        if(all(err_stats < 1.0e-4)){
            print('PASS')
        }else{
            print('FAIL')
        }
    }else if(time_ind == time_inds[2]){
        # By now the models differ more significantly on some domains.
        ## These were about 10x larger than the err_stats at the time of writing.
        ## So this is really a 'regression test'.
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
