source('../../../plot.R')
# Since we call SWALS from inside R, let's pass the openmp commands.
omp_run_command = Sys.getenv('OMP_RUN_COMMAND')


numerical_methods = c('rk2', 'leapfrog_nonlinear')

for(numerical_method in numerical_methods){

    run_command = paste0(omp_run_command, ' ./BP1_testcases ', numerical_method, ' 1.0 0.019 > outfile.log')
    system(run_command)

    ERR_TOL = 1.0e-02

    # Get time-series at 2 locations
    canonical_ts = readLines('../test_repository/BP01-DmitryN-Single_wave_on_simple_beach/canonical_ts.txt')
    canonical_ts = canonical_ts[-(1:5)]
    ts_ab = read.table(text=canonical_ts, sep="\t")

    # Get the model results
    sink(tempfile())
    md = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])
    x = md[[1]] #get_all_recent_results()
    sink()
    d = 1.0 # This is set in the model code
    g = 9.8 

    #pdf('Model-vs-data.pdf', width=20, height=10)
    #pdf(paste0('Model-vs-data_', numerical_method, '.pdf'), width=20, height=10)
    png(paste0('Model-vs-data-at-two-sites_', numerical_method, '.png'), width=10, height=6, units='in', res=300)
    par(mfrow=c(2,1))
    par(mar=c(4,4,2,1))
    # x/d = 0.25
    plot(ts_ab[,1:2], xlab='Time x sqrt(d/g)', ylab='Stage (m)', t='l', 
         main='Timeseries at location x/d ~ 0.25; Analytical in black, SWALS model in red.')
    k = which.min(abs(x$gauges$lon - 0.25))
    points(x$time * sqrt(g/d), x$gauges$time_var$stage[k,], t='l', col='red')

    ## Quick test
    model_at_obs_t = approx(x$time * sqrt(g/d), x$gauges$time_var$stage[k,], xout=ts_ab[,1])
    # Ignore points in the quite-shallow part where the model goes dry but measurements do not
    k = which(ts_ab[,2] > -0.005)
    err_stat = mean(abs(model_at_obs_t$y[k] - ts_ab[k,2]))/diff(range(ts_ab[k,2]))
    if(err_stat < ERR_TOL){
        print('PASS')
    }else{
        print('FAIL')
    }

    # x/d = 9.95
    plot(ts_ab[,3:4], xlab='Time x sqrt(d/g)', ylab='Stage (m)', t='l', 
         main='Timeseries at location x/d ~ 9.95; Analytical in black, SWALS model in red.')
    k = which.min(abs(x$gauges$lon - 9.95))
    points(x$time * sqrt(g/d), x$gauges$time_var$stage[k,], t='l', col='red')

    ## Quick test
    model_at_obs_t = approx(x$time * sqrt(g/d), x$gauges$time_var$stage[k,], xout=ts_ab[,3], rule=2)
    # Ignore points in the quite-shallow part where the model goes dry but measurements do not
    k = 1:480
    err_stat = mean(abs(model_at_obs_t$y[k] - ts_ab[k,4]))/diff(range(ts_ab[k,4]))
    if(err_stat < ERR_TOL){
        print('PASS')
    }else{
        print('FAIL')
    }

    dev.off()

    ## Profiles
    canonical_profiles = read.table('../test_repository/BP01-DmitryN-Single_wave_on_simple_beach/canonical_profiles.txt', skip=5,header=FALSE)
    profile_times = seq(35, 70, by=5)
    png(paste0('Model-vs-data-canonical-profiles_', numerical_method, '.png'), width=10, height=6, units='in', res=300)
    par(mfrow=c(3,3))
    par(mar=c(4,4,3,1))
    for(i in 1:length(profile_times)){
        tind = which.min(abs(x$time * sqrt(g/d) - profile_times[i]))
        model_time_rescaled = x$time[tind] * sqrt(g/d)
        plot(canonical_profiles[,1], canonical_profiles[,i+1],t='l', cex.main=1.2, xlab='x (m)', ylab='Stage (m)', 
             main=paste0( 'T = ', profile_times[i], 'sqrt(d/g);  Numerical T = ', round(model_time_rescaled, 2), 'sqrt(d/g)'))
        points(x$gauges$lon, x$gauges$time_var$stage[,tind], t='l', col='red')
    }
    legend('bottomright', c('Analytical', 'SWALS model'), col=c('black', 'red'), lty=c(1,1), pch=c(NA, NA), cex=1.2)

    dev.off()
}

