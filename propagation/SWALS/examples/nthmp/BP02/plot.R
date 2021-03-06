source('../../../plot.R')
# Since we call SWALS from inside R, let's pass the openmp commands.
omp_run_command = Sys.getenv('OMP_RUN_COMMAND')

ERR_TOL = 1.0e-02

test_cases = c('a', 'b', 'c')
test_cases_caps = c('A', 'B', 'C')

timestepping_method = 'linear' 
#timestepping_method = 'leapfrog_nonlinear'
#timestepping_method = 'rk2'

show_obs = TRUE #FALSE

for(i in 1:length(test_cases)){

    run_command = paste0(omp_run_command, ' ./BP2_testcases case', test_cases_caps[i], ' ', timestepping_method, ' > outfile.log')
    system(run_command)

    x = get_all_recent_results(quiet=TRUE)

    analytical = read.table(paste0('../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', 
        test_cases[i], '_analytical.txt'), skip=6)
    obs = read.table(paste0('../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', test_cases[i], 
        '.txt'), skip=6)

    gauge_names = c(paste0('G', 4:10), 'Wall')

    # Plot the geometry
    if(i == 1){
        png(paste0('solution_geometry_caseA_', timestepping_method, '.png'), width=8, height=6, units='in', res=300)
        l = length(x$xs)
        plot(x$xs[1:(l-1)], x$elev0[1:(l-1),10], t='l', xlab='Along-beach distance (m)', ylab='Vertical coordinate (m)',
            ylim=c(-0.25, 0.25), lwd=2)
        abline(h=0, col='red', lty='dashed')
        legend('topleft', c('Beach profile', 'Initial water surface', 'Gauge location'),
            lwd=c(2,1,1), lty=c('solid', 'dashed', 'dotted'), col=c('black', 'red', 'orange'),
            bg='white')
        abline(v=x$gauges$lon[4:11], col='orange', lty='dotted')
        text(x$gauges$lon[4:11], c(rep(0.1, 7), 0.15), gauge_names)
        dev.off()
    }


    # Adjust time offset
    x$time = x$time + analytical[1,1]

    #pdf(paste0('solution', test_cases_caps[i], '.pdf'), width=15, height=8)
    png(paste0('solution', test_cases_caps[i], '_', timestepping_method, '.png'), width=15, height=8, units='in', res=300)
    par(mfrow=c(3,3))
    par(mar=c(5,5,3,2))
    plot_ylim = range(analytical[,2:9])
    plot_ylim = range(c(plot_ylim, range(x$gauges$time_var$stage)))
    for(ii in 2:ncol(analytical)){
        plot(analytical[,1], analytical[,ii],t='l', ylim=plot_ylim, 
            xlab='Time (s)', ylab='Stage (m)', lwd=4, cex.axis=1.7, cex.lab=2)
        if(show_obs & ii < ncol(analytical)) points(obs[,1], obs[,ii], t='l', col='green') 
        points(x$time, x$gauges$time_var$stage[ii+2,],t='l', col='red', lwd=1)
        if(ii == 2){
            if(show_obs){
                legend('topright', c('Analytical', 'Model', 'Experiment'), col=c('black', 'red', 'green'), lwd=c(4,1,1), cex=2)
            }else{
                legend('topright', c('Analytical', 'Model'), col=c('black', 'red'), lwd=c(4,1), cex=3, bty='n')
            }
        }
        title(gauge_names[ii-1], cex.main=3)

        # Quick test that error are small
        model_at_analytical = approx(x$time, x$gauges$time_var$stage[ii+2,], xout=analytical[,1])
        err_stat = mean(abs(model_at_analytical$y - analytical[,ii]))/diff(range(analytical[,ii]))
        if(err_stat < ERR_TOL){
            print('PASS')
        }else{
            print('FAIL')
        }
    }
    dev.off()
}
