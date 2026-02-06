source('../../../plot.R')
# Since we call SWALS from inside R, let's pass the openmp commands.
omp_run_command = Sys.getenv('OMP_RUN_COMMAND')

ERR_TOL = 1.0e-02

test_cases = c('a', 'b', 'c')
test_cases_caps = c('A', 'B', 'C')

timestepping_method = commandArgs(trailingOnly=TRUE)[1]
stopifnot(any(timestepping_method == c('midpoint', 'rk2', 'leapfrog_nonlinear', 'cliffs')))

show_obs = TRUE

for(i in 1:length(test_cases)){

    run_command = paste0(omp_run_command, ' ./BP5_testcases case', test_cases_caps[i], ' ', timestepping_method, ' > outfile.log')
    system(run_command)

    md = get_multidomain(sort(Sys.glob('OUTPUTS/RUN*'), decreasing=TRUE)[1])
    x = md[[2]]

    # Hack to make 'x' include all the gauges
    tmp = merge_multidomain_gauges(md)
    x$gauges = tmp

    analytical = read.table(paste0('../../nthmp/test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', 
        test_cases[i], '_analytical.txt'), skip=6)
    obs = read.table(paste0('../../nthmp/test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', test_cases[i], 
        '.txt'), skip=6)

    gauge_names = c(paste0('G', 4:10), 'Wall')

    yind  = floor(dim(x$elev0)[2] / 2) + 1

    # Plot the geometry
    if(i == 1){
        png(paste0('solution_geometry_caseA_', timestepping_method, '.png'), width=8, height=6, units='in', res=300)
        l = length(x$xs)
        plot(x$xs[1:(l-1)], x$elev0[1:(l-1),yind], t='l', xlab='Along-beach distance (m)', ylab='Vertical coordinate (m)',
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
    x$gauges$time = x$gauges$time + analytical[1,1]

    #pdf(paste0('solution', test_cases_caps[i], '.pdf'), width=15, height=8)
    png(paste0('solution', test_cases_caps[i], '_', timestepping_method, '.png'), width=15, height=8, units='in', res=300)
    par(mfrow=c(3,3))
    par(mar=c(5,5,3,2))
    #plot_ylim = range(analytical[,2:9])
    #plot_ylim = range(c(plot_ylim, range(x$gauges$time_var$stage)))
    if(test_cases_caps[i] == 'A'){
        # Case A plots better with a smaller yrange 
        plot_ylim = c(-0.01, 0.03)
    }else{
        plot_ylim = c(-0.05, 0.2)
    }
    for(ii in 2:ncol(analytical)){
        plot(analytical[,1], analytical[,ii],t='l', ylim=plot_ylim, 
            xlab='Time (s)', ylab='Stage (m)', lwd=4, cex.axis=1.7, cex.lab=2, col='white')
        include_obs = (show_obs & ii < ncol(analytical))
        if(include_obs) points(obs[,1], obs[,ii], t='l', col='green') 

        model_time = x$gauges$time
        model_stage = x$gauges$time_var$stage[ii+2,]
        points(model_time, model_stage,t='l', col='red', lwd=1)
        if(ii == 2){
            if(show_obs){
                legend('topright', 
                       c(paste0('Model (', timestepping_method, ')'), 'Experiment'), 
                       col=c('red', 'green'), lwd=c(1,1), cex=2, bty='n')
            }else{
                legend('topright', 
                       c(paste0('Model (', timestepping_method, ')')), 
                       col=c('red'), lwd=c(1), cex=3, bty='n')
            }
        }
        title(gauge_names[ii-1], cex.main=3)

        # Include some PASS/FAIL tests (not a substitute for looking at the results!)
        if(include_obs){
            # Quick model-vs-data stat
            ts = seq(270, 290, by=0.01)
            model_ts = approx(model_time, model_stage, xout=ts)$y
            obs_ts = approx(obs[,1], obs[,ii], xout=ts)$y
            err_stat = mean(abs(model_ts - obs_ts))

            # Ad-hoc error tolerance for cases A, B, C
            err_thresh = c(0.0015, 0.005,  0.01)[i]
            if(err_stat < err_thresh){
                print('PASS')
            }else{
                print(c('FAIL', timestepping_method, test_cases_caps[i], ii, err_stat, err_thresh))
            }
        }
    }
    dev.off()
}
