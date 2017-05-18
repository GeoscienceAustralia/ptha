test_cases = c('a', 'b', 'c')
test_cases_caps = c('A', 'B', 'C')

timestepping_method = 'linear' #'rk2'

for(i in 1:length(test_cases)){

    run_command = paste0('./BP2_testcases case', test_cases_caps[i], ' ', timestepping_method)
    system(run_command)

    source('../../../plot.R')

    x = get_all_recent_results()

    analytical = read.table(paste0('../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', 
        test_cases[i], '_analytical.txt'), skip=6)
    obs = read.table(paste0('../test_repository/BP02-DmitryN-Solitary_wave_on_composite_beach_analytic/ts3', test_cases[i], 
        '.txt'), skip=6)

    # Adjust time offset
    x$time = x$time + analytical[1,1]

    pdf(paste0('solution', test_cases_caps[i], '.pdf'), width=15, height=8)
    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2))
    plot_ylim = range(analytical[,2:9])
    plot_ylim = range(c(plot_ylim, range(x$gauges$time_var$stage)))
    for(ii in 2:ncol(analytical)){
        plot(analytical[,1], analytical[,ii],t='l', ylim=plot_ylim, xlab='x', ylab='y', lwd=4)
        if(ii < ncol(analytical)) points(obs[,1], obs[,ii], t='l', col='green') 
        points(x$time, x$gauges$time_var$stage[ii+2,],t='l', col='red', lwd=1)
        if(ii == 2){
            legend('topright', c('Analytical', 'Model', 'Experiment'), col=c('black', 'red', 'green'), lwd=c(4,1,1))
        }
    }
    dev.off()
}
