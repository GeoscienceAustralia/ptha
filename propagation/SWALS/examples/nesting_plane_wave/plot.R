source('../../plot.R')

# Since we call SWALS from inside R, let's pass the openmp commands.
omp_run_command = Sys.getenv('OMP_RUN_COMMAND')

ts_methods = c('rk2', 'midpoint', 'linear', 'leapfrog_nonlinear')

for(ts_method in ts_methods){

    # Run the model
    system(paste0(omp_run_command, ' ./nesting_reflection ', ts_method, ' > outfile.log'))

    # Read the recent data
    my_dir = rev(Sys.glob('OUTPUTS/RUN*'))[1]
    x = lapply(as.list(Sys.glob(paste0(my_dir, '/RUN*'))), get_all_recent_results)

    # When should the wave reach the start again?
    # When it has travelled the domain length
    wave_speed = sqrt(9.8 * 100)
    domain_length = 100000
    time_to_cycle = domain_length/wave_speed
    nd = length(x)
    tind = which.min(abs(x[[nd]]$time - time_to_cycle))
    yind = round(dim(x[[nd]]$stage)[2]/2)
    # Given the write-out time won't be exactly "time_to_cycle", in theory
    # we can correct for the time-offset using a space offset
    x_offset = (x[[nd]]$time[tind] - time_to_cycle)*wave_speed


    # Plot the solutions
    #pdf(paste0('cycle_solution_', ts_method, '.pdf'), width=7, height=4)
    png(paste0('cycle_solution_', ts_method, '.png'), width=7, height=4, units='in', res=300)

    plot_ylim = range(c(range(x[[nd]]$stage[,yind,1]), range(x[[nd]]$stage[,yind,tind])))

    plot(x[[nd]]$xs, x[[nd]]$stage[,yind,tind], t='o', 
         main='Waveform at start (red) and after ~ 1 traverse of full domain (black). \n Plot shows inner domain only',
         xlab='x', ylim=plot_ylim,
         ylab='stage (m)')
    points(x[[nd]]$xs + x_offset, x[[nd]]$stage[,yind,1], t='l', col='red')

    dev.off()

    amp_loss_fraction = diff(range(x[[nd]]$stage[,yind,tind]))/diff(range(x[[nd]]$stage[,yind,1]))

    png(paste0('cycle_solution_relative_', ts_method, '.png'), width=7, height=4, units='in', res=300)
    plot(x[[nd]]$xs, x[[nd]]$stage[,yind,tind], t='o', ylim=plot_ylim,
         main=paste0('Waveform range at end as percent of start (black) = ', signif(amp_loss_fraction, 4)),
         xlab='x', 
         ylab='stage (m)')
    points(x[[nd]]$xs + x_offset, x[[nd]]$stage[,yind,1], t='l', col='red')
    dev.off()

    # Check the error
    s1 = x[[nd]]$stage[,yind,tind]
    s0 = x[[nd]]$stage[,yind, 1]

    k = which(abs(x[[nd]]$xs) < 10000)

    # Interpolate the initial waveform, but with a time-offset
    initial_s_shifted = approx(x[[nd]]$xs + x_offset, s0, xout=x[[nd]]$xs[k])

    err = sum( (s1[k]-initial_s_shifted$y)^2)/sum(s1[k]^2 + initial_s_shifted$y^2)

    png(paste0('cycle_solution_error_', ts_method, '.png'), width=7, height=4, units='in', res=300)
    plot(x[[nd]]$xs[k], s1[k] - initial_s_shifted$y, main='Error in the central part of the domain')
    dev.off()
    
    err_tol = 0.01

    if(err < err_tol){
        print(c('PASS', err))
    }else{
        print(c('FAIL', err))
    }
}
