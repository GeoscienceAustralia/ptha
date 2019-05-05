source('../../plot.R')

ts_methods = c('rk2', 'linear')

for(ts_method in ts_methods){

    # Run the model
    system(paste0('./nesting_reflection ', ts_method, ' > outfile.log'))

    # Read the recent data
    my_dir = rev(Sys.glob('OUTPUTS/RUN*'))[1]
    x = lapply(as.list(Sys.glob(paste0(my_dir, '/RUN*'))), get_all_recent_results)

    # When should the wave reach the start again?
    # When it has travelled the domain length
    wave_speed = sqrt(9.8 * 100)
    domain_length = 100000
    time_to_cycle = domain_length/wave_speed
    tind = which.min(abs(x[[2]]$time - time_to_cycle))
    yind = round(dim(x[[2]]$stage)[2]/2)
    # Given the write-out time won't be exactly "time_to_cycle", in theory
    # we can correct for the time-offset using a space offset
    x_offset = (x[[2]]$time[tind] - time_to_cycle)*wave_speed

    # Plot the solutions
    pdf(paste0('cycle_solution_', ts_method, '.pdf'), width=7, height=5)

    plot(x[[2]]$xs, x[[2]]$stage[,yind,tind], t='o', 
         main='Waveform at start (red) and after ~ 1 traverse of full domain (black). \n Plot shows inner domain only (3x refinement)',
         xlab='x', 
         ylab='stage (m)')
    points(x[[2]]$xs + x_offset, x[[2]]$stage[,yind,1], t='l', col='red')


    # Check the error
    s1 = x[[2]]$stage[,yind,tind]
    s0 = x[[2]]$stage[,yind, 1]

    k = which(abs(x[[2]]$xs) < 10000)

    # Interpolate the initial waveform, but with a time-offset
    initial_s_shifted = approx(x[[2]]$xs + x_offset, s0, xout=x[[2]]$xs[k])

    err = sum( (s1[k]-initial_s_shifted$y)^2)/sum(s1[k]^2 + initial_s_shifted$y^2)

    plot(x[[2]]$xs[k], s1[k] - initial_s_shifted$y, main='Error in the central part of the domain')

    dev.off()
    
    err_tol = 0.01

    if(err < err_tol){
        print(c('PASS', err))
    }else{
        print(c('FAIL', err))
    }
}
