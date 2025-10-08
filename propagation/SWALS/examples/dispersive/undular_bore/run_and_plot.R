#
# Run and plot the undular bore simulation at various resolutions
#
numerical_schemes = c('midpoint', 'leapfrog_nonlinear', 'rk2', 'cliffs')
resolutions = c(20, 10, 5, 2)

for(nm in numerical_schemes){
    
    # Run each resolution
    for(res in resolutions){
        run_command = paste0('OMP_NUM_THREADS=6 OMP_PROC_BIND=true ./undular_bore ', nm, ' ', res)
        system(run_command)
    }

    # Find each run
    matching_dirs = Sys.glob(paste0('OUTPUTS/', nm, '*/RUN*'))
    matching_dir_labels = basename(dirname(matching_dirs))

    # Get a good order
    o1 = order(as.numeric(sapply(matching_dir_labels, 
                                function(x){tmp = strsplit(x, '_')[[1]]; tmp[length(tmp)]}
            )), decreasing=TRUE)
    matching_dirs = matching_dirs[o1]
    matching_dir_labels = matching_dir_labels[o1]

    madsen_2008_bore = read.csv('madsen_2008_bore.csv')

    # Plot like second panel of Figure 10 in Madsen et al (2008)
    png(paste0('Time_1500_', nm, '.png'), width=8, height=6, units='in', res=200)
    library(ncdf4)
    for(i in 1:length(matching_dirs)){
        nc_file = Sys.glob(paste0(matching_dirs[i], '/RUN*/Grid*.nc'))
        fid = nc_open(nc_file)
        stg = ncvar_get(fid, 'stage')
        x = ncvar_get(fid, 'x')
        time = ncvar_get(fid, 'time')
        nc_close(fid)
        ind = which.min(abs(time-1500))
        stage = stg[,3,ind]
        if(abs(time[ind] - 1500) > 1) print('WARNING: output time varies from 1500 by more than 1 second')
        if(i == 1){
            plot(x, stage, t='l', ylim=c(-0.5, 3.5), xlim=c(19000, 22000), col=i, 
                main='Undular bore at different resolutions', xlab='x', ylab='Stage (m) @ t=1500')
            grid()
        }else{
            points(x, stage, t='l', col=i)
        }

        points(madsen_2008_bore[,1], madsen_2008_bore[,2], pch=1, cex=0.3, col='black')
    }
    legend('bottomleft', 
        c(paste0(matching_dir_labels, 'm'), 'Madsen et al. (2008, their Fig 10) with a more complex dispersive model.'), 
        col=c(1:length(matching_dirs), 1), 
        lty=c(rep(1, length(matching_dirs)), NA), 
        pch=c(rep(NA, length(matching_dirs)), 1),
        bty='n')
    dev.off()

    # Add in a PASS/FAIL test by comparison with Madsen
    test_xrange = seq(21250, 21600, by=10) # Leading part of bore
    test_SWALS = approx(x, stage, xout=test_xrange)$y # Due to plot ordering, this should be the highest-res solution
    test_Madsen = approx(madsen_2008_bore[,1], madsen_2008_bore[,2], xout=test_xrange)$y

    # Test 1 -- compare average abs deviation with Madsen.
    test_stat = mean(abs(test_SWALS - test_Madsen))
    err_tol = 0.25 # Ad-hoc choice
    #print(test_stat)
    if(test_stat < err_tol){
        print('PASS')
    }else{
        print(paste0('FAIL', nm, test_stat))
    }

    # Test 2 -- compare maxima with madsen
    test_stat_2 = max(test_SWALS) - max(test_Madsen)
    err_tol_2 = 0.25 # Ad-hoc choice
    if(abs(test_stat_2) < err_tol_2){
        print('PASS')
    }else{
        print(paste0('FAIL', nm, test_stat, test_stat_2))
    }
        
}
