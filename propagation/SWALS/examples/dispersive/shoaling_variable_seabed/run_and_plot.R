#
# Run and plot the "Shoaling with a variable seabed' problem 
#
numerical_schemes = c('midpoint', 'rk2', 'leapfrog_nonlinear', 'cliffs') #c('midpoint', 'leapfrog_nonlinear', 'rk2', 'cliffs')
resolutions = c(0.04, 0.01) #c(0.04, 0.02, 0.01) #c(0.1, 0.05)

data_store = list()
for(nm in numerical_schemes){
    plot_started = FALSE # After we start plotting this will be TRUE (and then we add to the existing plot)

    for(long_dimension in c('y', 'x')){
    
        # Run each resolution
        for(res in resolutions){
            run_command = paste0('OMP_NUM_THREADS=6 OMP_PROC_BIND=true ./shoaling_variable_seabed ', nm, ' ', res, ' ', long_dimension)
            system(run_command)
        }

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

    data_store[[nm]] = vector(mode='list', length=length(matching_dirs))
    names(data_store[[nm]]) = matching_dir_labels

    #
    # Plot like second panel of Figure 10 in Madsen et al (2008)
    #
    png(paste0('Solution_', nm, '.png'), width=8, height=4, units='in', res=200)
    library(ncdf4)
    for(i in 1:length(matching_dirs)){
        # Read outputs
        nc_file = Sys.glob(paste0(matching_dirs[i], '/RUN*/Grid*.nc'))
        fid = nc_open(nc_file)
        stg = ncvar_get(fid, 'stage')
        x = ncvar_get(fid, 'x')
        y = ncvar_get(fid, 'y')
        time = ncvar_get(fid, 'time')
        max_stage = ncvar_get(fid, 'max_stage')
        min_stage = ncvar_get(fid, 'max_stage')
        min_stage = ncvar_get(fid, 'min_stage')
        elev = ncvar_get(fid, 'elevation0')
        nc_close(fid)

        # Plot depends on whether the long dimension is 'y' or 'x'
        if(grepl('long_dimension_y', matching_dirs[i])){
            k = which(elev[3,] < max(elev)-1e-03) # Used to avoid walls (1.0 m)
            x_k = y[k]
            max_stg_k = max_stage[3,k]
            min_stg_k = min_stage[3,k]
            elev_k = elev[3,k]  
            initial_stage_k = stg[3,k,1]
            plotme = FALSE
        }else if(grepl('long_dimension_x', matching_dirs[i])){
            k = which(elev[,3] < max(elev)-1e-03) # Used to avoid walls (1.0 m)
            x_k = x[k]
            max_stg_k = max_stage[k,3]
            min_stg_k = min_stage[k,3]
            elev_k = elev[k,3]
            initial_stage_k = stg[k,3,1]
            plotme = TRUE
        }else{
            stop('unknown long_dimension')
        }

        if(plotme){
            if(!plot_started){
                plot_started = TRUE
                plot(x_k, max_stg_k, t='l', ylim=c(-0.25, 0.15), col=i, 
                    main=paste0('Water surface envelope: ', nm), xlab='x', ylab='Stage (m)')
                grid()
            }else{
                points(x_k, max_stg_k, t='l', col=i)
            }
            points(x_k, min_stg_k, t='l', col=i)
        }

        data_store[[nm]][[matching_dir_labels[i]]] = data.frame(x_k = x_k, max_stg_k = max_stg_k, min_stg_k = min_stg_k)

    }

    # Whispers 3D solution digitized from 2 figures, so 2 parts
    d1 = read.csv('solutions_up_to_30m_A2.csv')
    d2 = read.csv('solutions_beyond_30m.csv')
    d1 = rbind(d1, d2)
    w3d_col = 14
    w3d_lty = 'dotted'
    points(d1[,1], d1[,2], t='l', lty=w3d_lty, col=w3d_col)
    points(d1[,1], d1[,3], t='l', lty=w3d_lty, col=w3d_col)

    # Only plot 'x' solutions (since the tests will check that they are the same as 'y' solutions)
    kld = which(grepl('long_dimension_x', matching_dir_labels))

    include_ms_solution = FALSE
    if(include_ms_solution){
        # Also plot the MS solution from Coulaud et al (2025). The non-default SWALS dispersion that includes BOTH
        # "Peregrine" type terms in flux form gives similar results to MS, which is also based on that
        # flux form of the equations (but with extra terms). Filippini et al (2015) suggest that 
        # formulation in flux or velocity form may be the most important difference between schemes in the nonlinear regime,
        # which might explain the similarity
        m1 = read.csv('solutions_up_to_30m_MS.csv')
        m2 = read.csv('solutions_beyond_30m_MS.csv')
        ms = rbind(m1, m2)
        ms_col = 2
        ms_lty = 'dashed'
        points(ms[,1], ms[,2], t='l', lty=ms_lty, col=ms_col)
        points(ms[,1], ms[,3], t='l', lty=ms_lty, col=ms_col)
        legend('bottomright', 
            c(paste0(matching_dir_labels[kld], 'm'), 'Whispers3D (Coulaud et al 2025)', 'MS (Coulaud et al 2025)'), 
            col=c((1:length(matching_dir_labels))[kld], w3d_col, ms_col), 
            lty=c(rep('solid', length(kld)), w3d_lty, ms_lty), 
            pch=NA,
            bty='n', cex=1.0)
    }else{
        legend('bottomright', 
            c(paste0(matching_dir_labels[kld], 'm'), 'Whispers3D (Coulaud et al 2025)'), 
            col=c((1:length(matching_dir_labels))[kld], w3d_col), 
            lty=c(rep('solid', length(kld)), w3d_lty), 
            pch=NA,
            bty='n', cex=1.0)
    }

    dev.off()
    #
    # Check that the numerical solutions computed along x/y are consistent with each other
    #
    for(res in resolutions){
        # Find runs that are the same except for their long_dimension
        ii = which(endsWith(names(data_store[[nm]]), paste0('_', res))); stopifnot(length(ii) == 2)
        # Should have minimal solution change when adjusting x/y dim
        err_stat = max(abs(data_store[[nm]][[ ii[1] ]] - data_store[[nm]][[ ii[2] ]]))
        if(all(err_stat < 1e-06)){
            print('PASS')
        }else{
            print(paste0('FAIL (long-x vs long-y) ', nm, ' ', res, ' ', err_stat))
        }

        # Compare models with the observed envelope 
        model_x = data_store[[nm]][[ ii[1] ]]$x_k
        model_max = data_store[[nm]][[ ii[1] ]]$max_stg_k
        model_min = data_store[[nm]][[ ii[1] ]]$min_stg_k
        data_upper = suppressWarnings(approx(d1[,1], d1[,2], xout=model_x )$y) # Ties --> warnings
        data_lower = suppressWarnings(approx(d1[,1], d1[,3], xout=model_x )$y)
        
        # Only use values in the range of the data
        dr = which((model_x > min(d1[,1])) & (model_x < max(d1[,1])))

        # Check there are enough values
        if(any(is.na(data_upper[dr]) | is.na(data_lower[dr]))){
            print(paste0('FAIL (problem interpolating data) ', nm, ' ', res))
        }else{
            print('PASS')
        }
        
        err_stat_upper = mean(abs(data_upper - model_max)[dr]) 
        err_stat_lower = mean(abs(data_lower - model_min)[dr]) 

        # Test with ad-hoc error thresholds
        if(err_stat_upper < 0.013){
            print('PASS')
        }else{
            print(paste0('FAIL (err_stat_upper) ', err_stat_upper, ' ', nm, ' ', res))
        }

        if(err_stat_lower < 0.005){
            print('PASS')
        }else{
            print(paste0('FAIL (err_stat_lower)', err_stat_lower, ' ', nm, ' ', res))
        }

    }

    # Plot the initial condition
    if(nm == numerical_schemes[1]){
        png('Initial_condition.png', width=8, height=6, units='in', res=200)
        plot(x_k, elev_k, t='l', ylim=c(-1, 0.2), col='brown', xlab='x (m)', ylab='Stage (m)', main='Initial condition', cex.main=2)
        points(x_k, initial_stage_k, t='l', col='blue')
        legend('bottomright', c('Bed', 'Initial stage'), lty=c(1,1), col=c('brown', 'blue'), pch=NA, bty='n', cex=2)
        dev.off()
    }

}
