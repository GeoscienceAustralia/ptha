#
# Quick check that the event netcdf files have consistent rates
# among uniform/stochastic/variable-uniform slip
#

config = new.env()
source('config.R', local=config)

check_source<-function(uniform_slip_tsunami_file, stochastic_slip_tsunami_file, variable_uniform_slip_tsunami_file){

    fids = list()
    fids_t = list()

    # Uniform slip
    fids$unif = nc_open(gsub('_tsunami', '', uniform_slip_tsunami_file), readunlim=FALSE)
    fids_t$unif = nc_open(uniform_slip_tsunami_file, readunlim=FALSE)

    # Stochastic slip
    fids$stoc = nc_open(gsub('_tsunami', '', stochastic_slip_tsunami_file), readunlim=FALSE)
    fids_t$stoc = nc_open(stochastic_slip_tsunami_file, readunlim=FALSE)

    # Variable uniform slip
    fids$vuni = nc_open(gsub('_tsunami', '', variable_uniform_slip_tsunami_file), readunlim=FALSE)
    fids_t$vuni = nc_open(variable_uniform_slip_tsunami_file, readunlim=FALSE)


    # Check that the event rates in the _tsunami and regular files are the
    # same, up to floating point storage
    for(var in c('rate_annual', 'rate_annual_lower_ci', 'rate_annual_upper_ci')){
        for(nme in c('unif', 'stoc', 'vuni')){
            r1 = ncvar_get(fids[[nme]], var)
            r2 = ncvar_get(fids_t[[nme]], paste0('event_', var))
            if(!all(abs(r1 - r2) < 1.0e-09)){
                stop(paste0('rates do not match ', var, ' ', nme ))
            }
        }
    }

    # Check that index strings match
    for(nme in c('unif', 'stoc', 'vuni')){
        eis1 = ncvar_get(fids[[nme]], 'event_index_string')
        eis2 = ncvar_get(fids_t[[nme]], 'event_index_string')
        if(!all(eis1 == eis2)){
            stop(paste0('event index strings do not match ', nme))
        }
    }

    # Check how much rates have changed
    for(var in c('rate_annual', 'rate_annual_lower_ci', 'rate_annual_upper_ci')){
        for(nme in c('unif', 'stoc', 'vuni')){
            r1 = ncvar_get(fids[[nme]], var)
            r2 = ncvar_get(fids[[nme]], paste0('variable_mu_', var))

            # Note for quantiles, arbitrary changes are possible in rare cases. 
            # For instance, slight changes in the weights can shift which curve is 
            # selected in the quantile computation -- and the difference in the selected
            # curves could be arbitrarly large (e.g. if one curve assigns a rate of zero)
            # This is not super-common but does happen
            ratios = (r1+1e-200)/(r2+1e-200)
            print(paste0(c('Event rate change due to variable mu, ', var, ' ', nme, ' :', 
                round(c(min(ratios), median(ratios), max(ratios)),4)), collapse=" "))
        }
    }

    # Check that when appropriately summed, the variable-uniform and stochastic rates
    # are the same as the uniform rates
    for(var in c('rate_annual', 'rate_annual_lower_ci', 'rate_annual_upper_ci')){
        r_unif = ncvar_get(fids$unif, var)
        for(nme in c('stoc', 'vuni')){
            # Stochastic    
            unif_event_row = ncvar_get(fids[[nme]], 'uniform_event_row')
            r1 = ncvar_get(fids[[nme]], var)
            rate_by_row = aggregate(r1, list(unif_event_row), sum)
            rate_by_row = rate_by_row$x

            if(!all(abs(r_unif-rate_by_row) <= 1e-06*r_unif)){
                stop(paste0('rates when summed do not match ', var, ' ', nme))
            }
        }
    }

    for(i in 1:length(fids)) nc_close(fids[[i]])
    for(i in 1:length(fids_t)) nc_close(fids_t[[i]])
}

for(i in 1:length(config$all_source_uniform_slip_tsunami)){
    print(' ')
    print('######################')
    print(basename(dirname(dirname(config$all_source_uniform_slip_tsunami[i]))))
    print(' ')
    check_source(config$all_source_uniform_slip_tsunami[i], config$all_source_stochastic_slip_tsunami[i], config$all_source_variable_uniform_slip_tsunami[i])
}

