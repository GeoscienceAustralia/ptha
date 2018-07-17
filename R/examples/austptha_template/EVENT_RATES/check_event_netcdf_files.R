#
# Quick check that the event netcdf files have consistent rates
# among uniform/stochastic/variable-uniform slip
#
library(rptha)

config = new.env()
source('config.R', local=config)

check_source<-function(uniform_slip_tsunami_file, stochastic_slip_tsunami_file, 
    variable_uniform_slip_tsunami_file, rates_should_be_zero=FALSE){

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

    # Unit source statistics
    uss_file = Sys.glob(paste0(dirname(uniform_slip_tsunami_file), '/unit_source_statistics_*.nc'))
    if(length(uss_file) != 1) stop('matching multiple unit-source files')
    uss = read_table_from_netcdf(uss_file)


    # Check that the event rates in the _tsunami and regular files are the
    # same, up to floating point storage
    for(var in c('rate_annual', 'rate_annual_lower_ci', 'rate_annual_upper_ci')){
        for(nme in c('unif', 'stoc', 'vuni')){
            r1 = ncvar_get(fids[[nme]], var)
            r2 = ncvar_get(fids_t[[nme]], paste0('event_', var))

            if(rates_should_be_zero){
                if(!all(r1==0 & r2==0)){
                    stop(paste0('Rates should be zero, but are not, in ', var, ' ', nme))
                }
            }else{
                # Allow for both 'relative' and 'absolute' errors
                err = abs(r1 - r2)
                if(!all( (err < 1.0e-06*r1)|(err < 1e-12))){
                    stop(paste0('rates do not match ', var, ' ', nme ))
                }
            }
        }
    }
    rm(r1, r2)

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
        # Skip if rates are all zero
        if(rates_should_be_zero) next

        for(nme in c('unif', 'stoc', 'vuni')){
            r1 = ncvar_get(fids[[nme]], var)
            r2 = ncvar_get(fids[[nme]], paste0('variable_mu_', var))

            # Note for quantiles, arbitrary changes are possible in rare cases. 
            # For instance, slight changes in the weights can shift which curve is 
            # selected in the quantile computation -- and the difference in the selected
            # curves could be arbitrarly large (e.g. if one curve assigns a rate of zero)
            # This is not super-common but does happen
            # Furthermore, we apply different bias-correction methods for variable_mu and
            # regular events (for stochastic and variable_uniform), which adds
            # greatly to the variability on the ratios for individual events
            ratios = (r1+1e-200)/(r2+1e-200)
            print(paste0(c('Event rate change due to variable mu, ', var, ' ', nme, ' :', 
                round(c(min(ratios), median(ratios), max(ratios)),4)), collapse=" "))
        }
    }
    rm(r1, r2, ratios)

    # Check that when appropriately summed, the variable-uniform and stochastic rates
    # are the same as the uniform rates
    for(var in c('rate_annual', 'rate_annual_lower_ci', 'rate_annual_upper_ci')){
        # Skip if all rates are zero
        if(rates_should_be_zero) next

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
    rm(r_unif, unif_event_row, r1, rate_by_row)

    # Sanity check on the ordering quantiles
    # Note the quantiles refer to the cumulative distribution function -- while
    # the 'event specific' values refer to the derivative of this. The
    # derivatives need not have a strict ordering (although there is certainly
    # a tendency for that.
    for(variable_mu in c('', 'variable_mu_')){
        for(nme in c('unif', 'stoc', 'vuni')){
            # Read the lowest quantile, and compare with the second lowest. Then
            # compare second lowest to third lowest. etc.
            var1 = paste0(variable_mu, 'rate_annual_lower_ci')
            oldvar_vals = ncvar_get(fids[[nme]], var1)
            other_vars = paste0(variable_mu, # Note the following variables are ordered
                c('rate_annual_16pc', 'rate_annual_median', 'rate_annual_84pc', 'rate_annual_upper_ci'))
            for(newvar in other_vars){
                newvar_vals = ncvar_get(fids[[nme]], newvar)
                # Since the values should be ordered by Mw, and the quantiles
                # are initially defined in terms of exceedance rates, a
                # cumulative sum of the rates should retain ordering.
                n1 = cumsum(rev(newvar_vals))
                n2 = cumsum(rev(oldvar_vals))
                if(!all(n1 >= n2)){
                    stop(paste0('ordering problem in quantiles near ', newvar, ' ', nme))
                }
                oldvar_vals = newvar_vals
                rm(n1, n2)
            }
        }
    }
    rm(oldvar_vals, newvar_vals)

    # Check weight with nonzero rate is sensible
    for(nme in c('unif', 'stoc', 'vuni')){

        weight_with_nonzero_rate = ncvar_get(fids[[nme]], 'weight_with_nonzero_rate')
        event_rate = ncvar_get(fids[[nme]], 'rate_annual')
        mw = round(ncvar_get(fids[[nme]], 'Mw'), 3) # Deal with imperfect netcdf floating point storage

        print(paste0('Range weight_with_nonzero_rate: ',
                    min(weight_with_nonzero_rate), ' ', 
                    max(weight_with_nonzero_rate), collapse=" "))

        if(rates_should_be_zero){
            if(!all(weight_with_nonzero_rate == 0)){
                stop(paste0('weights_with_nonzero_rate should all be zero, but they are not, in ', nme))
            }
            next
        }

        k = which(weight_with_nonzero_rate == 0)
        k2 = which(event_rate == 0)
        if(!all(k == k2)){
            stop('conflict between weight_with_nonzero_rate and rate_annual')
        }
 
        # We should not have events with Mw > config$MAXIMUM_ALLOWED_MW_MAX
        if(uss$rake[1] == 90){
            mw_max_limit = config$MAXIMUM_ALLOWED_MW_MAX
        }else if(uss$rake[1] == -90){
            mw_max_limit = config$MAXIMUM_ALLOWED_MW_MAX_NORMAL
        }else{
            stop('Rake is neither pure thrust nor pure normal')
        } 
        mw_min_limit = config$MINIMUM_ALLOWED_MW_MAX

        # Mw-max is respected
        k = which(mw > mw_max_limit)
        if(length(k) > 0){
            if(!all(weight_with_nonzero_rate[k] == 0)){
                stop('weight_with_nonzero_rate is > 0 for events beyond the upper limit')
            }
        } 

        # Check Mw-min is respected
        # Most events with mw < mw_min should haev weight = 1. However, it is
        # possible to have events with mw<mw_min that have zero probability of
        # occurring, if they are near a divergent section of the plate boundary
        # based on Bird's model. For example, this occurs in a couple of places
        # in western Aleutians. 
        k = which(mw < mw_min_limit & (event_rate > 0))
        stopifnot(length(k) > 0)
        if(!all(abs(weight_with_nonzero_rate[k] - 1) < 1.0e-06)){
            stop('weight_with_nonzero_rate is not equal to 1 below the mw_max_lower_limit')
        }
    }


    # Check that the quantile adjustment is sensible (conforms with
    # interpolation functions in config)
    for(nme in c('stoc', 'vuni')){

        if(rates_should_be_zero) next      
 
        uniform_event_row = ncvar_get(fids[[nme]], 'uniform_event_row')
        event_slip_string = ncvar_get(fids[[nme]], 'event_slip_string')
        peak_slip = sapply(event_slip_string, f<-function(x) max(as.numeric(strsplit(x,'_')[[1]])))
        event_rate = ncvar_get(fids[[nme]], 'rate_annual')
        event_rate_mu_vary = ncvar_get(fids[[nme]], 'variable_mu_rate_annual')

        unique_uniform_event_row = unique(uniform_event_row)

        # Get the bias adjustment functions
        if(uss$rake[1] == 90){
            if(nme == 'stoc'){
                bias_adjuster_fixed_mu = config$peak_slip_bias_adjustment_stochastic
                bias_adjuster_vary_mu = config$peak_slip_bias_adjustment_stochastic_mu_vary
            }else if(nme == 'vuni'){
                bias_adjuster_fixed_mu = config$peak_slip_bias_adjustment_variable_uniform
                bias_adjuster_vary_mu = config$peak_slip_bias_adjustment_variable_uniform_mu_vary
            }else{
                stop(paste0('invalid nme ', nme))
            }
        }else{
            # For non-thrust events, we do not do bias adjustment or treat
            # variable shear modulus
            bias_adjuster_fixed_mu<-function(x) 1 + x*0
            bias_adjuster_vary_mu<-function(x) 1 + x*0
        }

        # Do the test
        for(uuer in unique_uniform_event_row){
            k = which(uniform_event_row == uuer)
            quantiles = rank(peak_slip[k], ties='first')/(length(k)+1)

            expected_ratios_fixed_mu = bias_adjuster_fixed_mu(quantiles)
            expected_ratios_vary_mu = bias_adjuster_vary_mu(quantiles)

            overall_rate_fixed_mu = sum(event_rate[k])
            overall_rate_vary_mu = sum(event_rate_mu_vary[k])

            err1 = max(abs(overall_rate_fixed_mu * expected_ratios_fixed_mu/sum(expected_ratios_fixed_mu) - event_rate[k]))
            err2 = max(abs(overall_rate_vary_mu * expected_ratios_vary_mu/sum(expected_ratios_vary_mu) - event_rate_mu_vary[k]))

            if(!(err1 <= 1.0e-06*overall_rate_fixed_mu)){
                stop(paste0('issue with bias adjustment fixed mu ', uuer))
            }
    
            if(!(err2 <= 1.0e-06*overall_rate_vary_mu)){
                stop(paste0('issue with bias adjustment variable mu ', uuer))
            }
        }
        
    }

    for(i in 1:length(fids)) nc_close(fids[[i]])
    for(i in 1:length(fids_t)) nc_close(fids_t[[i]])
}

for(i in 1:length(config$all_source_uniform_slip_tsunami)){
    print(' ')
    print('######################')
    source_name = basename(dirname(dirname(config$all_source_uniform_slip_tsunami[i])))
    print(source_name)

    # For 'defunct' sources, we just check that all their rates are zero
    if(source_name %in% c('timor', 'flores', 'macquarienorth')){
        rates_should_be_zero=TRUE
        print('   rates should be zero!')
    }else{
        rates_should_be_zero=FALSE
        print('   regular source')
    }
    print(' ')
    check_source(config$all_source_uniform_slip_tsunami[i], 
        config$all_source_stochastic_slip_tsunami[i], 
        config$all_source_variable_uniform_slip_tsunami[i],
        rates_should_be_zero=rates_should_be_zero)
}

