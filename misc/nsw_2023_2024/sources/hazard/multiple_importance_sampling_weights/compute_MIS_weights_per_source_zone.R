#
# Investigations of the "optimal weights" for multiple importance samples.
#
# We have N separate Monte Carlo samples, using different importance functions, which
# are in principle optimised differently for different sites. We can estimate the hazard
# by combining the samples IS1, IS2, IS3 ... with 
#    multiple_importance_sampling_exrate = sum_j{ w_j * exrate(ISj) } 
# where exrate(ISj) is an exceedance-rate of interest for the j'th importance sample
# for a scenario frequency model of interest, and
#    sum_j { w_j } = 1
#
# If the weights are determined without "peeking" at the samples, then this will be an unbiased
# estimate of the exceedance-rate of interest.
#
# In the NSW study context, we can use different weights for different source-zones.
#
#

# Functions to access the PTHA18
find_ptha_results_file<-function(){
    get_ptha_results_file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
    get_ptha_results_file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
    ifelse(file.exists(get_ptha_results_file_nci), get_ptha_results_file_nci, get_ptha_results_file_home)
}
ptha18 = new.env()
source(find_ptha_results_file(), local=ptha18, chdir=TRUE)

# Importance sampling utilities
isu = new.env()
source('../importance_sampling_utilities.R', local=isu)

#source('fit_weights_by_latitude.R')

# Main worker function
compute_MIS_weights_per_source_zone<-function(ptha18, isu, 
    all_sample_dirs, all_samples, 
    model_ptha18_point_gauges,
    stage_exceedance_rate_thresholds){

    library(ncdf4)

    stopifnot(length(all_sample_dirs) == length(all_samples))

    # Get the RDS data for each importance sample
    all_samples_by_sources = lapply(all_sample_dirs, 
        function(x) readRDS(paste0(x, '/source_zone_data_backup.RDS')))
    names(all_samples_by_sources) = all_samples
        # Key variables: 
        #    all_samples_by_sources$source_zone$sampling_prob
        #    all_samples_by_sources$source_zone$event_rate_logic_tree_mean
        #    all_samples_by_sources$source_zone$N_MC
        # We can then use isu$analytical_exrate_and_MC_variance to get the variance
        # if we have the event_stage at a PTHA18 point.
        #
        # Do that calculation for every importance sample, for a good range of PTHA18
        # points offshore of NSW and the islands, and see what weights are suggested
        # at different return periods. Then use the results to figure out a weighting
        # approach.

    source_zones = names(all_samples_by_sources[[1]])
    # Check all samples have the same source zones
    # The code below assumes they do (although could be modified to a more general case)
    for(i in 1:length(all_samples_by_sources)){
        stopifnot(all( names(all_samples_by_sources[[i]]) == source_zones ))
    }

    ## Get exceedance-rate curves for every PTHA18 point -- sum over source zones.
    ## This is for the SUM OVER SOURCE ZONES. We use it purely to decide at which stage
    ## values the weights should be computed (aiming to keep a consistent exceedance-rate).
    #fid = nc_open(stage_exceedance_rate_curves_file)
    #ptha18_gauge_IDs_in_stage_exrate_file = ptha18$clean_gaugeID(ncvar_get(fid, 'gaugeID'))
    #ptha18_gauge_elev_in_stage_exrate_file = ncvar_get(fid, 'elev')
    #stochastic_slip_stage_exrate_table = ncvar_get(fid, 'stochastic_slip_rate_84pc') # 84th percentile curve
    #stochastic_slip_stages = ncvar_get(fid, 'stage') 
    #nc_close(fid)

    # Comptue the exrate and its MC variance for each source_zone, site, importance sample, stage threshold
    szd = vector(mode='list', length=length(source_zones)) # Storage for "source zone data"
    names(szd) = source_zones
    for(j in 1:length(source_zones)){
        #browser()

        source_zone = source_zones[j]
        print(source_zone)

        # Get exceedance-rate curves for the source zone (individually).
        # Using these exceedance-rate curves it should be easier to ensure
        # PTHA18 has a reasonable effective number of scenarios, compared to
        # using the sum-of-all-source-zones curve [where many source zones
        # won't really have scenarios matching rare return periods].
        stage_exceedance_rate_curves_file = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
            'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/revised1_', 
            'tsunami_stage_exceedance_rates_', source_zone, '.nc')
        fid = nc_open(stage_exceedance_rate_curves_file)
        ptha18_gauge_IDs_in_stage_exrate_file = ptha18$clean_gaugeID(ncvar_get(fid, 'gaugeID'))
        #ptha18_gauge_elev_in_stage_exrate_file = ncvar_get(fid, 'elev')
        stochastic_slip_stage_exrate_table = ncvar_get(fid, 'stochastic_slip_rate_84pc') # 84th percentile curve FOR THE SOURCE ZONE
        stochastic_slip_stages = ncvar_get(fid, 'stage') 
        nc_close(fid)

        # Many variables are held for each source zone
        # Most are stored for each site, at multiple stage values (one per
        # exceedance-rate threshold) and for each importance sample.
        szd[[source_zone]]$mean = array(NA, 
            dim=c(nrow(model_ptha18_point_gauges), 
                  length(stage_exceedance_rate_thresholds), 
                  length(all_samples))) 
            # Exact mean exrate.
        szd[[source_zone]]$var = szd[[source_zone]]$mean 
            # Exact Monte Carlo variance
        szd[[source_zone]]$best_weights = szd[[source_zone]]$mean
            # Optimal weights of each importance sample (analytical solution)
        szd[[source_zone]]$ptha18_effective_population_size = szd[[source_zone]]$mean
            # An estimate of the effective number of scenarios above the stage
            # threshold in PTHA18. This does not consider the Monte Carlo
            # sample size, it's purely about how many scenarios contribute to
            # the EXACT solution in PTHA18. 
            # It's of interest because when there is a small effective population size,
            # we should not be surprised to see erratic optimal results
        szd[[source_zone]]$IS_effective_population_size = szd[[source_zone]]$mean
            # An estimate of the effective number of scenarios above the stage
            # threshold, with scenario weights reflecting the importance
            # sampling. This does not consider the Monte Carlo sample size,
            # it's purely about how many scenarios contribute to the EXACT
            # solution when the PTHA18 weights are modified to reflect the
            # importance sampling.
            # It's of interest because when there is a small effective population size,
            # we should not be surprised to see erratic optimal results
        szd[[source_zone]]$variance_reduction = array(NA, dim=c(nrow(model_ptha18_point_gauges), length(stage_exceedance_rate_thresholds)))
            # How much is the variance improved using the optimal weights,
            # rather than just equal weighting for all sceanrios. This has one
            # less dimension than earlier variables (uses all the importance samples).
        #szd[[source_zone]]$modelled_weight_variance_reduction = szd[[source_zone]]$variance_reduction
            # As above, but using the statistical model fit to the weights, rather than the optimal weights
        #szd[[source_zone]]$modelled_weight_100_variance_reduction = szd[[source_zone]]$variance_reduction
            # As above, but using the statistical model fit to the weights (rate of 1/100), rather than the optimal weights

        szd[[source_zone]]$fitted_weights = vector(mode='list', length=length(stage_exceedance_rate_thresholds))
        names(szd[[source_zone]]$fitted_weights) = as.character(round(1/stage_exceedance_rate_thresholds))
            # Statistical Model of the weights
        # Logic tree mean rate of every scenario
        event_rate_ltm = all_samples_by_sources[[1]][[source_zone]]$event_rate_logic_tree_mean
        # Every sample is using the same logic-tree-mean rates
        for(n in 1:length(all_samples_by_sources)){
            stopifnot(all(abs(all_samples_by_sources[[n]][[source_zone]]$event_rate_logic_tree_mean - event_rate_ltm) < 1e-10))
        }

        # File containing max-stage for every event and every point on this source zone
        event_max_stage_netcdf_file = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION, 'SOURCE_ZONES/',
            source_zone, '/TSUNAMI_EVENTS/MAX_STAGE_ONLY_', 
            'all_stochastic_slip_earthquake_events_tsunami_', 
            source_zone, '_MAX_STAGE_ONLY.nc')

        # Find the point indices of interest in event_max_stage_netcdf_file
        model_ptha18_point_gauges_index = ptha18$get_netcdf_gauge_index_matching_ID(
            event_max_stage_netcdf_file, 
            ptha18$clean_gaugeID(model_ptha18_point_gauges$gaugeID))
        #if(!('elev' %in% names(model_ptha18_point_gauges))){
        #    # Append the elev to the gauges data.frame
        #    model_ptha18_point_gauges = cbind(model_ptha18_point_gauges, 
        #        data.frame(elev=ptha18_gauge_elev_in_stage_exrate_file[model_ptha18_point_gauges_index]))
        #}

        # A few sanity checks
        stopifnot(sum(is.na(model_ptha18_point_gauges_index)) == 0) # No missing matches
        stopifnot(length(model_ptha18_point_gauges_index) == length(unique(model_ptha18_point_gauges_index))) # No repeated sites
        # Gauges in event_max_stage_netcdf file are aligned with those in stage_exeedance_rate_curves_file
        stopifnot(all(ptha18_gauge_IDs_in_stage_exrate_file[model_ptha18_point_gauges_index] ==
            ptha18$clean_gaugeID(model_ptha18_point_gauges$gaugeID)))

        # Read max stage values. 
        # Read row-by-row to prevent violation of NCI's max download size
        fid = nc_open(event_max_stage_netcdf_file)
        all_max_stage = matrix(NA, nrow=length(event_rate_ltm), ncol=length(model_ptha18_point_gauges_index))
        for(ri in 1:length(model_ptha18_point_gauges_index)){
            all_max_stage[,ri] = ncvar_get(fid, 'max_stage', start=c(1, model_ptha18_point_gauges_index[ri]), count=c(fid$dim$event$len, 1))
        }
        nc_close(fid)

        # Loop over every point, and get the true exceedance-rate, and the variance
        # of the Monte Carlo exceedance-rates, for every source zone.
        for(i in 1:nrow(model_ptha18_point_gauges)){
            if((i-1)%%100 == 0) print(i)

            site_index = model_ptha18_point_gauges_index[i]
            max_stage = all_max_stage[,i] #all_max_stage[,site_index]

            # Get the stage thresholds corresponding to
            # stage_exceedance_rate_thresholds at this site. 
            # Note this is based on the "all source zones" stage exceedance-rate
            # curve, not the one for this particular source zone. That's desirable
            # because in general we are trying to approximate the 'all source
            # zones' result.
            if(sum(!is.na(stochastic_slip_stage_exrate_table[,site_index])) > 1){
                stage_thresholds = approx(stochastic_slip_stage_exrate_table[,site_index],
                    stochastic_slip_stages, xout=stage_exceedance_rate_thresholds, ties='min', rule=2)$y
            }else{
                # Sites with tiny waves that never have non-missing values in the table. Just set the stage
                # to the second lowest value of the lookup table, which will be small, and should ensure
                # zero scenarios 
                stage_thresholds = rep(stochastic_slip_stages[2], length(stage_exceedance_rate_thresholds))
            }

            for(st in 1:length(stage_thresholds)){
                for(k in 1:length(all_samples)){

                    # analytical_exrate_and_MC_variance<-function(event_rate, event_stage, stage_threshold, sampling_prob, N_MC)
                    analytical_MC_values = isu$analytical_exrate_and_MC_variance(
                        event_rate_ltm, 
                        max_stage, 
                        stage_thresholds[st], 
                        all_samples_by_sources[[k]][[source_zone]]$sampling_prob,
                        all_samples_by_sources[[k]][[source_zone]]$N_MC
                    )
                    szd[[source_zone]]$mean[i, st, k] = analytical_MC_values[1] # Note this is the rate FOR THE SOURCE ZONE, not for the sum-over-source-zones.
                    szd[[source_zone]]$var[i, st, k] = analytical_MC_values[2]

                    # Get the effective sample size for scenarios exceeding the threshold in PTHA18
                    # Equation for effective sample size, Owen (2013) Eqn 9.13
                    wt_ptha18 = event_rate_ltm*(max_stage > stage_thresholds[st])
                    wt_ptha18 = wt_ptha18/sum(wt_ptha18) # Probability of each scenario for this threshold
                    sol = sum(wt_ptha18)**2 / sum(wt_ptha18**2) # For equal weights, this would be sum(max_stage > stage_thresholds[st])
                    if(!is.finite(sol)) sol=0
                    szd[[source_zone]]$ptha18_effective_population_size[i, st,k] = sol

                    # Now get the effective sample size for scenarios exceeding the threshold in the Monte Carlo sample
                    # Equation for effective sample size, Owen (2013) Eqn 9.13
                    IS_sampling_prob = all_samples_by_sources[[k]][[source_zone]]$sampling_prob * (max_stage > stage_thresholds[st])
                    IS_sampling_prob = IS_sampling_prob / sum(IS_sampling_prob) # Probability of each scenario with max-stage > threshold
                    wt_IS = wt_ptha18 / IS_sampling_prob # Ratio in Owen (2013) Eqn 9.13
                    wt_IS[IS_sampling_prob == 0] = 0
                    sol = sum(wt_IS)**2 / sum(wt_IS**2) # For equal weights, this would be sum(max_stage > stage_thresholds[st])
                    if(!is.finite(sol)) sol=0
                    szd[[source_zone]]$IS_effective_population_size[i,st,k] = sol

                } # each importance sample 

                # Analytical weights that will minimise the sum of the variances
                vars = szd[[source_zone]]$var[i,st,]
                optimal_weights = (1/vars) * 1/(sum(1/vars))
                szd[[source_zone]]$best_weights[i,st,] = optimal_weights
                # Optimal variance divided by equal-weight variance
                # The optimal variance is sum_i{ w_i**2 \sigma_i }, where w_i are the optimal_weights (with sum_i{ w_i } == 1)
                szd[[source_zone]]$variance_reduction[i,st] = sum(optimal_weights**2 * vars)/sum((1/length(all_samples))**2 * vars)
            } # each stage threshold
        } # model_ptha18_point_gauges

        ## Fit models to the weight variation
        #for(st in 1:length(stage_exceedance_rate_thresholds)){

        #    # Might fail if the weights are NA or zero
        #    fit_wts = try(model_spatially_varying_weights(
        #        model_ptha18_point_gauges$lat, szd[[source_zone]]$best_weights[,st,], poly_degree=4))

        #    szd[[source_zone]]$fitted_weights[[st]] = fit_wts # Store the model

        #    # Compute the variance we'll get with the statistical model of the weights.
        #    vars = szd[[source_zone]]$var[,st,]
        #    # Ratio of the variance above, to what we'd get with equal weights.
        #    equal_wt_variance = rowSums( (1/3)**2 * vars)

        #    if(!is(fit_wts, 'try-error')){
        #        fit_wt_variance = rowSums(fit_wts$fitted_weights**2 * vars)
        #        szd[[source_zone]]$modelled_weight_variance_reduction[,st] = fit_wt_variance/equal_wt_variance
        #    }else{
        #        szd[[source_zone]]$modelled_weight_variance_reduction[,st] = NA
        #    }

        #    # Now repeat the above calculation, using the weight model derived at the 1/100 exceedance rate.
        #    # In practice we probably cannot vary the weights by exceedance-rate.
        #    fit_wt_100_variance = rowSums(szd[[source_zone]]$fitted_weights[[1]]$fitted_weights**2 * vars)
        #    szd[[source_zone]]$modelled_weight_100_variance_reduction[,st] = fit_wt_100_variance/equal_wt_variance

        #}

    } # each source zone

    return(environment())
}

