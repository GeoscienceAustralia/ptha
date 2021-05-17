Rprof(line.profiling=TRUE)

# Read key codes
source('R/sum_tsunami_unit_sources.R', local=TRUE)
source('R/config.R', local=TRUE)

# Make a file on NCI to use to test
# Find a file that contains hazard points. Easiest way is to read them from a tide gauge file.
unit_source_stats_puysegur = paste0(.GDATA_OPENDAP_BASE_LOCATION, 
    'SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/unit_source_statistics_puysegur2.nc')
fid = nc_open(unit_source_stats_puysegur)
gauge_netcdf_file = ncvar_get(fid, 'tide_gauge_file', start=c(1, 1), count=c(fid$dim$max_nchar$len,1))[1]
nc_close(fid)
gauge_netcdf_file = adjust_path_to_gdata_base_location(gauge_netcdf_file)

# Run the unit tests in sum_tsunami_unit_sources.R
test_sum_tsunami_unit_sources(gauge_netcdf_file)

# Make a regression test to check that we get the same results
# for Puysegur. Obviously this will have to be updated each time
# the database changes.
test_puysegur2<-function(){
    source('./get_PTHA_results.R', local=TRUE)

    # Events data [stochastic slip]. We will use this for tests here AND further down
    # in this function
    puysegur = get_source_zone_events_data('puysegur2', include_potential_energy=TRUE)

    # Check the potential energy by comparison with a few independently calculated values
    # (which will not be exactly the same because of discretization/projection/DEM differences
    # , but are very similar)
    inds = c(1643    , 1767    ) 
    pes  = c(6.93e+12, 8.86e+12)
    if(all(abs(puysegur$events$initial_potential_energy[inds] - pes) < 0.01*pes)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # As above, but check uniform slip with some subset of rows
    inds = c(107     , 117     )
    pes  = c(2.62e+12, 1.13e+13)
    puysegur_U = get_source_zone_events_data('puysegur2', slip_type='uniform', 
        desired_event_rows = inds, include_potential_energy=TRUE) 
    if(all(abs(puysegur_U$events$initial_potential_energy - pes) < 0.01*pes)){
        print('PASS')
    }else{
        print('FAIL')
    }
    rm(inds, pes, puysegur_U)

    # As above, but check variable_uniform slip with some subset of rows
    inds = c(1815    , 1835    )
    pes  = c(2.77e+13, 1.39e+13)
    puysegur_VU = get_source_zone_events_data('puysegur2', slip_type='variable_uniform', 
        desired_event_rows = inds, include_potential_energy=TRUE) 
    if(all(abs(puysegur_VU$events$initial_potential_energy - pes) < 0.01*pes)){
        print('PASS')
    }else{
        print('FAIL')
    }
    rm(inds, pes, puysegur_VU)

    model_3051 = get_flow_time_series_at_hazard_point(puysegur, 3051, 
        hazard_point_ID = c(1.1, 10.1, 22.1, 55015.4, 55042.4))

    max_55015 = max(model_3051$flow[['55015.4']][1,,1])
    max_55042 = max(model_3051$flow[['55042.4']][1,,1])

    # Values known from previous checks
    er1 = abs(max_55015 - 0.1771731) #0.050279)
    er2 = abs(max_55042 - 0.3976239) #0.138759)

    m1 = which.max(model_3051$flow[['55015.4']][1,,1])
    m2 = which.max(model_3051$flow[['55042.4']][1,,1])

    if((er1 < 1.0e-05) & (er2 < 1.0e-05) & (m1 == 59) & (m2 == 61)){
        print('PASS')
    }else{
        print('FAIL: test_puysegur giving different results')

    }


    # Check location info is right -- this is based on 'dput' of an earlier result
    test_data = structure(list(lon = structure(c(129.089294433594, 129.18830871582, 
        127.597557067871, 160.256103515625, 161.841110229492), .Dim = 5L), 
        lat = structure(c(-8.40988159179688, -8.29211044311523, -8.10336685180664, 
        -46.8297233581543, -44.897777557373), .Dim = 5L), elev = structure(c(-775.006286621094, 
        -716.009460449219, -455.951629638672, -5052.998046875, -4819.00732421875
        ), .Dim = 5L), gaugeID = structure(c(1.1, 10.1, 22.1, 55015.4, 
        55042.4), .Dim = 5L)), class = "data.frame", row.names = c(NA, 
        -5L))

    m1 = max(abs(model_3051$location - test_data))
    if(m1 > 1.0e-010){
        print('FAIL -- location information haz changed')
    }else{
        print('PASS')
    }

    # Check if we unpack to gauge-based list, it still works
    model_3051b = get_flow_time_series_at_hazard_point(puysegur, 3051, 
        hazard_point_ID = c(1.1, 10.1, 22.1, 55015.4, 55042.4),
        store_by_gauge=FALSE)


    if(length(model_3051b$flow) == dim(model_3051$flow[[1]])[1]){
        print('PASS')
    }else{
        print('FAIL: store_by_gauge not working as desired')
    }

    FAILED=FALSE
    for(i in 1:length(model_3051b$flow)){
        for(j in 1:length(model_3051$flow)){
            if(! all(model_3051b$flow[[i]][j,,] == model_3051$flow[[j]][i,,])){
                FAILED=TRUE
            }
        }
    }
    if(FAILED){
        print('FAIL: store_by_gauge not ordered as desired')
    }else{
        print('PASS')
    }
 
}
# Run the puysegur regression test
t1 = system.time(test_puysegur2())


# Test that we can read a scenario with many unit sources -- i.e. one that
# would fail if the netcdf install is not sufficiently up to date.

test_large_event_read<-function(){

    source('./get_PTHA_results.R', local=TRUE)

    expected_event_index_string = "324-325-329-330-333-334-338-344-345-348-349-350-351-352-354-355-358-359-361-362-363-365-366-367-369-370-372-373-376-377-381-384-385-386-387-388-389-390-391-392-393-394-395-396-397-398-399-400-401-402-403-404-405-406-407-408-409-410-411-412-413-414-415-416-417-418-419-420-421-422-423-424-425-426-427-428-430-431-434-435-436-438-439-440-443-454-461-462-"

    expected_event_slip_string = "2.428_0.4631_6.739_1.026_6.742_8.005_7.856_2.007_2.8_7.307_1.599_19.78_15.27_4.94_21.14_29.31_25.13_18.47_10.37_19.73_20.71_26.1_4.098_3.605_35.83_6.34_6.517_34.34_3.984_30.11_41.04_0.7955_57.45_32.88_25.11_35.7_77.66_84.67_64.68_56.23_87.42_117_91.92_68.23_114.9_161.2_142.1_109.2_162.3_196.5_198.3_153.7_186_196.9_206.4_182.5_169.6_163_187.9_181.4_127.7_133.5_155.5_144_85.38_102.2_105.1_73.65_30.45_56.9_59.6_36.4_8.488_38.66_80.98_16.74_46.75_70.94_51.39_50.19_24.28_20.31_24.2_14.35_7.434_1.672_2.327_2.264_"

    x = get_source_zone_events_data('southamerica', slip_type='stochastic', desired_event_rows=150000)

    test_result = (x$events$event_slip_string == expected_event_slip_string &
        x$events$event_index_string == expected_event_index_string)

    if(test_result){
        print('PASS')
    }else{
        print('FAIL -- incorrect read of a long string of event metadata. This suggests your netcdf install is not sufficiently up to date (there was a bug in the remote reading of long character strings in netcdf-c, whcih was fixed in the bleeding edge source in late 2017)')
    }

}

test_large_event_read()


# Test that we can work with the detailed scenario information. ONLY do it if
# compute_rates_all_sources_session.RData exists locally [since otherwise the
# download will be slow]
if(file.exists('compute_rates_all_sources_session.RData')){
    detailed_ptha18 = new.env()
    source('get_detailed_PTHA18_source_zone_info.R', local=detailed_ptha18)
    detailed_ptha18$test_get_detailed_PTHA18_source_zone_info()
}


test_random_scenario_sampling<-function(){
    #
    # Check various functions for random scenario sampling (including importance sampling)
    #

    ptha18 = new.env()
    source('./get_PTHA_results.R', local=ptha18, chdir=TRUE)

    #
    # Get scenario data for the testing
    #

    # Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in PTHA18)
    source_zone = 'kermadectonga2'
    kt2_scenarios = ptha18$get_source_zone_events_data(source_zone,  slip_type='stochastic')

    # Convenient shorthand for the magnitudes and rates in the event table
    event_Mw = kt2_scenarios$events$Mw 
    event_rates = kt2_scenarios$events$rate_annual
    event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
        target_point = c(185.1239, -21.0888), # Known location of PTHA18 hazard point
        slip_type='stochastic',
        all_source_names=source_zone)
    # Convenient shorthand
    event_peak_stage = event_peak_stage_at_refpoint$kermadectonga2$max_stage


    # Make the test reproducible despite our random scenarios
    set.seed(123)

    random_scenarios = vector(mode='list', length=3)
    names(random_scenarios) = c('simple', 'mw_stratified', 'stage_mw_weighted')
    # Sample random scenarios
    # -- simple random sampling
    random_scenarios[['simple']] = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
        event_rates=event_rates,
        event_Mw=event_Mw,
        samples_per_Mw=function(Mw){ 12 }, 
        mw_limits=c(7.15, 9.85)
        )
    # -- stratified sampling, with more samples at large magnitudes
    random_scenarios[['mw_stratified']] = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
        event_rates=event_rates,
        event_Mw=event_Mw,
        samples_per_Mw=function(Mw){ round( 6 + 12 * (Mw - 7.15)/(9.65 - 7.15) ) },
        mw_limits=c(7.15, 9.85)
        )
    # -- combined importance-sampling + stratified sampling, with more samples at large magnitudes
    random_scenarios[['stage_mw_weighted']] = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
        event_rates=event_rates,
        event_Mw=event_Mw,
        event_importance_weighted_sampling_probs = (event_peak_stage*event_rates),
        samples_per_Mw=function(Mw){ round( 6 + 12 * (Mw - 7.15)/(9.65 - 7.15) ) },
        mw_limits=c(7.15, 9.85)
        )

    #
    # Test of mean/variance estimators of exceedance-rates, for different importance-sampling types
    #
    importance_sampling_types = c('basic', 'self_normalised', 'control_variate')

    # Make empty data-structure with results for each sampling-type and importance-sampling-type
    exrate_uncertainty = vector(mode='list', length=length(random_scenarios))
    names(exrate_uncertainty) = names(random_scenarios)
    for(rsn in names(random_scenarios)){
        exrate_uncertainty[[rsn]] = vector(mode='list', length=length(importance_sampling_types))
        names(exrate_uncertainty[[rsn]]) = importance_sampling_types
    }

    threshold_stage = 2.3
    # Make estimates of the mean/variance of threshold_stage exceedances for
    # each Mw-bin, for each sampling-type and importance-sampling-type
    for(rsn in names(random_scenarios)){
        for(ist in importance_sampling_types){

            exrate_uncertainty[[rsn]][[ist]] = ptha18$get_exrate_uncertainty_at_stage(
                random_scenarios[[rsn]], event_peak_stage, threshold_stage, 
                importance_sampling_type=ist, return_per_Mw_bin=TRUE)
        }
    }


    # For simple random sampling, and mw-stratified-sampling, the 3
    # importance_sampling_types should lead to identical results.
    for(sampling_type in c('simple', 'mw_stratified')){
        for(ist in importance_sampling_types[2:3]){
            the_test = isTRUE(all.equal(
                exrate_uncertainty[[sampling_type]][[ist]], 
                exrate_uncertainty[[sampling_type]][['basic']]))
            if(the_test){
                print('PASS')
            }else{
                print(paste0('FAIL -- all exrate estimators should be the same for ', sampling_type, 
                             ' random sampling'))
            }
        }
    }
    
    # For simple random sampling, and mw-stratified random sampling, the
    # mean/variance of each technique should be equal to the typical estimates 
    for(sampling_type in c('simple', 'mw_stratified')){

        # Prepare to compute the mean/variance from standard theory

        rate_with_this_mw = aggregate(random_scenarios[[sampling_type]]$rate_with_this_mw, 
            by=list(random_scenarios[[sampling_type]]$mw), function(x) x[1])
        inv_count_with_this_mw = aggregate(random_scenarios[[sampling_type]]$rate_with_this_mw, 
            by=list(random_scenarios[[sampling_type]]$mw), function(x) 1/length(x))
        pb = aggregate(event_peak_stage[random_scenarios[[sampling_type]]$inds] > threshold_stage, 
            by=list(random_scenarios[[sampling_type]]$mw), mean)

        # Naive-estimate of mean/variance
        typical_rate_estimate = rate_with_this_mw$x * pb$x 
        typical_rate_variance_estimate = rate_with_this_mw$x^2 * pb$x * (1-pb$x) * inv_count_with_this_mw$x

        # The trailing NA's can be replaced with 0 (they occur because mw 9.7, 9.8 are impossible)
        n = length(typical_rate_estimate)
        typical_rate_estimate[(n-1):n] = 0
        typical_rate_variance_estimate[(n-1):n] = 0
        
        # Make a data.frame matching the data-structure in exrate_uncertainty[[sampling_type]][[ist]]
        typical_estimates = exrate_uncertainty[[sampling_type]][[1]]
        typical_estimates$exrate = typical_rate_estimate
        typical_estimates$exrate_variance = typical_rate_variance_estimate

        for(ist in importance_sampling_types[1:3]){
            if(isTRUE(all.equal(exrate_uncertainty[[sampling_type]][[ist]], typical_estimates))){
                print('PASS')
            }else{
                print('FAIL -- in this situation all exrate estimators should match the typical estimate')
            }
        }
    }

    #
    # Test the coverage of approximate confidence intervals derived from the
    # importance-sampling variance estimators. Unlike the tests above, here we
    # are considering 'real' importance sampling
    #
    true_exrate = sum(event_rates * (event_peak_stage > threshold_stage))
    Nsam = 1000

    # Store the coverage of confidence intervals
    does_interval_cover_normal = vector(mode='list', length=3)
    names(does_interval_cover_normal) = importance_sampling_types
    for(i in 1:3) does_interval_cover_normal[[i]] = rep(NA, Nsam)
    does_interval_cover_beta = does_interval_cover_normal

    # Store the parameter estimates for later analysis
    est_sd_store = vector(mode='list', length=3)
    names(est_sd_store) = importance_sampling_types
    for(i in 1:3) est_sd_store[[i]] = matrix(NA, ncol=2, nrow=Nsam)

    for(i in 1:Nsam){

        random_scenarios_repeated = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
            event_rates=event_rates,
            event_Mw=event_Mw,
            event_importance_weighted_sampling_probs = (event_peak_stage*event_rates),
            samples_per_Mw=function(Mw){ round( 6 + 12 * (Mw - 7.15)/(9.65 - 7.15) ) },
            mw_limits=c(7.15, 9.85)
            )

        # Get approximate confidence intervals for each sampling type
        for(is_type in importance_sampling_types){

            exrate_uncertainty_repeated = ptha18$get_exrate_uncertainty_at_stage(
                random_scenarios_repeated, event_peak_stage, threshold_stage, 
                importance_sampling_type=is_type)

            est = exrate_uncertainty_repeated[1]
            sd_est = sqrt(exrate_uncertainty_repeated[2])

            est_sd_store[[is_type]][i,1] = est
            est_sd_store[[is_type]][i,2] = sd_est

            # 95% confidence interval (normal approximation).
            # (-1 -- under-estimate)
            # ( 0 --> covered )
            # (+1 --> over-estimate)
            normal_limits = qnorm(c(0.025, 0.975), mean=est, sd=sd_est)
            does_interval_cover_normal[[is_type]][i] = # 95% 2-sided interval (2.5% above, 2.5% below)
                -1*( true_exrate < normal_limits[1]) +
                 1*( true_exrate > normal_limits[2])

            # Alternative confidence interval (beta approximation, scaled to span from 
            #    [0, rate_of_any_event] )

            # Get the 'true' rate of any event (use the self-normalised
            # rates, because they are forced to match the true rate)
            rate_of_any_event = sum(
                random_scenarios_repeated$importance_sampling_scenario_rates_self_normalised, 
                na.rm=TRUE)
                
            beta_params_from_mu_var <- function(mu, var) {
                alpha = ((1 - mu) / var - 1 / mu) * mu ^ 2
                beta = alpha * (1 / mu - 1)
                params = list(alpha = alpha, beta = beta)
                return(params)
                }
    
            beta_par = beta_params_from_mu_var(est/rate_of_any_event, (sd_est/rate_of_any_event)^2)
            # 95% 2-sided confidence-interval
            beta_limits =  rate_of_any_event * qbeta(c(0.025, 0.975), shape1=beta_par$alpha, 
                shape2=beta_par$beta)

            ## 95% confidence interval coverage
            # (-1 -- under-estimate)
            # ( 0 --> covered )
            # (+1 --> over-estimate)
            does_interval_cover_beta[[is_type]][i] = 
                -1*( true_exrate < beta_limits[1]) +
                 1*( true_exrate > beta_limits[2])

        }
    }

    # Tests of the estimates
    error_thresholds = list()
    # Store the thresholds as tolerances for 
    #    c(median_relative_error, 
    #      sd(relative_error), 
    #      empirical_coverage_of_normal_95%_interval, 
    #      empirical_coverage_of_beta_95%_interval)
    # The thresholds for 'basic' importance sampling are more stringent because
    # that method seems to work better.
    error_thresholds$basic           = c(1.0e-02, 0.15, 0.9, 0.9) 
    error_thresholds$self_normalised = c(0.15   , 0.25, 0.8, 0.8)
    error_thresholds$control_variate = c(0.15   , 0.25, 0.8, 0.8)

    for(is_type in importance_sampling_types){

        # Bias of the estimator
        is_errors = (est_sd_store[[is_type]][,1] - true_exrate)/true_exrate
        if( (abs(median(is_errors)) < error_thresholds[[is_type]][1]) & 
            (sd(is_errors)          < error_thresholds[[is_type]][2])){
            print('PASS')
        }else{
            print(paste0('FAIL -- ', is_type, ' importance sampling estimate'))
        }

        # 95% confidence interval coverage
        if((mean(does_interval_cover_normal[[is_type]] == 0) > error_thresholds[[is_type]][3]) &
           (mean(does_interval_cover_beta[[is_type]]   == 0) > error_thresholds[[is_type]][4])){
            print('PASS')
        }else{
            print(paste0('FAIL -- confidence interval for ', is_type, ' importance sampling'))
        }
    }


    #
    # Test our estimates of the optimal number of scenarios to sample in each Mw bin
    #

    # Simple check of variance_numerator by back-calculation
    t0 = ptha18$get_optimal_number_of_samples_per_Mw(event_Mw, event_rates, event_peak_stage,
        stage_threshold=threshold_stage, total_samples=nrow(random_scenarios_repeated))
    rate_with_this_Mw = aggregate(event_rates, by=list(event_Mw), sum)$x
    rate_exceeding_with_this_Mw = aggregate(event_rates * (event_peak_stage > threshold_stage), 
        by=list(event_Mw), sum)$x
    pb = rate_exceeding_with_this_Mw/rate_with_this_Mw
    pb[rate_with_this_Mw == 0] = 0 # Fix NaN values
    variance_numerator = rate_with_this_Mw^2 * pb * (1-pb)

    if(isTRUE(all.equal(variance_numerator, t0$variance_numerator, tol=1.0e-14))){
        print('PASS')
    }else{
        print('FAIL -- variance_numerator is not as expected')
    }

    #
    # Check that the theoretical variance (assuming we sample per Mw-bin like
    # in random_scenarios_repeated) is close to the empirical variances (from
    # the above repeated random sampling). This exercises importance-sampling as well
    # as unequal magnitude stratification.
    #
    t0 = ptha18$get_optimal_number_of_samples_per_Mw(event_Mw, event_rates, event_peak_stage,
        stage_threshold=threshold_stage, total_samples=nrow(random_scenarios_repeated),
        event_importance_weighted_sampling_probs=(event_peak_stage*event_rates))
    repeated_random_scenario_counts = as.numeric(table(random_scenarios_repeated$mw))
    expected_variance = sum(t0$variance_numerator/repeated_random_scenario_counts)
    empirical_variance = var(est_sd_store[['basic']][,1])
    # Check they agree within a few percent
    if(abs(empirical_variance - expected_variance) < 0.03*expected_variance){
        print('PASS')
    }else{
        print('FAIL -- the theoretical variance is not close enough to the sampling variance')
    }

    # Check that the typical 'sample-sd-estimate' is reasonably close to the expected value
    typical_sd_estimate = median(est_sd_store[['basic']][,2])
    expected_sd = sqrt(expected_variance)
    if( abs(typical_sd_estimate - expected_sd) < 0.05*expected_sd ){
        print('PASS')
    }else{
        print('FAIL -- the typical estimate of the sd is not close enough to the theoretical sd')
    }

    # Check the optimality of the variance estimate in a range of cases. For
    # this check, we do not round the sample size (although that is necessary
    # in practice, we lose the guarentee of optimality, and the latter is
    # convenient for the testing.

    t1 = ptha18$get_optimal_number_of_samples_per_Mw(event_Mw, event_rates, event_peak_stage,
        stage_threshold=1, total_samples=1200)
    t2 = ptha18$get_optimal_number_of_samples_per_Mw(event_Mw, event_rates, event_peak_stage,
        stage_threshold=2, total_samples=1200)
    t5 = ptha18$get_optimal_number_of_samples_per_Mw(event_Mw, event_rates, event_peak_stage,
        stage_threshold=5, total_samples=1200)

    # A good compromise sampling effort
    test_N = (t1$Nsamples + t2$Nsamples + t5$Nsamples)/3

    a1 = sqrt( sum(t1$variance_numerator/t1$Nsamples, na.rm=TRUE) )
    #[1] 0.000280247
    a2 = sqrt( sum(t1$variance_numerator/(1200/nrow(t1)), na.rm=TRUE) ) # Equal sampling
    #[1] 0.0003843171
    a3 = sqrt( sum(t1$variance_numerator/(test_N), na.rm=TRUE) ) # Alternative
    if(a1 < 0.75*a2 & a1 < a3){
        print('PASS')
    }else{
        print('FAIL -- number of samples is not optimal at stage_threshold=1')
    }

    a1 = sqrt( sum(t2$variance_numerator/t2$Nsamples, na.rm=TRUE) )
    # [1] 0.0001268824
    a2 = sqrt( sum(t2$variance_numerator/(1200/nrow(t1)), na.rm=TRUE) )
    # [1] 0.0001814382
    a3 = sqrt( sum(t2$variance_numerator/(test_N), na.rm=TRUE) ) # Alternative
    if(a1 < 0.7*a2 & a1 < a3){
        print('PASS')
    }else{
        print('FAIL -- number of samples is not optimal at stage_threshold=2')
    }

    a1 = sqrt( sum(t5$variance_numerator/t5$Nsamples, na.rm=TRUE) )
    # [1] 2.825913e-05
    a2 = sqrt( sum(t5$variance_numerator/(1200/nrow(t5)), na.rm=TRUE) )
    # [1] 4.999452e-05
    a3 = sqrt( sum(t5$variance_numerator/(test_N), na.rm=TRUE) ) # Alternative
    if(a1 < 0.6*a2 & a1 < a3){
        print('PASS')
    }else{
        print('FAIL -- number of samples is not optimal at stage_threshold=5')
    }

}
test_random_scenario_sampling()


Rprof(NULL)
