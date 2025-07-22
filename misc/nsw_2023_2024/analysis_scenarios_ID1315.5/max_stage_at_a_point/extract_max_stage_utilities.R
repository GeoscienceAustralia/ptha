library(stars)

# Importance sampling utilities
isu = new.env()
source('../../sources/hazard/importance_sampling_utilities.R', local=isu)

# Get the indices of the target_point in our raster file
get_indices_of_target_point_on_raster_file<-function(tarred_raster_file, raster_filename, target_point){
    # Get the indices of interest by checking one tif
    # Here I assume all tifs with the same raster_filename have the same extent (but that will be checked)
    r1 = read_stars(paste0('/vsitar/', tarred_raster_file, '/', raster_filename))

    # Get XIND,YIND such that the value of interest is r1[[1]][XIND, YIND]
    # Beware: If the raster had been read with terra then these indices are reversed
    #   i.e. r1_stars[[1]][XIND, YIND] = r1_terra[YIND, XIND]
    XIND = which.min(abs( st_get_dimension_values(r1, 'x', where='center') - target_point[1]))
    YIND = which.min(abs( st_get_dimension_values(r1, 'y', where='center') - target_point[2]))
    XCOORD = st_get_dimension_values(r1, 'x', where='center')[XIND]
    YCOORD = st_get_dimension_values(r1, 'y', where='center')[YIND]

    r1_dimensions = st_dimensions(r1)
    # Gauge should be within dx/2 of the cell centre location.
    r1_errtol = abs(c(r1_dimensions$x$delta, r1_dimensions$y$delta))/2
    rm(r1); gc()

    if(any(abs(c(XCOORD, YCOORD) - target_point) > r1_errtol)) stop('Coordinate error')

    return(list(XIND=XIND, YIND=YIND, coordinate_errtol = r1_errtol ))
}

# Get the max-stage at the desired point from a single file.
get_max_stage_at_target_point<-function(tarred_raster_file, raster_filename, XIND, YIND, target_point, coordinate_errtol){

    # Read the file
    raster_file = paste0('/vsitar/', tarred_raster_file, '/', raster_filename)
    r1 = read_stars(raster_file)
    coordinate_vals = c(
        st_get_dimension_values(r1, 'x', where='center')[XIND],
        st_get_dimension_values(r1, 'y', where='center')[YIND])

    segfault_workaround=TRUE # Avoid some segfaults with stars
    if(!segfault_workaround){
        # This was throwing segfaults on NCI on 19/02/2024 for NSW run. 
        # Didn't see this issue previously (many WA model runs).
        val = r1[[1]][XIND, YIND]
    }else{
        # Use terra to avoid segfaults
        library(terra)
        r2 = rast(raster_file)
        val = extract(r2, matrix(coordinate_vals, ncol=2, nrow=1))[1,1]
        rm(r2)
    }

    rm(r1); gc()

    # Check the coordinates are as expected
    if(any(abs(coordinate_vals - target_point) > coordinate_errtol)){
        return(NA)
    }else{
        # Return the value of interest
        return(val)
    }

}

plot_tsunami_maxima_in_nonlinear_model_and_PTHA18<-function(target_point, scenario_max_stages){
    # Compare nonlinear model and PTHA18 (36 hours, frictionless). They should generally agree
    # at deep water points. Typically PTHA18 tends to be larger due to lack of friction and 
    # longer run-time, but differences are smaller very offshore. 
    #
    sourceID = as.integer(as.factor(scenario_max_stages$nonlinear_source_name))
    COLZ = colorRampPalette(rainbow(max(sourceID)+2))(max(sourceID)+2) # For plot

    # Get a mapping between the source ID and the source name, to make a plot legend
    source_name_ID = aggregate(sourceID, 
        by=list(source_name=scenario_max_stages$nonlinear_source_name), 
        function(x) x[1])

    out_file = paste0('Tsunami_maxima_in_PTHA18_and_Nonlinear_Model_', target_point[1], '_', target_point[2], '.png')
    png(out_file, width=12, height=9, units='in', res=300)
    plot(scenario_max_stages$ptha18_max_stage, scenario_max_stages$nonlinear_max_stage_minus_MSL,
         xlab='PTHA18 Max-stage (36hours, linear frictionless model)',
         ylab='Highres model max-stage',
         cex.lab=1.5, cex.axis=1.5, pch=19, asp=1,
         col=COLZ[sourceID])
    grid(col='orange')
    abline(0, 1, col='red')
    abline(0, median(scenario_max_stages$nonlinear_max_stage_minus_MSL/scenario_max_stages$ptha18_max_stage), col='blue')
    title(main='Tsunami maxima in PTHA18 vs High-resolution model',
        cex.main=1.7)
    legend('bottomright', source_name_ID$source_name, pch=19, col=COLZ[source_name_ID$x], bty='n') 
    dev.off()

    # Another way of displaying this
    out_file = paste0('Tsunami_maxima_in_PTHA18_and_Nonlinear_Model_relative_', target_point[1], '_', target_point[2], '.png')
    mean_max_stage = scenario_max_stages$ptha18_max_stage 
    diff_max_stage = (scenario_max_stages$ptha18_max_stage - scenario_max_stages$nonlinear_max_stage_minus_MSL)
    reldiff =  diff_max_stage/mean_max_stage
    png(out_file, width=12, height=9, units='in', res=300)
    scatter.smooth(mean_max_stage, reldiff, xlab='PTHA18 max-stage',
         ylab='[PTHA18 - Nonlinear]/PTHA18',
         col=COLZ[sourceID],
        cex.lab=1.5, cex.axis=1.5, pch=19)
    abline(h=0, col='red', lty='dashed')
    title(main='PTHA18 max-stage vs Relative difference in nonlinear model', cex.main=1.7)
    legend('bottomright', source_name_ID$source_name, pch=19, col=COLZ[source_name_ID$x], bty='n') 
    dev.off()
}

# Plot logic-tree-mean exceedance-rate curves from PTHA18 and nonlinear model,
# separately for each source-zone, and in combination. The differences reflect
# both Monte Carlo sampling and differences in the hydrodynamic models. We
# also make plots where the nonlinear-model results are derived using ptha18
# max-stages, which is useful for testing since in this situation the
# differences only reflect Monte Carlo sampling.
#
plot_logic_tree_mean_exrate_curves_from_nonlinear_model_and_PTHA18<-function(
    target_point, 
    nonlinear_model_curves, nonlinear_model_combined_curve, nonlinear_model_MSL, 
    ptha18_curves, ptha18_combined_curve, EXRATE_PLOT_YLIM){

    # Plot the nonlinear model exceedance rates for each source zone.
    out_file = paste0('Exceedance_Rates_in_PTHA18_and_Nonlinear_Model_', 
        target_point[1], '_', target_point[2], '.png')
    NPANEL = length(names(nonlinear_model_curves))
    HT = 7.5
    WD = HT * NPANEL
    png(out_file, width=WD, height=HT, units='in', res=300)
    par(mfrow=c(1,NPANEL))
    options(scipen=5)
    for(nm in names(nonlinear_model_curves)){
        plot(ptha18_curves[[nm]], t='l', log='xy', ylim=EXRATE_PLOT_YLIM, col='skyblue', lwd=3,
            xlab='Max-stage above MSL (m)', ylab='Exceedance-rate (events/year)', cex.lab=1.4, cex.axis=1.4)
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               nonlinear_model_curves[[nm]]$exrate, col='black', t='l', lwd=2)
        title(paste0('Exceedance-rate curves: ', nm), cex.main=1.5)

        # Importance-sampling CIs
        lower_CI = nonlinear_model_curves[[nm]]$exrate + qnorm(0.025)*sqrt(nonlinear_model_curves[[nm]]$exrate_variance)
        upper_CI = nonlinear_model_curves[[nm]]$exrate + qnorm(0.975)*sqrt(nonlinear_model_curves[[nm]]$exrate_variance)
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               pmax(lower_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               pmax(upper_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot

        grid(col='orange')

        legend('topright', c('PTHA18', 'High-res nonlinear model', 'High-res 95% CI'),
               col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
               pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
    }
    dev.off()

    # Make a plot like the above, but replace the nonlinear model max-stages
    # with the PTHA18 max-stages (adjusted for MSL). This is useful for
    # testing, because with this approach, all differences in results are due
    # to Monte Carlo sampling.
    out_file = paste0('Exceedance_Rates_in_PTHA18_and_Sampled_Scenarios_using_PTHA18_max_stages_', 
        target_point[1], '_', target_point[2], '.png')
    NPANEL = length(names(nonlinear_model_curves))
    HT = 7.5
    WD = HT * NPANEL
    png(out_file, width=WD, height=HT, units='in', res=300)
    par(mfrow=c(1,NPANEL))
    options(scipen=5)
    for(nm in names(nonlinear_model_curves)){
        plot(ptha18_curves[[nm]], t='l', log='xy', ylim=EXRATE_PLOT_YLIM, col='skyblue', lwd=3,
            xlab='Max-stage above MSL (m)', ylab='Exceedance-rate (events/year)', cex.lab=1.4, cex.axis=1.4)
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               nonlinear_model_curves[[nm]]$exrate_using_ptha18_stages_plus_nonlinear_model_MSL, col='black', t='l', lwd=2)
        title(paste0('Exceedance-rate curves: ', nm), cex.main=1.5)

        # Importance-sampling CIs
        lower_CI = nonlinear_model_curves[[nm]]$exrate_using_ptha18_stages_plus_nonlinear_model_MSL + 
            qnorm(0.025)*sqrt(nonlinear_model_curves[[nm]]$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL)
        upper_CI = nonlinear_model_curves[[nm]]$exrate_using_ptha18_stages_plus_nonlinear_model_MSL + 
            qnorm(0.975)*sqrt(nonlinear_model_curves[[nm]]$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL)
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               pmax(lower_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot
        points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
               pmax(upper_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot

        grid(col='orange')

        legend('topright', c('PTHA18', 'sampled PTHA18', 'sampled PTHA18 95% CI'),
               col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
               pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
    }
    dev.off()

    # Plot the nonlinear model exceedance rates, summed over all source zones.
    out_file = paste0('Exceedance_Rates_Summed_in_PTHA18_and_Nonlinear_Model_',
        target_point[1], '_', target_point[2], '.png')
    png(out_file, width=HT, height=HT, units='in', res=300)
    options(scipen=5)
    plot(ptha18_combined_curve, t='l', log='xy', ylim=EXRATE_PLOT_YLIM, col='skyblue', lwd=3,
        xlab='Max-stage above MSL (m)', ylab='Exceedance-rate (events/year)', cex.lab=1.4, cex.axis=1.4)
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           nonlinear_model_combined_curve$exrate, col='black', t='l', lwd=2)
    title(paste0('Exceedance-rate curves'), cex.main=1.5)
    # Importance-sampling CIs
    lower_CI = nonlinear_model_combined_curve$exrate + qnorm(0.025)*sqrt(nonlinear_model_combined_curve$exrate_variance)
    upper_CI = nonlinear_model_combined_curve$exrate + qnorm(0.975)*sqrt(nonlinear_model_combined_curve$exrate_variance)
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           pmax(lower_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           pmax(upper_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot

    grid(col='orange')

    legend('topright', c('PTHA18', 'High-res nonlinear model', 'High-res 95% CI'),
           col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
           pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
    dev.off()


    # Plot like above, but replacing the nonlinear-model max-stages with the
    # PTHA18 max-stages (adjusted for MSL). This is useful for testing, because
    # with this approach, all differences in results are due to Monte Carlo
    out_file = paste0('Exceedance_Rates_Summed_in_PTHA18_and_Sampled_scenarios_using_PTHA18_max_stages_',
        target_point[1], '_', target_point[2], '.png')
    png(out_file, width=HT, height=HT, units='in', res=300)
    options(scipen=5)
    plot(ptha18_combined_curve, t='l', log='xy', ylim=EXRATE_PLOT_YLIM, col='skyblue', lwd=3,
        xlab='Max-stage above MSL (m)', ylab='Exceedance-rate (events/year)', cex.lab=1.4, cex.axis=1.4)
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL, col='black', t='l', lwd=2)
    title(paste0('Exceedance-rate curves'), cex.main=1.5)
    # Importance-sampling CIs
    lower_CI = nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL + 
        qnorm(0.025)*sqrt(nonlinear_model_combined_curve$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL)
    upper_CI = nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL + 
        qnorm(0.975)*sqrt(nonlinear_model_combined_curve$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL)
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           pmax(lower_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot
    points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
           pmax(upper_CI, 1e-200), col='black', t='l', lty='dashed') # Lower bound for log-log plot

    grid(col='orange')

    legend('topright', c('PTHA18', 'sampled PTHA18', 'sampled PTHA18 95% CI'),
           col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
           pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
    dev.off()

}

# Compute stage exceedance-rate curves using importance sampling for all
# logic-tree-branches, ASSUMING heterogeneous slip scenarios with constant rigidity.
# 
# The function output lists have entries similar enough to those of a routine
# for stratified-importance sampling
#   get_detailed_PTHA18_source_zone_info.R::random_scenario_exceedance_rates_all_logic_tree_branches
# that we can run the code (previously used for stratified-importance-sampling)
#   get_PTHA_results.R::compute_exceedance_rate_percentiles_with_random_sampling
# to compute percentile uncertainties.
#
estimate_stage_exrate_curves_with_IS_for_all_logic_tree_branches<-function(
    ptha18_detailed,
    source_zone,
    source_zone_segment,
    importance_sampling_scenarios_logic_tree_mean_on_source_zone,
    sampled_scenarios_max_stage,
    threshold_stage_values){

    # Chance of sampling each sampled scenario. This was determined at the time of sampling
    sampling_prob_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$sampling_prob
    # Indices of sampled scenarios in the PTHA18 event tables - may have repeated values
    row_ind_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$inds 
    # Scenario magnitude
    mw_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean_on_source_zone$event_Mw

    # Logic-tree-mean rates of the sampled scenarios (weighted unseg+seg).
    # We won't use this directly but it might be useful for checking
    # rates_sampled_scenarios_LTM = importance_sampling_scenarios_logic_tree_mean_on_source_zone$event_rate_logic_tree_mean


    # Get information on all rate models for this source zone and segment (or unsegmented model) unsegmented model. 
    is_segment = (source_zone_segment != "")
    if(is_segment){
        source_seg_name = paste0(source_zone, '_', source_zone_segment)
    }else{
        source_seg_name = source_zone
    }
    stopifnot(source_seg_name %in% names(ptha18_detailed$crs_data$source_envs))

    # Get rate models with constant rigidity for all logic tree branches
    # (There are different posterior weights for models with variable rigidity,
    # here we will only use constant rigidity).
    rate_models_all_branches = ptha18_detailed$crs_data$source_envs[[source_seg_name]]$mw_rate_function(
            NA, return_all_logic_tree_branches=TRUE)

    # Get the PTHA18 scenario conditional probability (for ALL SCENARIOS)
    # according to the logic-tree-mean model over the segment. This also
    # includes information on PTHA18 LTM rates, but we don't need to use that.
    ptha18_conditional_probability_and_rates = 
        ptha18_detailed$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
            source_zone=source_zone, segment=source_zone_segment)

    # The scenario conditional probability does not change between
    # logic-tree-branches on the source_zone and segment.
    # NOTE WE ASSUME CONSTANT RIGIDITY HS SCENARIOS
    ptha18_conditional_prob_given_mw_sampled_scenarios = 
        ptha18_conditional_probability_and_rates$HS_prob_given_Mw[row_ind_sampled_scenarios]

    # Range of Mw values in PTHA18 = 7.2, 7.3, 7.4, .... 9.6, 9.7, 9.8
    ptha18_unique_mw_values = sort(unique(ptha18_conditional_probability_and_rates$HS_mw))
    stopifnot(isTRUE(all(abs(diff(ptha18_unique_mw_values) - 0.1) < 1e-10)))
    # Mw bin boundaries in PTHA18 = 7.15, 7.25, ... 9.65, 9.75, 9.85
    ptha18_unique_mw_bb = c(ptha18_unique_mw_values - 0.05, max(ptha18_unique_mw_values) + 0.05)

    # Mw bin corresponding to sampled scenarios. Here we use rounding to avoid the possibility of tiny
    # floating point differences affecting the result
    sampled_scenario_Mw_bin = match(round(mw_sampled_scenarios, 3), round(ptha18_unique_mw_values, 3))
    stopifnot(all(!is.na(sampled_scenario_Mw_bin)))

    # Make matrix to store each exceedance-rate curve
    NTS = length(threshold_stage_values) # Number of stage thresholds
    NLTB = nrow(rate_models_all_branches$all_par) # Number of logic tree branches
    exrate_by_stage_threshold = matrix(NA, nrow=NTS, ncol=NLTB)

    # For each logic-tree branch
    for(i in 1:NLTB){

        # Exceedance-rate at mw-bin boundaries
        # Must interpolate with rule=1 because even at the highest value of Mw_seq, all_rate_matrix[i,] 
        # may not be zero.
        mw_bb_exrate = approx(rate_models_all_branches$Mw_seq, rate_models_all_branches$all_rate_matrix[i,],
            xout=ptha18_unique_mw_bb, rule=1)$y
        k = is.na(mw_bb_exrate)
        if(any(k)) mw_bb_exrate[k] = 0

        # Rate of scenarios in each mw bin (corresponding to ptha18_unique_mw_values)
        mw_bin_rate = -diff(mw_bb_exrate)

        # map onto sampled scenarios
        sampled_scenario_mw_bin_rate = mw_bin_rate[sampled_scenario_Mw_bin]

        # r_i(e) for sampled scenrios e
        #    = (rate-of-bin-containing-scenario-e) * (conditional-probability-of-scenario-e-in-that-bin)
        sampled_scenario_rate = sampled_scenario_mw_bin_rate * ptha18_conditional_prob_given_mw_sampled_scenarios

        # Get the importance-sampling-based exceedance rates, and an approximate variance
        #exrates_with_uncertainty = isu$estimate_exrate_and_variance_sampled_scenarios(
        #        rates_sampled_scenarios=sampled_scenario_rate,
        #        sampling_prob_sampled_scenarios=sampling_prob_sampled_scenarios,
        #        stage_sampled_scenarios = sampled_scenarios_max_stage,
        #        desired_stage_thresholds = threshold_stage_values)
        #exrate_by_stage_threshold[,i] = exrates_with_uncertainty[,1]

        for(j in 1:length(threshold_stage_values)){
            exrate_by_stage_threshold[j,i] = sum(
                sampled_scenario_rate * (sampled_scenarios_max_stage > threshold_stage_values[j])/
                (length(sampled_scenario_rate) * sampling_prob_sampled_scenarios))
        }
    }

    output = list(logic_tree_branch_exceedance_rates = exrate_by_stage_threshold, 
        threshold_stages = threshold_stage_values,
        logic_tree_branch_posterior_prob = rate_models_all_branches$all_par_prob)
}

# Plot the exceedance-rates for each source-zone separately in the PTHA18 and nonlinear model.
plot_epistemic_uncertainties_in_PTHA18_and_nonlinear_model_by_source_zone<-function(
    target_point,
    percentile_uncertainty_results,
    nonlinear_model_MSL,
    ptha18_percentile_uncertainty_results,
    max_stage_type = 'nonlinear_model'
    ){

    #
    # Plot epistemic uncertainties in PTHA18 and nonlinear model, by source-zone
    #
    if(max_stage_type == 'nonlinear_model'){
        # Real case
        out_file = paste0('Exceedance_Rates_with_epistemic_uncertainty_in_PTHA18_and_Nonlinear_Model_', 
            target_point[1], '_', target_point[2], '.png')
    }else if(max_stage_type == 'ptha18'){
        # Useful test case
        out_file = paste0('Exceedance_Rates_with_epistemic_uncertainty_in_PTHA18_and_sampled_ptha18_scenarios_', 
            target_point[1], '_', target_point[2], '.png')
    }else{
        stop(paste0('unknown max_stage_type: ', max_stage_type))
    }


    NPANEL = length(names(percentile_uncertainty_results))
    HT = 7.5
    WD = HT * NPANEL
    png(out_file, width=WD, height=HT, units='in', res=300)
    par(mfrow=c(1,NPANEL))
    options(scipen=5)
    nr = nrow(percentile_uncertainty_results[[source_zone]]$percentile_exrate)

    for(source_zone in names(percentile_uncertainty_results)){

        plot(percentile_uncertainty_results[[source_zone]]$threshold_stages - nonlinear_model_MSL, 
             percentile_uncertainty_results[[source_zone]]$percentile_exrate[nr,], 
             ylim=c(1.0e-05, 0.1), t='o', log='xy', col='white')
        for(i in 1:nr){
             points(percentile_uncertainty_results[[source_zone]]$threshold_stages - nonlinear_model_MSL, 
                    percentile_uncertainty_results[[source_zone]]$percentile_exrate[i,], 
                    t='o', col=i)
        }

        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate_lower_ci, 
               t='l', lty='dotted')
        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate_16pc, 
               t='l', lty='dotted')
        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate_median, 
               t='l', lty='dotted')
        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate_84pc, 
               t='l', lty='dotted')
        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate_upper_ci, 
               t='l', lty='dotted')

        if(max_stage_type == 'nonlinear_model'){
            mytitle = 'Max-stage exceedance-rate with epistemic uncertainty \n PTHA18 & nonlinear model: '
        }else if(max_stage_type == 'ptha18'){
            mytitle = 'Max-stage exceedance-rate with epistemic uncertainty \n PTHA18 & sampled ptha18: '
        }else{
            stop(paste0('unknown max_stage_type: ', max_stage_type))
        }

        title(paste0(mytitle, source_zone), cex.main=1.6)
        legend('topright', c('2.5%', '16%', '50%', '84%', '97.5%', 'PTHA18'), 
               col=c(1:nr, 1), lty=c(rep('solid', nr), 'dotted'), pch=c(rep(1, nr), NA))
    }
    dev.off()
}

# Plot the exceedance-rates summed over all source zones in the PTHA18 and the nonlinear model.
plot_summed_exrates_with_epistemic_uncertainty_in_PTHA18_and_nonlinear_model<-function(
    target_point,
    percentile_uncertainty_results,
    nonlinear_model_percentile_uncertainty,
    nonlinear_model_MSL,
    epistemic_uncertainty_threshold_stage_values,
    ptha18_percentile_uncertainty_results,
    ptha18_vals,
    max_stage_type
    ){

    if(max_stage_type == 'nonlinear_model'){
        # Real case
        out_file = paste0('Exceedance_Rates_Summed_with_epistemic_uncertainty_in_PTHA18_and_Nonlinear_Model_', 
            target_point[1], '_', target_point[2], '.png')
    }else if(max_stage_type == 'ptha18'){
        # Useful test case
        out_file = paste0('Exceedance_Rates_Summed_with_epistemic_uncertainty_in_PTHA18_and_sampled_ptha18_scenarios_', 
            target_point[1], '_', target_point[2], '.png')
    }else{
        stop(paste0('unknown max_stage_type: ', max_stage_type))
    }


    HT = 7.5
    png(out_file, width=HT*1.2, height=HT, units='in', res=300)
    options(scipen=5)
    nr = nrow(percentile_uncertainty_results[[source_zone]]$percentile_exrate)

    plot(percentile_uncertainty_results[[source_zone]]$threshold_stages - nonlinear_model_MSL, 
         nonlinear_model_percentile_uncertainty[nr,], 
         xlim=c(0.02, max(epistemic_uncertainty_threshold_stage_values)), ylim=c(1.0e-05, 0.1), 
         t='o', log='xy', col='white',
         xlab='Tsunami maxima (m above ambient sea-level)', 
         ylab='Exceedance-rate (events/year)', cex.axis=1.5, cex.lab=1.5)
    for(i in 1:nr){
         points(percentile_uncertainty_results[[source_zone]]$threshold_stages - nonlinear_model_MSL, 
                nonlinear_model_percentile_uncertainty[i,], 
                t='o', col=i)
    }
    for(nm in names(ptha18_vals)){
        points(ptha18_percentile_uncertainty_results[[1]]$stage,
               ptha18_vals[[nm]], t='l', lty='dotted')
    }

    if(max_stage_type == 'nonlinear_model'){
        # Real case
        mytitle = 'Max-stage exceedance-rate with epistemic uncertainties \n PTHA18 & nonlinear model, sum over source-zones'
    }else if(max_stage_type == 'ptha18'){
        mytitle = 'Max-stage exceedance-rate with epistemic uncertainties \n PTHA18 & sampled pth18 scenarios, sum over source-zones'
    }else{
        stop(paste0('unknown max_stage_type: ', max_stage_type))
    }

    title(mytitle, cex.main=1.7)
    legend('topright', c('2.5%', '16%', '50%', '84%', '97.5%', 'PTHA18'), 
           col=c(1:nr, 1), lty=c(rep('solid', nr), 'dotted'), pch=c(rep(1, nr), NA))

    dev.off()
}
