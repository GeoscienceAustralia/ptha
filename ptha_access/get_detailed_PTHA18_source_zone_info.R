#
# The code here works with the PTHA18 "compute_rates_all_sources_session.RData"
# R session, and so one can interrogate details of the source-zone rate models
# that are unavailable via the netcdf files.
#
# It looks for a local copy of "compute_rates_all_sources_session.RData" and
# downloads it if it is unavailable. Hence the first time you use it, it may
# take some time to source().
#


#
# Get the PTHA18 access codes, needed for functions below
#
get_ptha_script = './get_PTHA_results.R'
ptha18 = new.env()
source(get_ptha_script, local=ptha18, chdir=TRUE)

#
# Get the R session resulting from "compute_rates_all_sources.R", needed for
# functions below
#

# Get the saved R-session associated with the source-zone event-rate
# computation from PTHA18. 
compute_rates_session = './compute_rates_all_sources_session.RData'
if(!file.exists(compute_rates_session)){
    # If we are on NCI we might be able to get a copy from here
    compute_rates_session_NCI = paste0('/g/data/fj6/PTHA/AustPTHA_1/', 
        'EVENT_RATES/compute_rates_all_sources_session.RData')
    if(file.exists(compute_rates_session_NCI)){
        file.copy(compute_rates_session_NCI, compute_rates_session)
    }
}
if(!file.exists(compute_rates_session)){
    # If the file wasn't found in the above locations, then download it locally
    compute_rates_session_download = paste0('https://thredds.nci.org.au/',
        'thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/',
        'compute_rates_all_sources_session.RData')
    download.file(compute_rates_session_download, compute_rates_session, mode='wb')
}
crs_data = new.env()
load(compute_rates_session, envir=crs_data)

#' PTHA18 rates and conditional probabilities on a segment-by-segment basis
#'
#' For a given source-zone and segment, get the PTHA18 scenario rates and their
#' conditional probabilities [conditional on Mw]. Note in regular ptha_access
#' outputs we provide a mixture of segmented and unsegmented treatments,
#' whereas this function can isolate each segment, as well as the unsegmented
#' branch
#'
#' @param source_zone name of source-zone in PTHA18
#' @param segment the name of a segment associated with the source-zone, or ''
#' if referring to the unsegmented source.
#' @param mu_type 'constant' or 'variable_mu' as per PTHA18. Note that with both
#' methods, the event rates are always a function of the constant-rigidity-Mw, but
#' with different posterior GR curve weights and scenario conditional probabilities
get_PTHA18_scenario_conditional_probability_and_rates_on_segment<-function(
    source_zone, segment='', mu_type='constant'){

    if(segment != ''){
        source_zone_segment = paste0(source_zone, '_', segment)
    }else{
        source_zone_segment = source_zone
    }

    # Convenient shorthand to refer to the environment with the target
    # source-zones data
    sz = crs_data$source_envs[[source_zone_segment]]

    # Get available Mw values [constant rigidity case]. Should be
    #     7.2, 7.3, ...., 9.6, 9.7, 9.8 
    # but not all of these will have a non-zero rate.
    mws = sort(unique(sz$event_table$Mw))
    if(!all(abs(diff(mws) - 0.1) < 1.0e-06))
        stop('discretization of mws are not as expected')
    if(!all(abs(range(mws) - c(7.2, 9.8)) < 1.0e-06)) 
        stop('range of mws not as expected')


    # Get the 'uniform slip' scenario probabilities conditional on Mw.
    # Constant rigidity.
    if(mu_type == 'variable_mu'){
        FAUS_event_rates = sz$event_rates_mu_vary # Logic-tree mean, variable rigidity
    }else if(mu_type == 'constant'){
        FAUS_event_rates = sz$event_rates # Logic-tree mean
    }else{
        stop(paste0('Unknown mu_type= ', mu_type))
    }
    FAUS_conditional_probabilities = rep(0, nrow(sz$event_table))
    FAUS_mw = sz$event_table$Mw
    for(i in 1:length(mws)){
        k = which(FAUS_mw == mws[i])
        prob_with_Mw = sum(FAUS_event_rates[k])
        if(prob_with_Mw > 0){
            FAUS_conditional_probabilities[k] = FAUS_event_rates[k]/prob_with_Mw
        }
    }

    # For each 'parent' uniform slip event, we have a set (N=15) of 'child'
    # heterogeneous-slip events (HS) and variable_area_uniform_slip events (VAUS).
    # Their unequal rates will add up to give the associated uniform_slip (FAUS)
    # rates. 
    # If we divide by the associated FAUS rate (so the child-scenario
    # numbers add to one), then those numbers will not vary depending on the
    # segment [but beware NA values caused by division by zero, associated with
    # impossible scenarios]
    get_rates_and_uniform_event_row<-function(slip_type){

        if(!any(slip_type %in% c('stochastic', 'variable_uniform'))){
            stop('slip_type must be either "stochastic" or "variable_uniform"')
        }

        # Read the data from NCI
        nc_file  = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
            'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', slip_type,
            '_slip_earthquake_events_', source_zone, '.nc')
        fid = nc_open(nc_file, readunlim=FALSE)

        if(mu_type == 'variable_mu'){
            rates_full_source = ncvar_get(fid, 'variable_mu_rate_annual')
        }else if(mu_type == 'constant'){
            rates_full_source = ncvar_get(fid, 'rate_annual')
        }
        uniform_event_row = ncvar_get(fid, 'uniform_event_row')
        mw = ncvar_get(fid, 'Mw')
        nc_close(fid)

        unique_uniform_event_row = sort(unique(uniform_event_row))
        child_conditional_prob = rates_full_source*0
        parent_uniform_scenario_rate = rates_full_source*0
        for(i in 1:length(unique_uniform_event_row)){
            k = which(uniform_event_row == unique_uniform_event_row[i])
            parent_uniform_scenario_rate[k] = 
                sum(rates_full_source[k]) # Deliberately all the same
            if(parent_uniform_scenario_rate[k[1]] > 0){
                child_conditional_prob[k] = 
                    rates_full_source[k]/parent_uniform_scenario_rate[k]
            }
        }

        return(data.frame(
            rates=rates_full_source, 
            uniform_event_row=uniform_event_row,
            parent_uniform_scenario_rate=parent_uniform_scenario_rate,
            child_conditional_prob=child_conditional_prob,
            Mw = mw))
    }

    HS_data = get_rates_and_uniform_event_row('stochastic')
    VAUS_data = get_rates_and_uniform_event_row('variable_uniform')

    # Get scenario probabilities on our source/segment combination, conditional
    # on a scenario with the same magnitude having occurred.
    HS_prob_given_Mw = HS_data$child_conditional_prob * 
        FAUS_conditional_probabilities[HS_data$uniform_event_row]
    VAUS_prob_given_Mw = VAUS_data$child_conditional_prob *
        FAUS_conditional_probabilities[VAUS_data$uniform_event_row]
    # Get the logic-tree mean rates applied to the individual scenarios.
    HS_event_rates = HS_data$child_conditional_prob * 
        FAUS_event_rates[HS_data$uniform_event_row]
    VAUS_event_rates = VAUS_data$child_conditional_prob * 
        FAUS_event_rates[VAUS_data$uniform_event_row]

    output = list(FAUS_prob_given_Mw = FAUS_conditional_probabilities,
                  FAUS_event_rates = FAUS_event_rates,
                  FAUS_mw = FAUS_mw,
                  HS_prob_given_Mw = HS_prob_given_Mw,
                  HS_event_rates = HS_event_rates,
                  HS_mw = HS_data$Mw,
                  HS_uniform_event_row = HS_data$uniform_event_row,
                  VAUS_prob_given_Mw = VAUS_prob_given_Mw,
                  VAUS_event_rates = VAUS_event_rates,
                  VAUS_mw = VAUS_data$Mw,
                  VAUS_uniform_event_row = VAUS_data$uniform_event_row)

    return(output)

}

#' Get information on unsegmented/segmented models on a source zone
#'
#' Convenience function to get names of unsegmented/segmented source representations,
#' and their weights. This can be helpful for epistemic uncertainty calculations
#' where we want to work with logic-tree-branches in crs_data$source_envs[[NAME]].
#' 
#'
#' @param source_zone The name of a source zone in the PTHA18
#' @return a list with the unsegmented source name, it's weight, the segmented
#' source names and the weight assigned to the union-of-segments model
#' 
get_unsegmented_and_segmented_source_names_on_source_zone<-function(source_zone){

    all_source_names_unseg_and_seg = names(crs_data$source_envs)

    unsegmented_source_zone = source_zone
    # Ensure unique match
    stopifnot(sum(unsegmented_source_zone == all_source_names_unseg_and_seg) == 1)

    k = startsWith(all_source_names_unseg_and_seg, paste0(source_zone, '_'))
    segmented_source_zones = all_source_names_unseg_and_seg[k]

    # If there are segments then PTHA18 assigns 50:50 weight to unsegmented:segments
    # Otherwise PTHA18 assigns 100% weight to the unsegmented model.
    has_segments = (length(segmented_source_zones) > 0)
    if(has_segments){
        unsegmented_weight = 0.5
        union_of_segments_weight = 0.5
    }else{
        unsegmented_weight = 1
        union_of_segments_weight = 0
    } 

    stopifnot(isTRUE(all.equal(unsegmented_weight + union_of_segments_weight, 1.0)))

    output = list(unsegmented_name = unsegmented_source_zone,
        unsegmented_weight = unsegmented_weight,
        segments = segmented_source_zones,
        union_of_segments_weight=union_of_segments_weight)
    return(output)
}


#' Random scenario exceedance-rates over all logic-tree branches.
#'
#' Given random scenarios and their stage values on a source-zone/segment,
#' compute the exceedance-rates for all logic tree branches, and the sample variance,
#' given at set of threshold stage values. This can be done relatively
#' efficiently by noting that within a magnitude bin, the PTHA18 scenario
#' conditional probabilities do not vary between logic tree branches. (That
#' fact enables the math to be simplified, and is exploited below).
#'
#' @param source_zone name of the source zone (e.g. "kermadectonga2")
#' @param segment name of the segment (e.g. "tonga"). This should be "" for
#' unsegmented sources.
#' @param random_scenarios The random scenarios that result from
#' ptha18$randomly_sample_scenarios_by_Mw_and_rate, (applied to the current
#' source_zone). 
#' @param all_scenario_stage The stage values for ALL SCENARIOS THAT COULD HAVE
#' BEEN SAMPLED, such that all_scenario_stage[random_scenarios$inds] gives the
#' stage values for the randomly sampled scenarios. If the stage-values for
#' non-sampled scenarios are not known, they can be set to NA without a
#' problem.
#' @param threshold_stages A vector of stage values for which we compute
#' exceedance-rates
#' @param check_consistency_random_scenarios_rate_and_PTHA18_mean_rate logical
#' if TRUE we check that random_scenarios$rate_with_this_mw is consistent with
#' the logic-tree mean exceedance rates from PTHA18 results on this source-zone
#' and segment. The check can be turned off because it is possible to use this
#' routine in a valid way for which that would not be true (any choice is OK so
#' long as the within-magnitude-bin importance-sampling weights would be
#' unchanged). However for my normal usage this test should pass, so it is
#' included by default.
#' @param mu_type 'constant' or 'variable_mu' as per PTHA18. Note that with both
#' methods, the event rates are always a function of the constant-rigidity-Mw, but
#' with different posterior GR curve weights and scenario conditional probabilities
#' @return a list with the following entries:
#' - source_segment_name (combined source_zone + segment name, which is in
#'   names(crs$source_envs)), 
#  - unique_mw (magnitude values found in random_scenarios), 
#  - threshold_stages (same as input), 
#' - logic_tree_branch_exceedance_rates (matrix with an exceedance rate for each
#'   threshold stage and logic-tree branch),
#' - logic_tree_branch_exceedance_rates_var (matrix with a "single-sample
#'   Monte-Carlo variance" of the exceedance rate for each threshold stage and
#'   logic-tree branch),
#' - logic_tree_mean_exceedance_rates (vector with the "weighted mean over all
#'   logic-tree branches" of the exceedance-rates).
#' - logic_tree_branch_mw_bin_rates (matrix with the rate of events in each
#'   mw_bin, for each logic-tree branch),
#' - logic_tree_branch_posterior_prob (vector with the PTHA18 posterior
#'   probabilities or weights for each logic-tree branch)
#' - conditional_prob_exceed_stage_mw (matrix with the conditional probability
#'   of exceeding each threshold stage, given the mw-bin)
#' 
random_scenario_exceedance_rates_all_logic_tree_branches<-function(
    source_zone, 
    segment,
    random_scenarios,
    all_scenario_stage,
    threshold_stages,
    check_consistency_random_scenarios_rate_and_PTHA18_rates=TRUE,
    mu_type='constant'){

    if(segment == ''){
        source_segment_name = source_zone
    }else{
        source_segment_name = paste0(source_zone, '_', segment)
    }

    # Check that the source_segment_name is valid
    if(!(source_segment_name %in% names(crs_data$source_envs))){
        msg = paste0('Cannot find this source_zone and segment in the PTHA18 source-zone data: ',
                     '"', source_segment_name, '" .',
                     'Likely there is an error in "source_zone" and/or "segment". Note that for ',
                     'unsegmented results, segment should be "" (i.e. the length=0 character) ')
        stop(msg)
    }

    # Get mw-exceedance-rate information for all logic tree branches in the PTHA18
    all_branches = crs_data$source_envs[[source_segment_name]]$mw_rate_function(NA, 
        return_all_logic_tree_branches=TRUE)

    # Unique magnitudes that we have random scenarios for, and the magnitude-bin boundaries
    unique_mw = ptha18$unique_sorted_with_check_for_even_spacing(random_scenarios$mw)
    dMw = (unique_mw[2] - unique_mw[1])
    unique_mw_bin_boundaries = c(unique_mw - (dMw/2), max(unique_mw) + dMw/2)

    # Get the logic-tree-mean rates within each magnitude bin (conveniently stored in the random_scenarios)
    unique_mw_rates_source = random_scenarios$rate_with_this_mw[match(unique_mw, random_scenarios$mw)]

    if(check_consistency_random_scenarios_rate_and_PTHA18_rates){
        # Check it is consistent with the PTHA18 results in all_branches
        if(mu_type == 'variable_mu'){
            unique_mw_rates_source_check = -diff(crs_data$source_envs[[source_segment_name]]$mw_rate_function(unique_mw_bin_boundaries, account_for_mw_obs_error=TRUE))
        }else if(mu_type == 'constant'){
            unique_mw_rates_source_check = -diff(crs_data$source_envs[[source_segment_name]]$mw_rate_function(unique_mw_bin_boundaries))
        }else{
            stop(paste0('unknown mu_type ', mu_type))
        }
        if(any(abs( unique_mw_rates_source_check - unique_mw_rates_source) > 1.0e-05*unique_mw_rates_source)){
            stop('inconsistency between random_scenarios$rate_with_this_mw and the PTHA18 logic-tree mean mw-bin rates')
        }
    }

    inv_unique_mw_rates_source = 1/unique_mw_rates_source
    inv_unique_mw_rates_source[unique_mw_rates_source == 0] = 0 # Avoid zero-division

    # Compute the conditional probability of exceeding each threshold_stage, given the magnitude bin
    # IN PTHA18 THIS IS INDEPENDENT OF THE LOGIC-TREE-BRANCH
    conditional_prob_exceed_stage_mw = matrix(NA, ncol=length(unique_mw), nrow=length(threshold_stages))
    # Also compute the 'within-bin-variance'/'within-bin-rate^2' -- with can likewise be adjusted to give
    # variances for other mw-frequency models with the same scenario conditional probabilities.
    conditional_var_exceed_stage_mw = matrix(NA, ncol=length(unique_mw), nrow=length(threshold_stages))
    for(i in 1:length(threshold_stages)){

        exrate_by_mw_bin = ptha18$estimate_exrate_uncertainty(
            random_scenarios, 
            all_scenario_stage,
            threshold_stage = threshold_stages[i], 
            return_per_Mw_bin=TRUE)

        conditional_prob_exceed_stage_mw[i,] = 
            exrate_by_mw_bin$exrate*inv_unique_mw_rates_source

        conditional_var_exceed_stage_mw[i,] = 
            exrate_by_mw_bin$exrate_variance * (inv_unique_mw_rates_source**2)
    }

    # For each logic-tree branch, make a matrix to store the rates in each mw bin.
    # - columns correspond to different logic-tree branches,
    # - rows correspond to unique_mw
    num_logictree_branches = length(all_branches$all_par_prob)
    logic_tree_branch_mw_bin_rates = matrix(NA, 
        ncol=num_logictree_branches, nrow=length(unique_mw_rates_source))

    for(i in 1:num_logictree_branches){
        # Interpolate the exceedance-rate curve the magnitude-bin boundaries
        # (typically 7.15, 7.25, .... 9.55, 9.65)
        exrates_tmp = approx(all_branches$Mw_seq, all_branches$all_rate_matrix[i,], 
            xout=unique_mw_bin_boundaries, rule=1)$y

        # The above interpolation will be NA if any Mw value exceeds
        # all_branches$Mw_seq. The latter has a varying range depending on the
        # source representation. But we know the exceedance-rate is zero is
        # this case
        k = is.na(exrates_tmp)
        if(any(k)) exrates_tmp[k] = 0
        # Store the individual bin rates
        logic_tree_branch_mw_bin_rates[,i] = -diff(exrates_tmp)
    }

    logic_tree_branch_exceedance_rates = conditional_prob_exceed_stage_mw%*%logic_tree_branch_mw_bin_rates
    logic_tree_branch_exceedance_rates_var = conditional_var_exceed_stage_mw%*%(logic_tree_branch_mw_bin_rates^2)

    #
    # Quick checks
    #
    logic_tree_rates_of_any_event = colSums(logic_tree_branch_mw_bin_rates)  

    if(mu_type == 'variable_mu'){
        post_prob = all_branches$all_par_prob_with_Mw_error
        expected_val = (
            crs_data$source_envs[[source_segment_name]]$mw_rate_function(min(unique_mw_bin_boundaries), account_for_mw_obs_error=TRUE) -
            crs_data$source_envs[[source_segment_name]]$mw_rate_function(max(unique_mw_bin_boundaries), account_for_mw_obs_error=TRUE) )
    }else if(mu_type == 'constant'){
        post_prob = all_branches$all_par_prob
        expected_val = (
            crs_data$source_envs[[source_segment_name]]$mw_rate_function(min(unique_mw_bin_boundaries)) -
            crs_data$source_envs[[source_segment_name]]$mw_rate_function(max(unique_mw_bin_boundaries)) )
    }else{
        stop(paste0('unknown mu_type ', mu_type))
    }

    mean_rate_any_event = sum(logic_tree_rates_of_any_event * post_prob)
    # This should be the same (up to floating point) as 
    # For testing it is useful to compute the mean-curve implied by the logic-tree branches.
    # This should be the same as previous results (but errors in interpolation
    # and so on lead to differences -- this helped catch some errors while writing the code)
    back_computed_mean_curve = apply(logic_tree_branch_exceedance_rates, 1, 
        function(x) weighted.mean(x,w=post_prob) )

    stopifnot(abs(mean_rate_any_event - expected_val) <= 1.0e-06*expected_val)

    outputs = list(source_segment_name = source_segment_name,
                   unique_mw = unique_mw,
                   threshold_stages = threshold_stages,
                   logic_tree_branch_exceedance_rates = logic_tree_branch_exceedance_rates,
                   logic_tree_branch_exceedance_rates_var = logic_tree_branch_exceedance_rates_var,
                   logic_tree_mean_exceedance_rates = back_computed_mean_curve,
                   logic_tree_branch_mw_bin_rates = logic_tree_branch_mw_bin_rates,
                   logic_tree_branch_posterior_prob = post_prob,
                   conditional_prob_exceed_stage_mw = conditional_prob_exceed_stage_mw
                   )

    # Explicitly remove the large variables defined above (which in practice
    # can help R manage memory, even though documentation suggests we shouldn't need to)
    rm(logic_tree_branch_exceedance_rates, logic_tree_branch_mw_bin_rates, post_prob,
       all_branches, conditional_prob_exceed_stage_mw, conditional_var_exceed_stage_mw,
       logic_tree_branch_exceedance_rates_var)
    gc()

    return(outputs)
}

.test_kermadectonga2<-function(mu_type='constant'){

    # Read unsegmented and segmented sources from PTHA18 kermadectonga2 source
    source_zone = 'kermadectonga2'
    kt_full      = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, '', mu_type=mu_type)
    kt_tonga     = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'tonga', mu_type=mu_type)
    kt_kermadec  = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'kermadec', mu_type=mu_type)
    kt_hikurangi = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'hikurangi', mu_type=mu_type)

    # Double check consistency with PTHA18 files.

    # Weighted sum of 'unsegmented' and 'union of segments' rates as used in PTHA18
    back_calculated_HS_rates_combination = 
        0.5*(kt_full$HS_event_rates) + 
        0.5*(kt_tonga$HS_event_rates + 
             kt_kermadec$HS_event_rates + 
             kt_hikurangi$HS_event_rates)

    # Compare with PTHA18 file scenaro rates -- it should be the same to within floating
    # point 
    nc_file  = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
        'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', 'stochastic',
        '_slip_earthquake_events_', source_zone, '.nc')
    fid = nc_open(nc_file, readunlim=FALSE)
    if(mu_type == 'variable_mu'){
        rates_full_source = ncvar_get(fid, 'variable_mu_rate_annual')
    }else if(mu_type == 'constant'){
        rates_full_source = ncvar_get(fid, 'rate_annual')
    }else{
        stop(paste0('unknown mu_type ', mu_type))
    }
    nc_close(fid)

    err = back_calculated_HS_rates_combination-rates_full_source
    if(all(abs(err) < 1.0e-16)){
        print('PASS')
    }else{
        print('FAIL')
    }
    return(invisible(0))
}

.test_unsegmented_source<-function(mu_type='constant'){

    # Read unsegmented sources from PTHA18 kermadectonga2 source
    source_zone = 'puysegur2'
    sz_full = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, '', mu_type=mu_type)

    # The conditional probability of a given magnitude should either sum to
    # '1', or '0' for impossible magnitudes, up to floating point
    t1 = aggregate(sz_full$HS_prob_given_Mw, by=list(sz_full$HS_mw), sum)
    tester = (t1$x == 0) | (abs(t1$x - 1) < 1.0e-14)
    if(all(tester)){
        print('PASS')
    }else{
        print('FAIL')
    }
    # As above for VAUS scenarios
    t2 = aggregate(sz_full$VAUS_prob_given_Mw, by=list(sz_full$VAUS_mw), sum)
    tester = (t2$x == 0) | (abs(t2$x - 1) < 1.0e-14)
    if(all( tester )){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Consistent HS and FAUS rates
    t1 = aggregate(sz_full$HS_event_rates, by=list(sz_full$HS_uniform_event_row), sum)
    t2 = sz_full$FAUS_event_rates[t1[,1]]
    err = abs(t1$x - t2)
    if(all(err < 1.0e-16)){
        print('PASS')
    }else{
        print('FAIL')
    }

    # Consistent VAUS and FAUS rates
    t1 = aggregate(sz_full$VAUS_event_rates, by=list(sz_full$VAUS_uniform_event_row), sum)
    t2 = sz_full$FAUS_event_rates[t1[,1]]
    err = abs(t1$x - t2)
    if(all(err < 1.0e-16)){
        print('PASS')
    }else{
        print('FAIL')
    }
}

.test_random_scenario_exceedance_rate_percentile_calculation<-function(mu_type='constant'){
    #
    # Test the random scenario percentile uncertainty calculation code
    # We compute the stage-vs-exceedance-rate curve at a specified point, using
    # random scenarios, and compare with the original PTHA18 results (calculated
    # with a different technique, as available on the NCI THREDDS Server)
    #

    source_zone = 'kermadectonga2'
    target_point = c(185.1239, -21.0888) # Known location of PTHA18 hazard point 

    # Get all the scenarios, and make a convenient shorthand for the magnitudes
    # and rates in the event table
    kt2_scenarios = ptha18$get_source_zone_events_data(source_zone=source_zone, 
        slip_type='stochastic') 
    event_Mw = kt2_scenarios$events$Mw 
    if(mu_type == 'variable_mu'){
        event_rates = kt2_scenarios$events$variable_mu_rate_annual
    }else if(mu_type == 'constant'){
        event_rates = kt2_scenarios$events$rate_annual
    }else{
        stop(paste0('unknown mu_type ', mu_type))
    }

    # Get the exceedance-rate info for the unsegmented + segmented branches
    source_segment_names = paste0(source_zone, 
        c('_unsegmented', '_tonga_segment', '_kermadec_segment', '_hikurangi_segment'))
    segment_names = c('', 'tonga', 'kermadec', 'hikurangi')
    all_rate_models = vector(mode='list', length=length(source_segment_names))
    names(all_rate_models) = source_segment_names
    for(i in 1:length(all_rate_models)){
        all_rate_models[[source_segment_names[i]]] = 
            get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
                source_zone, segment_names[i], mu_type)
    }

    # Get the event peak stage at the target point
    event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
        target_point = target_point, 
        slip_type='stochastic',
        all_source_names=source_zone)
    # Convenient shorthand
    event_peak_stage = event_peak_stage_at_refpoint[[source_zone]]$max_stage

    # Get the PTHA18 stage exceedance-rate curve, for this source zone only, at the target point
    # This contains percentile information that we will test against
    stage_exrate_curve = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
        target_index=event_peak_stage_at_refpoint[[source_zone]]$target_index,
        source_name=source_zone)


    # Compute random scenarios for each source representation
    event_importance_weighted_sampling_probs = event_rates * event_peak_stage
    samples_per_Mw = function(Mw){ 4000 }
    mw_limits = c(7.15, 9.65)

    # Make a list that stores the random scenarios for unsegmented + each
    # segment.  We want to have the same scenarios in all cases (because for
    # inundation hazard studies, this is one way to gain computational
    # efficiency).
    random_scenarios = vector(mode='list', length=length(all_rate_models))
    names(random_scenarios) = names(all_rate_models)
    for(nm_i in names(random_scenarios)){
        set.seed(1234) # Same scenarios for all cases (keep resetting this)
        random_scenarios[[nm_i]] = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
            event_rates=all_rate_models[[nm_i]]$HS_event_rates,
            event_Mw = event_Mw,
            event_importance_weighted_sampling_probs = event_importance_weighted_sampling_probs,
            samples_per_Mw = samples_per_Mw,
            mw_limits = mw_limits
        )
    }

    # Choose stage values at which we will evaluate the exceedance-rate percentiles
    threshold_stage_values = stage_exrate_curve$stage

    # Compute exceedance-rates from random scenarios for all logic-tree branches: Unsegmented model
    unsegmented_scenario_exrates_logic_tree = random_scenario_exceedance_rates_all_logic_tree_branches(
        source_zone=source_zone,
        segment='',
        random_scenarios = random_scenarios$kermadectonga2_unsegmented,
        all_scenario_stage = event_peak_stage, 
        threshold_stages = threshold_stage_values,
        mu_type=mu_type)

    # Compute exceedance-rates from random scenarios for all logic-tree branches: all segment models
    segment_inds = which(!grepl('unsegmented', names(random_scenarios))) # Indices of segments in the lists above
    segmented_scenario_exrates_logic_tree = vector(mode='list', length=length(segment_inds))
    names(segmented_scenario_exrates_logic_tree) = source_segment_names[segment_inds]
    for(i in segment_inds){
        nm_i = names(random_scenarios)[i]
        segmented_scenario_exrates_logic_tree[[nm_i]] = random_scenario_exceedance_rates_all_logic_tree_branches(
            source_zone=source_zone,
            segment=segment_names[i],
            random_scenarios = random_scenarios[[nm_i]],
            all_scenario_stage = event_peak_stage, 
            threshold_stages = threshold_stage_values,
            mu_type=mu_type)
    }

    percentile_uncertainty_results = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
        unsegmented_scenario_exrates_logic_tree,
        segmented_scenario_exrates_logic_tree,
        N = 1e+06,
        unsegmented_wt=0.5,
        union_of_segments_wt=0.5,
        segments_copula_type = 'comonotonic',
        print_progress_every_nth_threshold=99999)

    check_similar<-function(p1, p2, reltol, which_inds){
        if(all(abs(p1[which_inds]-p2[which_inds]) < reltol*abs(p2[which_inds]))){
            print('PASS')
        }else{
            print('FAIL -- tolerance exceeded')
        }
    }
    
    # Test that the values are similar to the original PTHA18.
    # A plot shows the agreement is very good. 
    # At some points there are small differences. I expect in part these
    # reflect different calculation techniques in the original PTHA18 (where
    # the combination of the unsegmented/segmented distributions was done with
    # a numerical interpolation procedure, rather than the sampling herein).

    #plot(percentile_uncertainty_results$threshold_stages, percentile_uncertainty_results$mean_exrate, t='p', log='xy')
    #points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate, t='l', col='red')
    nm = ifelse(mu_type == 'constant', 'stochastic_slip_rate', 'variable_mu_stochastic_slip_rate')
    check_similar(percentile_uncertainty_results$mean_exrate, 
                  stage_exrate_curve[[nm]], 
                  reltol=0.06, 
                  which_inds = which(stage_exrate_curve[[nm]] > 1.0e-05) )

    #points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_84pc, t='l', col='blue')
    #points(percentile_uncertainty_results$threshold_stages, percentile_uncertainty_results$percentile_exrate[4,], col='blue')
    nm = ifelse(mu_type == 'constant', 'stochastic_slip_rate_84pc', 'variable_mu_stochastic_slip_rate_84pc')
    check_similar(percentile_uncertainty_results$percentile_exrate[4,], 
                  stage_exrate_curve[[nm]], 
                  reltol=0.05, 
                  which_inds = which(stage_exrate_curve[[nm]] > 1.0e-05) )

    #points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_16pc, t='l', col='blue')
    #points(percentile_uncertainty_results$threshold_stages, percentile_uncertainty_results$percentile_exrate[2,], col='blue')
    nm = ifelse(mu_type == 'constant', 'stochastic_slip_rate_16pc', 'variable_mu_stochastic_slip_rate_16pc')
    check_similar(percentile_uncertainty_results$percentile_exrate[2,], 
                  stage_exrate_curve[[nm]], 
                  reltol=0.12, 
                  which_inds = which(stage_exrate_curve[[nm]] > 1.0e-05) )

    #points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_upper_ci, t='l', col='brown')
    #points(percentile_uncertainty_results$threshold_stages, percentile_uncertainty_results$percentile_exrate[5,], col='brown')
    nm = ifelse(mu_type == 'constant', 'stochastic_slip_rate_upper_ci', 'variable_mu_stochastic_slip_rate_upper_ci')
    check_similar(percentile_uncertainty_results$percentile_exrate[5,], 
                  stage_exrate_curve[[nm]], 
                  reltol=0.14, 
                  which_inds = which(stage_exrate_curve[[nm]] > 1.0e-05) )

    #points(stage_exrate_curve$stage, stage_exrate_curve$stochastic_slip_rate_lower_ci, t='l', col='brown')
    #points(percentile_uncertainty_results$threshold_stages, percentile_uncertainty_results$percentile_exrate[1,], col='brown')
    nm = ifelse(mu_type == 'constant', 'stochastic_slip_rate_lower_ci', 'variable_mu_stochastic_slip_rate_lower_ci')
    check_similar(percentile_uncertainty_results$percentile_exrate[1,], 
                  stage_exrate_curve[[nm]], 
                  reltol=0.16, 
                  which_inds = which(stage_exrate_curve[[nm]] > 1.0e-05) )

    # Alternative approach to the computation
    percentile_uncertainty_results2 = ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
        unsegmented_scenario_exrates_logic_tree,
        segmented_scenario_exrates_logic_tree,
        N = 1e+06,
        unsegmented_wt=0.5,
        union_of_segments_wt=0.5,
        segments_copula_type = 'comonotonic',
        print_progress_every_nth_threshold=99999,
        use_numerical_probs=FALSE)

    # Check that the alternative percentile uncertainty calculation is very close to the other one
    check_similar(percentile_uncertainty_results$percentile_exrate,
                  percentile_uncertainty_results2$percentile_exrate, 
                  reltol=0.02,
                  which_inds=which(percentile_uncertainty_results$percentile_exrate > 1e-20))
}


test_get_detailed_PTHA18_source_zone_info<-function(){
    .test_kermadectonga2(mu_type='constant')
    .test_unsegmented_source(mu_type='constant')
    .test_random_scenario_exceedance_rate_percentile_calculation(mu_type='constant')

    .test_kermadectonga2(mu_type='variable_mu')
    .test_unsegmented_source(mu_type='variable_mu')
    .test_random_scenario_exceedance_rate_percentile_calculation(mu_type='variable_mu')
}
