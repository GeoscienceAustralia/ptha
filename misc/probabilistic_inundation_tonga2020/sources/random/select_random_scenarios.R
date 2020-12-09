#
# Pick random scenarios for a source zone with an importance-sampling approach.
#
# The idea is to sample scenarios by magnitude -- with a specified number
# scenarios per magnitude-bin, which may vary with Mw -- using sampling
# probabilities proportional to:
#     scenario_probability_conditional_on_Mw * scenario_size_indicator
# where:
#    - scenario_probability_conditional_on_Mw is based on PTHA18, and
#    - scenario_size_indicator (> 0) is user defined value used to give preferential
#      sampling to some scenarios.
#
# If scenario_size_indicator=1 for all scenarios, then we have simple random
# sampling, and for a given Mw we can assign an equal rate to all sampled
# scenarios:
#    nominal_scenario_rate = [ rate_of_event_with_the_same_Mw / 
#                              Number_of_scenarios_in_magnitude_bin ]
# However, if scenario_size_indicator is non-uniform, we can assign a nominal rate to
# each scenario as:
#    nominal_scenario_rate = [ rate_of_event_with_the_same_Mw * 
#                              1/scenario_size_indicator ] / sum(1/scenario_size_indicator)
#
# Herein we use the peak-stage near the site of interest as a scenario_size_indicator, which
# leads to more sampling of scenarios that are likely to produce significant inundation.




#
# Get the PTHA18 access codes, needed for functions below
#
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

#
# Get the R session resulting from "compute_rates_all_sources.R", needed for functions below
#

# Get the saved R-session associated with the source-zone event-rate computation from PTHA18. 
compute_rates_session = './compute_rates_all_sources_session.RData'
if(!file.exists(compute_rates_session)){
    # If we are on NCI we might be able to get a copy from here
    compute_rates_session_NCI = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES/compute_rates_all_sources_session.RData'
    if(file.exists(compute_rates_session_NCI)){
        file.copy(compute_rates_session_NCI, compute_rates_session)
    }
}
if(!file.exists(compute_rates_session)){
    # If the file wasn't found in the above locations, then download it locally
    compute_rates_session_download = 
        'http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/AustPTHA_1/EVENT_RATES/compute_rates_all_sources_session.RData'
    download.file(compute_rates_session_download, compute_rates_session)
}
crs_data = new.env()
load(compute_rates_session, envir=crs_data)


#' For a given source-zone and segment, get the PTHA18 scenario rates and their conditional
#' probabilities [conditional on Mw]. Note in regular ptha_access outputs we provide a mixture
#' of segmented and unsegmented treatments, whereas this function can isolate each segment, as 
#' well as the unsegmented branch
#'
#' @param source_zone name of source-zone in PTHA18
#' @param segment the name of a segment associated with the source-zone, or '' if referring to the unsegmented source.
#'
#' Note that in PTHA18, the 'full' source is often combination of segmented and
#' unsegmented branches -- whereas here we separate them.
get_PTHA18_scenario_conditional_probability_and_rates_on_segment<-function(source_zone, segment=''){

    if(segment != ''){
        source_zone_segment = paste0(source_zone, '_', segment)
    }else{
        source_zone_segment = source_zone
    }

    # Convenient shorthand to refer to the environment with the target source-zones data
    sz = crs_data$source_envs[[source_zone_segment]]

    # Get available Mw values [constant rigidity case]. Should be
    #     7.2, 7.3, ...., 9.6, 9.7, 9.8 
    # but not all of these will have a non-zero rate.
    mws = sort(unique(sz$event_table$Mw))
    if(!all(abs(diff(mws) - 0.1) < 1.0e-06)) stop('discretization of mws are not as expected')
    if(!all(abs(range(mws) - c(7.2, 9.8)) < 1.0e-06)) stop('range of mws not as expected')


    # Get the 'uniform slip' scenario probabilities conditional on Mw.
    # Constant rigidity.
    FAUS_event_rates = sz$event_rates # Logic-tree mean
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
        rates_full_source = ncvar_get(fid, 'rate_annual')
        uniform_event_row = ncvar_get(fid, 'uniform_event_row')
        mw = ncvar_get(fid, 'Mw')
        nc_close(fid)

        unique_uniform_event_row = sort(unique(uniform_event_row))
        child_conditional_prob = rates_full_source*0
        parent_uniform_scenario_rate = rates_full_source*0
        for(i in 1:length(unique_uniform_event_row)){
            k = which(uniform_event_row == unique_uniform_event_row[i])
            parent_uniform_scenario_rate[k] = sum(rates_full_source[k]) # Deliberately all the same
            if(parent_uniform_scenario_rate[k[1]] > 0){
                child_conditional_prob[k] = rates_full_source[k]/parent_uniform_scenario_rate[k]
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

    ## These should either be '1' (or '0' for impossible magnitudes).
    # aggregate(HS_prob_given_Mw, by=list(HS_data$Mw), sum)
    # aggregate(VAUS_prob_given_Mw, by=list(VAUS_data$Mw), sum)

    # Below are some plots used for ad-hoc check of the code [should ONLY work
    # for TOTALLY UNSEGMENTED source-zones]. It will not work where PTHA18 puts
    # partial weight on segmentation.
    if(FALSE){

        # For a totally UNSEGMENTED source-zone [e.g. puysegur2, newguinea2 ],
        # the following should plot on the 1:1 line. I confirmed it does for
        # "puysegur2" and "newguinea2".
        plot(HS_data$parent_uniform_scenario_rate, FAUS_event_rates[HS_data$uniform_event_row])
        grid(); abline(0, 1, col='red')
        # If the source-zone includes partial weight on a segmented
        # interpretation (e.g. as do most of our large source-zones), then it
        # won't plot on the 1:1 line in general, because the "file-values" read
        # above represent a mixture of the different segment interpretations

        # Absolute error -- should be within round-off on totally UNSEGMENTED source-zones
        # I confirmed this holds for puysegur2 and newguinea2
        plot((HS_data$parent_uniform_scenario_rate - FAUS_event_rates[HS_data$uniform_event_row]))

        ## Relative error -- should be consistent with round-off on totally
        ## UNSEGMENTED source-zones, beware near-zero division.
        ## I confirmed this holds for puysegur2 and newguinea2
        plot((HS_data$parent_uniform_scenario_rate - FAUS_event_rates[HS_data$uniform_event_row])/
             HS_data$parent_uniform_scenario_rate)

        # Plots as above, for VAUS
        plot(VAUS_data$parent_uniform_scenario_rate, FAUS_event_rates[VAUS_data$uniform_event_row])
        grid(); abline(0, 1, col='red')

        plot((VAUS_data$parent_uniform_scenario_rate - FAUS_event_rates[VAUS_data$uniform_event_row]))
        plot((VAUS_data$parent_uniform_scenario_rate - FAUS_event_rates[VAUS_data$uniform_event_row])/
             VAUS_data$parent_uniform_scenario_rate)
    }

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

#
# Sample scenarios by splitting according to magnitude, 
#
randomly_sample_scenarios_by_Mw_and_size<-function(
    event_rates,
    event_Mw,
    event_size_indicator,
    samples_per_Mw=function(x){ 50 + 0*x} ,
    mw_limits=c(7.15, 9.85)){

    unique_Mw = sort(unique(event_Mw))

    diff_unique_Mw = diff(unique_Mw)

    if(any(abs(diff_unique_Mw - diff_unique_Mw[1]) > 1.0e-06) ){
        stop('event_Mw is not evenly spaced')
    }

    # Ignore small scenarios, and scenarios exceeding Mw-max
    unique_Mw = unique_Mw[(unique_Mw > mw_limits[1] & unique_Mw < mw_limits[2])] 

    random_scenario_info = lapply(unique_Mw,
        f<-function(mw){
            # Match Mw [careful with real numbers]
            k = which(abs(event_Mw - mw) < 1.0e-03)
            if(length(k) == 1) stop('error: Only one scenario -- need to be careful using sample below')

            # Here we have 'n' scenarios that can be treated as an equal-weight representation of the
            # distribution conditional on the magnitude being mw.
            #sample_of_k = sample(k, size=150, prob=event_rates[k], replace=TRUE)

            nsam = samples_per_Mw(mw)
            local_sample = sample(1:length(k), size=nsam, 
                prob=event_rates[k]*event_size_indicator[k], 
                replace=TRUE)
            sample_of_k = k[local_sample]

            # The original scenario conditional probability distribution
            dist_f = event_rates[k]/sum(event_rates[k])
            # The distribution we sampled from above
            dist_g = (event_rates[k]*event_size_indicator[k])/sum(event_rates[k]*event_size_indicator[k])
            # The exact importance-sampling correction -- while these weights do not sum to 1, they
            # make the estimators unbiased
            alternate_weights = (dist_f[local_sample]/dist_g[local_sample])/length(local_sample)

            # Get the rate of any event with this mw
            rate_with_this_mw = sum(event_rates[k])

            # Importance sampling correction -- this is the 'self-normalised' approach
            # Estimate a rate for each individual scenario such that the sampled
            # scenario weights add to the target weight, and reflect the original distribution
            # Their chance of being sampled was inflated by 'event_size_indicator', so we deflate
            # by 'event_size_indicator' here.
            random_scenario_rates = rate_with_this_mw * (1/event_size_indicator[sample_of_k]) / 
                sum((1/event_size_indicator[sample_of_k]))

            # Importance sampling correction -- although "alternate_weights" does not sum to 1,
            # the self normalised importance sampling leads to some bias, whereas this approach is
            # unbiased.
            alternate_random_scenario_rates = rate_with_this_mw * alternate_weights

            return(list(
                inds=sample_of_k, 
                weight_inflation=event_size_indicator[sample_of_k], 
                random_scenario_rates = random_scenario_rates,
                alternate_random_scenario_rates = alternate_random_scenario_rates,
                rate_with_this_mw = rate_with_this_mw,
                mw = mw))
        })

    names(random_scenario_info) = as.character(unique_Mw)

    return(random_scenario_info)
}


# Get the scenario rates for each model.
source_zone = 'kermadectonga2'
kt_full      = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, '')
kt_tonga     = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'tonga')
kt_kermadec  = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'kermadec')
kt_hikurangi = get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, 'hikurangi')

check_consistency_with_PTHA18<-function(){
    # Double check consistency with PTHA18 files.

    # Weighted sum of 'unsegmented' and 'union of segments' rates as used in PTHA18
    back_calculated_HS_rates_combination = 
        0.5*(kt_full$HS_event_rates) + 
        0.5*(kt_tonga$HS_event_rates + 
             kt_kermadec$HS_event_rates + 
             kt_hikurangi$HS_event_rates)

    # Compare with PTHA18 file scenaro rates -- it should be the same.
    nc_file  = paste0(ptha18$config_env$.GDATA_OPENDAP_BASE_LOCATION,
        'SOURCE_ZONES/', source_zone, '/TSUNAMI_EVENTS/all_', 'stochastic',
        '_slip_earthquake_events_', source_zone, '.nc')
    fid = nc_open(nc_file, readunlim=FALSE)
    rates_full_source = ncvar_get(fid, 'rate_annual')
    nc_close(fid)

    # All going well, these two things should be equal to within acceptable
    # floating-point differences. Yes, they are.
    plot(back_calculated_HS_rates_combination, rates_full_source)
    abline(0, 1, col='red')
    # Plot the difference -- should be super small -- yes it is.
    plot(back_calculated_HS_rates_combination - rates_full_source)
}
pdf(paste0('check_consistency_rates_ptha18_', source_zone, '.pdf'), width=9, height=6)
check_consistency_with_PTHA18()
dev.off()


# Get the peak-stage at a point near Tonga -- we will use this to measure the 
# scenario 'size' and influence sampling.
event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3458.3,
    slip_type='stochastic',
    all_source_names=source_zone)

# We will use these nearby points for testing -- can we reproduce the stage-vs-rate curve
# with the randomly sampled scenarios? 
alternative1_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3454.3,
    slip_type='stochastic',
    all_source_names=source_zone)
alternative2_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3456.3,
    slip_type='stochastic',
    all_source_names=source_zone)
alternative3_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3455.3,
    slip_type='stochastic',
    all_source_names=source_zone)


#
# Make a plot of max-stage vs exceedance-rate according to the random sample of scenarios,
# and compare with the full set of scenarios. 
#
# Do this while using event_size_indicator = event_peak_stage**POW_STAGE for a range of POW_STAGE values.
# 
# Aim is to reproduce the full results well, with many fewer scenarios.
# 
plot_scenario_stages_with_different_scenario_size_weightings<-function(
    event_peak_stage, 
    event_Mw, 
    event_rates,
    samples_per_Mw = function(mw){ 50 + 0*mw},
    event_peak_stage_alt1 = alternative1_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    event_peak_stage_alt2 = alternative2_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    event_peak_stage_alt3 = alternative3_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    pdf_output_name_extra = '',
    POW_STAGE_VALUES = c(0, 1, 1.3)
    ){

    pdf(paste0('Effect_of_weighting_power_on_samples_fig_', pdf_output_name_extra, '.pdf'), 
        width=12, height=9)

    # Here we generate scenarios from Mw=7.5 up to whatever is the upper-limit
    # with non-zero rate on the source-zone.
    max_mw_limit = max(event_Mw[event_rates > 0]) + 0.05
    min_mw_limit = 7.45

    # Try a range of event_size_indicators, based on a power-law transformation of the event_peak_stage
    for(POW_STAGE in POW_STAGE_VALUES){

        set.seed(123) # Reproducible random seed

        event_size_indicator = event_peak_stage**POW_STAGE
        #event_size_indicator = pmin(event_peak_stage, 10)**POW_STAGE

        random_scenario_info = randomly_sample_scenarios_by_Mw_and_size(
            event_rates,
            event_Mw,
            event_size_indicator,
            samples_per_Mw = samples_per_Mw,
            mw_limits = c(min_mw_limit, max_mw_limit))

        random_scenario_inds = lapply(random_scenario_info, f<-function(x) x$inds)
        # Look at the stage values at the gauge offshore of Tongatapu.
        random_scenario_stages = lapply(random_scenario_info, f<-function(x) event_peak_stage[x$inds])
        random_scenario_rates = lapply(random_scenario_info, f<-function(x) x$random_scenario_rates)
        alternate_random_scenario_rates = lapply(random_scenario_info, f<-function(x) x$alternate_random_scenario_rates)
        random_scenario_Mws = lapply(random_scenario_info, f<-function(x) x$inds * 0 + x$mw)

        # Plot the distribution of max-stage for each Mw
        par(mfrow=c(4,6))
        mapply(f<-function(random_scenario_peak_stage, random_scenario_rates, random_scenario_Mws){
            percentiles = seq(0.001, 0.999, by=0.001)
            stages = weighted_percentile(random_scenario_peak_stage, random_scenario_rates, percentiles)

            k = which(event_Mw == random_scenario_Mws[1])
            stages_ref = weighted_percentile(event_peak_stage[k], event_rates[k], percentiles)

            plot(percentiles, stages, t='o', main=random_scenario_Mws[1], 
                 ylim=c(0.01, max(0.02, max(stages_ref))), log='')
            grid()
            abline(h=c(1, 5, 10), col='red')
            points(seq(0.001, 0.999, len=length(random_scenario_peak_stage)), 
                   sort(random_scenario_peak_stage), t='l', col='grey')

            # Add the 'correct' weighted percentile using all PTHA18 scenarios
            points(percentiles, stages_ref, t='l', col='green')

            }, random_scenario_stages, random_scenario_rates, random_scenario_Mws)

        #
        # Check if we can reproduce max-stage-vs-exceedance-rate curves
        #
        exrate_plot<-function(stage, rate, add=FALSE, ...){

            all_stages = seq(min(stage)-0.01, max(stage)-0.01, len=10000)
            all_rates = sapply(all_stages, f<-function(x) sum(rate*(stage>x)))

            if(!add){
                plot(all_stages, all_rates, log='xy', xlim=c(0.01, 20), ylim=c(1e-06, 1), ...)
                abline(v=c(1, 5, 10), col='red')
                abline(h=1/c(100, 1000, 10000), col='red')
            }else{
                points(all_stages, all_rates, ...)
            }
        }

        # Check if we can reproduce the stage-vs-exceedance-rate curve
        par(mfrow=c(2,2))
        exrate_plot(event_peak_stage, event_rates, t='l', 
                    main=paste0('Reference point; POW_STAGE=', POW_STAGE, '; Nunique=', 
                                length(unique(unlist(random_scenario_inds)))))
        exrate_plot(unlist(random_scenario_stages), unlist(random_scenario_rates), 
                    add=TRUE, col='green', t='l')
        exrate_plot(unlist(random_scenario_stages), unlist(alternate_random_scenario_rates), 
                    add=TRUE, col='blue', t='l')

        # Check a nearby site.
        random_scenario_stages_alt1 = lapply(random_scenario_inds, f<-function(x) event_peak_stage_alt1[x])
        exrate_plot(event_peak_stage_alt1, event_rates, t='l', main='Alternative 1')
        exrate_plot(unlist(random_scenario_stages_alt1), unlist(random_scenario_rates), 
                    add=TRUE, col='green', t='l')
        exrate_plot(unlist(random_scenario_stages_alt1), unlist(alternate_random_scenario_rates), 
                    add=TRUE, col='blue', t='l')
        # An another nearby site
        random_scenario_stages_alt2 = lapply(random_scenario_inds, f<-function(x) event_peak_stage_alt2[x])
        exrate_plot(event_peak_stage_alt2, event_rates, t='l', main='Alternative 2')
        exrate_plot(unlist(random_scenario_stages_alt2), unlist(random_scenario_rates), 
                    add=TRUE, col='green', t='l')
        exrate_plot(unlist(random_scenario_stages_alt2), unlist(alternate_random_scenario_rates), 
                    add=TRUE, col='blue', t='l')
        # An another nearby site
        random_scenario_stages_alt3 = lapply(random_scenario_inds, f<-function(x) event_peak_stage_alt3[x])
        exrate_plot(event_peak_stage_alt3, event_rates, t='l', main='Alternative 3')
        exrate_plot(unlist(random_scenario_stages_alt3), unlist(random_scenario_rates), 
                    add=TRUE, col='green', t='l')
        exrate_plot(unlist(random_scenario_stages_alt3), unlist(alternate_random_scenario_rates), 
                    add=TRUE, col='blue', t='l')

        par(mfrow=c(2,2))
        stages = c(1, 2, 5, 10)
        for(i in 1:length(stages)){
            percentiles = seq(0.001, 0.999, by=0.001)
            # Find possible events that exceed stages[i]
            k = which(event_peak_stage > stages[i] & event_rates > 0)

            # Plot if we can
            if(length(k) > 0){
                mws = weighted_percentile(event_Mw[k], event_rates[k], percentiles)
                plot(percentiles, mws, t='l')

                title(paste0('Mw distribution with stage > ', stages[i]))
            }else{
                plot(percentiles, percentiles*0, col='white')
                title(paste0('No scenarios with stage > ', stages[i]))
            }

            # Add the random scenarios to the plot, if we can.
            k = which(unlist(random_scenario_stages) > stages[i])
            if(length(k) > 0){
                mws = weighted_percentile(unlist(random_scenario_Mws)[k], 
                    unlist(random_scenario_rates)[k], percentiles)
                points(percentiles, mws, t='l', col='red')
            }else{
                title(sub='no random scenarios were found', col.sub='red')
            }

        }

    }
    dev.off()
}

#
# Check how well we can reproduce the full 'max-stage vs exceedance-rate' curve
# on segmented/unsegmented branches respectively. This code will make plots.
#
event_peak_stage = event_peak_stage_at_refpoint[[source_zone]]$max_stage
event_Mw = event_peak_stage_at_refpoint[[source_zone]]$Mw
plot_scenario_stages_with_different_scenario_size_weightings(
    event_peak_stage, 
    event_Mw, 
    kt_full$HS_event_rates,
    # Mean of 30 samples per Mw, but more at high Mw
    samples_per_Mw=function(x){ 15 + 30*(x-7.45)/(9.65-7.45)},
    pdf_output_name_extra='unsegmented')
plot_scenario_stages_with_different_scenario_size_weightings(
    event_peak_stage, 
    event_Mw, 
    kt_tonga$HS_event_rates,
    # Mean of 20 samples per Mw, but more at high Mw
    samples_per_Mw=function(x){ 10 + 20*(x-7.45)/(9.65-7.45)},
    pdf_output_name_extra='tonga_segment')
plot_scenario_stages_with_different_scenario_size_weightings(
    event_peak_stage, 
    event_Mw, 
    kt_kermadec$HS_event_rates,
    # Mean of 12.5 samples per Mw, but more at high Mw
    samples_per_Mw=function(x){ 5 + 10*(x-7.45)/(9.65-7.45)},
    pdf_output_name_extra='kermadec_segment')
plot_scenario_stages_with_different_scenario_size_weightings(
    event_peak_stage, 
    event_Mw, 
    kt_hikurangi$HS_event_rates,
    # 3 samples per Mw
    samples_per_Mw=function(x){3+0*x},
    pdf_output_name_extra='hikurangi_segment')

#
# Make the samples, with more scenarios at higher magnitudes, and a "scenario_size"
# proportional to the event_peak_stage.
#
sampling_setup = list(
    'unsegmented' = list(event_rates = kt_full$HS_event_rates, 
                         samples_per_Mw = function(x){ 15 + 30*(x-7.45)/(9.65-7.45)}),
    'tonga' = list(event_rates = kt_tonga$HS_event_rates, 
                   samples_per_Mw = function(x){ 10 + 20*(x-7.45)/(9.65-7.45)}),
    'kermadec' = list(event_rates = kt_kermadec$HS_event_rates, 
                      samples_per_Mw = function(x){ 5 + 10*(x-7.45)/(9.65-7.45)}),
    'hikurangi' = list(event_rates = kt_hikurangi$HS_event_rates, 
                       samples_per_Mw = function(x) {3+0*x})
    )
random_scenario_info = vector(mode='list', length=length(sampling_setup))
names(random_scenario_info) = names(sampling_setup)
for(i in 1:length(random_scenario_info)){
    set.seed(123) # Reproducible random seed

    # Here are the magnitude limits where we bother sampling
    max_mw_limit = max(event_Mw[sampling_setup[[i]]$event_rates > 0]) + 0.05
    min_mw_limit = 7.45

    # Take the sample
    random_scenario_info[[i]] = randomly_sample_scenarios_by_Mw_and_size(
        sampling_setup[[i]]$event_rates,
        event_Mw,
        event_size_indicator=event_peak_stage**1.0,
        samples_per_Mw = sampling_setup[[i]]$samples_per_Mw,
        mw_limits = c(min_mw_limit, max_mw_limit))
}

# How many scenarios will we have to run?
all_scenarios = unlist(lapply(random_scenario_info, f<-function(x) unlist(lapply(x, f<-function(y) y$inds))))
length(unique(all_scenarios))
quantile(event_peak_stage[all_scenarios], probs=seq(0, 1, len=21))

#
# Write out the scenarios to a csv format
#
write_scenario_info<-function(random_scenario_info, filename){
    # Indices in the source-zone event table
    inds = unlist(lapply(random_scenario_info, f<-function(x) x$inds))        
    # The "size_factor" used for importance sampling
    weight_inflation = unlist(lapply(random_scenario_info, f<-function(x) x$weight_inflation))        
    # The scenario rates that should be used for calculations [corrected for
    # importance sampling, non-uniform conditional probs, etc]
    scenario_rates = unlist(lapply(random_scenario_info, f<-function(x) x$random_scenario_rates))        
    # Here are the 'unbiased' variants of the scenario_rates. They are not always exactly consistent
    # with the standard scenario rates [e.g. sum_of_rate_for_a_given_Mw is not exactly consistent with the PTHA18]
    # but this has the benefit of unbiased-ness when using them to estimate probability integrals.
    alternate_scenario_rates = unlist(lapply(random_scenario_info, f<-function(x) x$alternate_random_scenario_rates))        
    # The magnitudes
    mw = unlist(lapply(random_scenario_info, f<-function(x) x$inds*0 + x$mw)) 

    output = cbind(scenario_row=inds, mw=mw, scenario_rates=scenario_rates, 
                   weight_inflation=weight_inflation, alternate_scenario_rates=alternate_scenario_rates)

    write.csv(output, filename, row.names=FALSE)
}

write_scenario_info(random_scenario_info$unsegmented, 'random_scenarios_kermadectonga2_unsegmented_HS.csv')
write_scenario_info(random_scenario_info$tonga, 'random_scenarios_kermadectonga2_tonga_segment_HS.csv')
write_scenario_info(random_scenario_info$kermadec, 'random_scenarios_kermadectonga2_kermadec_segment_HS.csv')
write_scenario_info(random_scenario_info$hikurangi, 'random_scenarios_kermadectonga2_hukurangi_segment_HS.csv')


## With 800 KSU, and the model taking 1.5 hours on 96 CPU, I can run about this many scenarios:
#> 800*1000/(1.5 * 96 * 2)
#[1] 2777.778

random_musings_while_developing_the_above<-function(){

    # Get the HS scenarios
    kt_HS_scenarios = ptha18$get_source_zone_events_data(source_zone, slip_type='stochastic')
    stopifnot(nrow(kt_HS_scenarios$events) == length(kt_full$HS_prob_given_Mw))
    event_centroids = lapply(1:nrow(kt_HS_scenarios$events), 
        f<-function(i){
            get_event_slip_weighted_centroid(kt_HS_scenarios$events[i,], 
                                             kt_HS_scenarios$unit_source_statistics)
        })
    event_centroids = matrix(unlist(event_centroids), ncol=2, byrow=TRUE)

    LOWER_LON = -28
    UPPER_LON = -17
    # Group the scenarios based on their centroid latitude. The idea is that
    # we prioritize scenarios in the latitude range near Tonga, and decimate
    # them somehow outside of that.
    event_category = (event_centroids[,2] < LOWER_LON) + (event_centroids[,2] < UPPER_LON)

    event_is_possible = event_peak_stage_at_refpoint[[source_zone]]$scenario_rate_is_positive

    k = which(event_is_possible)
    plot(event_centroids[k,2], event_peak_stage[k], log='y', pch='.', cex=2, ylim=c(0.1, 20))
    abline(v=c(LOWER_LON, UPPER_LON))

    random_scenario_categories = lapply(random_scenario_inds, f<-function(x) event_category[x])
    mw_levels = c(7.45, 7.95, 8.45, 8.95, 9.45, 9.95)
    par(mfrow=c(length(mw_levels)-1, 1))
    par(mar=c(2,2,2,1))
    for(mm in 1:(length(mw_levels)-1)){
        lat_levels = seq(-44, -14, by=2)
        stages = c(0.25, 0.5, 1, 2, 3, 4, 5)
        count_big = matrix(NA, nrow=length(lat_levels), ncol=length(stages))
        si = unlist(random_scenario_inds)
        ss = unlist(random_scenario_stages)
        sc = unlist(random_scenario_categories)
        for(j in 1:length(stages)){
            for(i in 1:(length(lat_levels)-1)){
                count_big[i,j] = sum(ss > stages[j] & 
                                   event_centroids[si,2] > lat_levels[i] & 
                                   event_centroids[si,2] <= lat_levels[i+1] &
                                   kt_full$HS_mw[si] < mw_levels[mm+1] &
                                   kt_full$HS_mw[si] > mw_levels[mm])
            }
        }
        for(j in seq(1, length(stages))){
           if(j == 1){
              plot(lat_levels+diff(lat_levels)[1]/2, count_big[,j], t='o', col=j)
           }else{
              points(lat_levels+diff(lat_levels)[1]/2, count_big[,j], t='o', col=j)
           }
        }
        abline(v=c(LOWER_LON, UPPER_LON))
        title('Consider more random scenarios at higher Mw, as well as the spatial weighting. Also consider Mw-varying spatial dependence. \n Consider plots like this based on scenario rates [all scenarios]. See how well approaches reproduce stage-vs-rate at target gauge, and nearby gauges')
        legend('topleft', legend=as.character(stages), col=1:length(stages), lty='solid', pch=1,
               title='No. scenarios > stage @ target point')
    }
}
