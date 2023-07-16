#
# Pick random scenarios for a source zone with an importance-sampling approach.
#

#
# Get the PTHA18 access codes, needed for functions below
#
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
ptha18 = new.env()
if(!file.exists(file_here)){
    stop("You need to change file_home or file_nci to point to get_PTHA_results.R in the ptha/ptha_access directory")
}
source(file_here, local=ptha18, chdir=TRUE)

#
# Get the R session resulting from "compute_rates_all_sources.R", needed for functions below
#
ptha18_source_rate_env = new.env()
source(paste0(dirname(file_here), '/get_detailed_PTHA18_source_zone_info.R'),
       local=ptha18_source_rate_env, chdir=TRUE)

# Get the scenario rates for each model.
source_zone = 'outerrisesunda'

REPRODUCIBLE_SEED = 1234

sz_full      = ptha18_source_rate_env$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(source_zone, '')

# Weighted sum of 'unsegmented' and 'union of segments' rates as used in PTHA18
back_calculated_HS_rates_combination = 1.0*(sz_full$HS_event_rates) 

# Get the sampling effort in each magnitude-bin, as a function
samples_per_Mw_data = read.csv('Non_uniform_sampling_effort_compromise_stratifiedImportance.csv')
samples_per_Mw_FUNCTION = approxfun(samples_per_Mw_data[,1], ceiling(samples_per_Mw_data[,2]), method='constant')

check_consistency_with_PTHA18<-function(){
    # Double check consistency with PTHA18 files.

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


# Get the peak-stage at a point near Perth -- we will use this to measure the 
# scenario 'size' and influence sampling.
event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 1775.1,
    slip_type='stochastic',
    all_source_names=source_zone)
# Sampling will be done with the following weights (but they are stratified by
# Mw and converted to a pmf within the random sampling function)
EVENT_IMPORTANCE_WEIGHTED_SAMPLING_PROBS = 
    (event_peak_stage_at_refpoint[[source_zone]]$max_stage * 
     back_calculated_HS_rates_combination)

# We will use these nearby points for testing -- can we reproduce the stage-vs-rate curve
# with the randomly sampled scenarios? 
alternative1_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 1772.1, # Near Mandurah
    slip_type='stochastic',
    all_source_names=source_zone)
alternative2_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 3754.3, # Near Jurian Bay
    slip_type='stochastic',
    all_source_names=source_zone)
alternative3_event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    hazard_point_gaugeID = 1790.1, # Near Geraldton
    slip_type='stochastic',
    all_source_names=source_zone)


#
# Make a plot of max-stage vs exceedance-rate according to the random sample of scenarios,
# and compare with the full set of scenarios. 
#
# Aim is to reproduce the full results well, with many fewer scenarios.
# 
plot_scenario_stages_with_different_scenario_size_weightings<-function(
    event_peak_stage, 
    event_Mw, 
    event_rates,
    event_importance_weighted_sampling_probs = EVENT_IMPORTANCE_WEIGHTED_SAMPLING_PROBS,
    samples_per_Mw = samples_per_Mw_FUNCTION,
    event_peak_stage_alt1 = alternative1_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    event_peak_stage_alt2 = alternative2_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    event_peak_stage_alt3 = alternative3_event_peak_stage_at_refpoint[[source_zone]]$max_stage,
    pdf_output_name_extra = ''
    ){

    pdf(paste0('Effect_of_weighting_power_on_samples_fig_', pdf_output_name_extra, '.pdf'), 
        width=12, height=9)

    # Range of scenarios in PTHA18
    max_mw_limit = 9.05
    min_mw_limit = 7.15

    set.seed(REPRODUCIBLE_SEED) # Reproducible random seed

    random_scenario_info = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
        event_rates,
        event_Mw,
        event_importance_weighted_sampling_probs = event_importance_weighted_sampling_probs,
        samples_per_Mw = samples_per_Mw,
        mw_limits = c(min_mw_limit, max_mw_limit))

    random_scenario_inds = random_scenario_info$inds
    # Look at the stage values at the gauge offshore of Perth
    random_scenario_stages = event_peak_stage[random_scenario_inds]
    random_scenario_rates = random_scenario_info$importance_sampling_scenario_rates_basic
    random_scenario_Mws = random_scenario_info$mw

    #
    # Check if we can reproduce max-stage-vs-exceedance-rate curves
    #
    exrate_plot<-function(stage, rate, add=FALSE, ...){

        all_stages = seq(min(stage)-0.01, max(stage)-0.01, len=10000)
        all_rates = sapply(all_stages, function(x) sum(rate*(stage>x)))

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
                main=paste0('Reference point; Nunique=', 
                            length(unique(random_scenario_inds))))
    exrate_plot(random_scenario_stages, random_scenario_rates, 
                add=TRUE, col='green', t='l')

    # Check a nearby site.
    random_scenario_stages_alt1 = event_peak_stage_alt1[random_scenario_inds]
    exrate_plot(event_peak_stage_alt1, event_rates, t='l', main='Alternative 1')
    exrate_plot(random_scenario_stages_alt1, random_scenario_rates, 
                add=TRUE, col='green', t='l')
    # An another nearby site
    random_scenario_stages_alt2 = event_peak_stage_alt2[random_scenario_inds]
    exrate_plot(event_peak_stage_alt2, event_rates, t='l', main='Alternative 2')
    exrate_plot(random_scenario_stages_alt2, random_scenario_rates, 
                add=TRUE, col='green', t='l')
    # An another nearby site
    random_scenario_stages_alt3 = event_peak_stage_alt3[random_scenario_inds]
    exrate_plot(event_peak_stage_alt3, event_rates, t='l', main='Alternative 3')
    exrate_plot(random_scenario_stages_alt3, random_scenario_rates, 
                add=TRUE, col='green', t='l')

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
        k = which(random_scenario_stages > stages[i])
        if(length(k) > 0){
            mws = weighted_percentile(random_scenario_Mws[k], 
                random_scenario_rates[k], percentiles)
            points(percentiles, mws, t='l', col='red')
        }else{
            title(sub='no random scenarios were found', col.sub='red')
        }

    }
    dev.off()

    return(invisible(random_scenario_inds))
}

#
# Check how well we can reproduce the full 'max-stage vs exceedance-rate' curve
# on segmented/unsegmented branches respectively. This code will make plots.
#
event_peak_stage = event_peak_stage_at_refpoint[[source_zone]]$max_stage
event_Mw = event_peak_stage_at_refpoint[[source_zone]]$Mw
inds_for_checking = plot_scenario_stages_with_different_scenario_size_weightings(
    event_peak_stage, 
    event_Mw, 
    sz_full$HS_event_rates,
    pdf_output_name_extra='unsegmented')
#
# Make the samples, with more scenarios at higher magnitudes, and a "scenario_size"
# proportional to the event_peak_stage.
#
sampling_setup = list(
    'unsegmented' = list(event_rates = sz_full$HS_event_rates)
    )

random_scenario_info = vector(mode='list', length=length(sampling_setup))
names(random_scenario_info) = names(sampling_setup)
for(i in 1:length(random_scenario_info)){
    # Reproducible random seed -- note we DELIBERATELY SAMPLE THE SAME
    # SCENARIOS ON UNSEGMENTED/SEGMENT.
    set.seed(REPRODUCIBLE_SEED)

    # Here are the magnitude limits where we bother sampling
    max_mw_limit = 9.05
    min_mw_limit = 7.15

    # Take the sample
    random_scenario_info[[i]] = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
        sampling_setup[[i]]$event_rates,
        event_Mw,
        event_importance_weighted_sampling_probs=EVENT_IMPORTANCE_WEIGHTED_SAMPLING_PROBS,
        samples_per_Mw = samples_per_Mw_FUNCTION,
        mw_limits = c(min_mw_limit, max_mw_limit))

    # Note we DELIBERATELY SAMPLE THE SAME SCENARIOS ON UNSEGMENTED/SEGMENT.
    # Confirm that we did
    stopifnot(all(random_scenario_info[[i]]$inds == inds_for_checking))

}

# How many scenarios will we have to run?
length(unique(random_scenario_info[[1]]$inds))

# Check we got the expected number of scenarios in each Mw bin 
sampled_per_mw = aggregate(random_scenario_info[[1]]$inds, 
    by=list(random_scenario_info[[1]]$mw), function(x){length(x)})$x
NN = length(sampled_per_mw) # Skip Mw 9.7, 9.8
stopifnot(all(sampled_per_mw == ceiling(samples_per_Mw_data[1:NN,2])))

# How many unique scenarios in each Mw bin?
unique_sampled_per_mw = aggregate(random_scenario_info[[1]]$inds, 
    by=list(random_scenario_info[[1]]$mw), function(x){length(unique(x))})$x
unique_fraction = unique_sampled_per_mw/sampled_per_mw

##
## Write out the scenarios to a csv format
##

write.csv(random_scenario_info$unsegmented, 'random_scenarios_outerrisesunda_unsegmented_HS.csv', row.names=FALSE)

#
# Convenient to also store the logic-tree-mean results -- as above, reproducible randomness
#
set.seed(REPRODUCIBLE_SEED)
random_scenarios_mean_curve = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    back_calculated_HS_rates_combination,
    event_Mw,
    event_importance_weighted_sampling_probs=EVENT_IMPORTANCE_WEIGHTED_SAMPLING_PROBS,
    samples_per_Mw = samples_per_Mw_FUNCTION,
    mw_limits = c(min_mw_limit, max_mw_limit))
stopifnot(all(random_scenarios_mean_curve$inds == inds_for_checking))
write.csv(random_scenarios_mean_curve, 
          'random_scenarios_outerrisesunda_logic_tree_mean_curve_HS.csv', row.names=FALSE)
