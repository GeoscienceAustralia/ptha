#
# Sample scenarios on outerrisesunda. This script mirrors this tutorial:
#    - https://github.com/GeoscienceAustralia/ptha/blob/master/ptha_access/example_event_access_scripts/random_scenarios_non_uniform_and_importance_sampling/random_scenario_sampling.md
#

# Get the scripts to access the PTHA18
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
source(file_here, local=ptha18, chdir=TRUE)

REPRODUCIBLE_SEED = 12345
TOTAL_SAMPLES = 75
SAMPLES_PER_BIN = TOTAL_SAMPLES/25

# Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in PTHA18)
source_zone = 'outerrisesunda'
stored_scenarios_file = 'outerrisesunda_scenarios_store.RDS'
if(!file.exists(stored_scenarios_file)){
    source_zone_scenarios = ptha18$get_source_zone_events_data(source_zone,  slip_type='stochastic')
    saveRDS(source_zone_scenarios, stored_scenarios_file)
}else{
    source_zone_scenarios = readRDS(stored_scenarios_file)
}

# Get the peak-stage at our reference point, offshore of Perth in about 100m depth, slightly north of west
event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    target_point = c(115.234169, -31.6791019), # Known location of PTHA18 hazard point
    slip_type='stochastic', # Matches the slip-type used to create scenarios
    all_source_names=source_zone)

# Convenient shorthand
event_peak_stage_ref = event_peak_stage_at_refpoint[[source_zone]]$max_stage

# Convenient shorthand for the magnitudes and rates in the event table
event_Mw = source_zone_scenarios$events$Mw 
event_rates = source_zone_scenarios$events$rate_annual # Logic-tree mean model

# Define the event importance
event_importance = event_peak_stage_ref # Could be another site -- suggest "near site of interest"
event_importance_weighted_sampling_probs = (event_rates * event_importance) # Importance sampling

stage_seq = seq(0.1, 4, by=0.1)
stage_exrates_ptha18 = sapply(stage_seq, function(x) sum(event_rates*(event_peak_stage_ref > x)))

#
# Importance sampling with uniform sampling of magnitude-bins
#

set.seed(REPRODUCIBLE_SEED)
# Random scenarios -- importance sampling, 12 samples per magnitude-bin.
random_scenarios_stage_weighted = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    event_importance_weighted_sampling_probs = (event_rates * event_importance),
    samples_per_Mw=function(Mw){ SAMPLES_PER_BIN*(Mw < 9.05) }
    )

# Compute the max-stage exceedance-rates at the site offshore Perth
stage_exrates_rs_stage_weighted = sapply(stage_seq, 
    function(x){
        sum(random_scenarios_stage_weighted$importance_sampling_scenario_rates_basic * 
            (event_peak_stage_ref[random_scenarios_stage_weighted$inds] > x), na.rm=TRUE)
    })

# Plot importance sampling without focus on particular magnitudes
plot(stage_seq, stage_exrates_ptha18, log='xy')
points(stage_seq, stage_exrates_rs_stage_weighted, t='l', col='red')
grid()

#
# Importance-sampling with non-uniform sampling of magnitude-bins
#
unique_Mws = ptha18$unique_sorted_with_check_for_even_spacing(event_Mw)

# stage_threshold = 0.3
non_uniform_samples_IS_0p3 = ptha18$get_optimal_number_of_samples_per_Mw(
    event_Mw, event_rates, event_peak_stage_ref, stage_threshold=0.3, 
    event_importance_weighted_sampling_probs = (event_rates*event_importance), # Importance sampling
    total_samples=TOTAL_SAMPLES)

# stage_threshold = 0.75
non_uniform_samples_IS_0p75 = ptha18$get_optimal_number_of_samples_per_Mw(
    event_Mw, event_rates, event_peak_stage_ref, stage_threshold=0.75, 
    event_importance_weighted_sampling_probs = (event_rates*event_importance), # Importance sampling
    total_samples=TOTAL_SAMPLES)

# Uniform sampling up to magnitude 9.
uniform_sampling_effort = SAMPLES_PER_BIN * (unique_Mws < 9.05) 

# Weight the optimal solutions with different stage-thresholds
average_nonuniform_sampling_effort_IS = 
    0.5*(non_uniform_samples_IS_0p3$Nsamples + non_uniform_samples_IS_0p75$Nsamples)

# 25% uniform, 75% weighted non-uniform, rounded to integer values
chosen_sampling_effort_IS = round(
    0.25*uniform_sampling_effort + 
    0.75*(average_nonuniform_sampling_effort_IS))
samples_per_Mw_stratified_importance = approxfun(unique_Mws, chosen_sampling_effort_IS, method='constant')

write.csv(data.frame(Mws = unique_Mws, sampling_effort=chosen_sampling_effort_IS), 
    file='Non_uniform_sampling_effort_compromise_stratifiedImportance.csv', 
    row.names=FALSE)

# Make the random scenarios -- 
set.seed(REPRODUCIBLE_SEED)
random_scenarios_stage_mw_weighted = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    event_importance_weighted_sampling_probs = (event_rates * event_importance),
    samples_per_Mw=samples_per_Mw_stratified_importance
    )

# Compute the max-stage exceedance-rates
stage_exrates_rs_stage_mw_weighted = sapply(stage_seq,
    function(x){
        sum(random_scenarios_stage_mw_weighted$importance_sampling_scenario_rates_basic * 
            (event_peak_stage_ref[random_scenarios_stage_mw_weighted$inds] > x), na.rm=TRUE)
    })

points(stage_seq, stage_exrates_rs_stage_mw_weighted, t='l', col='green')



#
# Check the errors
#
mean_exrate = stage_seq * NA
var_exrate = stage_seq * NA
for(i in 1:length(stage_seq)){
    # Get the "exact" exceedance-rate and the analytical Monte-Carlo variance
    exrate_and_var_stratified_importance = ptha18$analytical_Monte_Carlo_exrate_uncertainty(
        event_Mw, event_rates, 
        event_peak_stage_ref, 
        stage_threshold=stage_seq[i], 
        samples_per_Mw=samples_per_Mw_stratified_importance,
        event_importance_weighted_sampling_probs = (event_rates * event_importance) # Importance sampling
        )

    mean_exrate[i] = exrate_and_var_stratified_importance[1]
    var_exrate[i] = exrate_and_var_stratified_importance[2]
}

points(stage_seq, mean_exrate, t='l', col='orange', lty='dashed')
points(stage_seq, mean_exrate + qnorm(0.025)*sqrt(var_exrate), t='l', col='orange', lty='dashed')
points(stage_seq, mean_exrate + qnorm(0.975)*sqrt(var_exrate), t='l', col='orange', lty='dashed')


