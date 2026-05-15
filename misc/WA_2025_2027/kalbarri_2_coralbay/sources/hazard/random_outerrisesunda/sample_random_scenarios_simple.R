#
# Sample scenarios on sunda2. 
#
# Note we use VARIABLE RIGIDITY but retain binning by the "constant rigidity magnitude".
#

target_point = c(112, -24.333333) # Known location of PTHA18 hazard point
mu_type = 'constant' #'variable_mu' # 'constant' # Shear modulus model 
REPRODUCIBLE_SEED = 98765
SAMPLING_MW_LIMIT = 9.05 # Make this in-between Mm values in PTHA18 (7.2, 7.3, ....9.6, 9.7, 9.8)
TOTAL_SAMPLES = 150 # 9 scenarios x 19 bins that are below Mw limit
SAMPLES_PER_BIN_BELOW_MW_LIMIT = TOTAL_SAMPLES/19 # Naive sampling option
fraction_uniform_sampling = 0.4 # How much sampling effort to spend on uniform sampling (the rest spread over 2 non-uniform sampling efforts)
stage_seq = seq(0.1, 5, by=0.1) # Range for plotting
stage_threshold_1 = 0.2 # Optimize sampling effort partly for this threshold
stage_threshold_2 = 0.8 # Optimize sampling effort partly for this threshold
sampling_effort_plot_title= paste0('Optimal number of outerrisesunda thrust scenarios @ ', round(target_point[1], 3), ',', 
    round(target_point[2], 3), ' \n using logic-tree mean scenario-frequency model.')
# Read all heterogeneous-slip scenario metadata (slip_type='stochastic' in PTHA18)
source_zone = 'outerrisesunda'
stored_scenarios_file = 'outerrisesunda_scenarios_store.RDS'

# Event rates depend on the shear modulus model.
# Note that in PTHA18, earthquake Mw-frequency is ALWAYS a function of the constant shear modulus magnitude,
# although the function differs depending on the shear modulus.
if(mu_type == 'variable_mu'){
    event_rate_variable = 'variable_mu_rate_annual'
}else if(mu_type == 'constant'){
    event_rate_variable = 'rate_annual'
}else{
    stop(paste0('unknown value of mu_type ', mu_type))
}

# Get the scripts to access the PTHA18
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
file_here = ifelse(file.exists(file_nci), file_nci, file_home)
source(file_here, local=ptha18, chdir=TRUE)

if(!file.exists(stored_scenarios_file)){
    source_zone_scenarios = ptha18$get_source_zone_events_data(source_zone,  slip_type='stochastic')
    saveRDS(source_zone_scenarios, stored_scenarios_file)
}else{
    source_zone_scenarios = readRDS(stored_scenarios_file)
}

# Get the peak-stage at our reference point, offshore of Perth in about 100m depth, slightly north of west
event_peak_stage_at_refpoint = ptha18$get_peak_stage_at_point_for_each_event(
    target_point = target_point, 
    slip_type='stochastic', # Matches the slip-type used to create scenarios
    all_source_names=source_zone)

# Convenient shorthand
event_peak_stage_ref = event_peak_stage_at_refpoint[[source_zone]]$max_stage

# Convenient shorthand for the magnitudes and rates in the event table
event_Mw = source_zone_scenarios$events$Mw 
#event_rates = source_zone_scenarios$events$rate_annual # Logic-tree mean model
event_rates = source_zone_scenarios$events[[event_rate_variable]] # Logic-tree mean model

# Define the event importance
event_importance = event_peak_stage_ref # Could be another site -- suggest "near site of interest"
event_importance_weighted_sampling_probs = (event_rates * event_importance) # Importance sampling

stage_exrates_ptha18 = sapply(stage_seq, function(x) sum(event_rates*(event_peak_stage_ref > x)))

set.seed(REPRODUCIBLE_SEED)

#
# Importance sampling with uniform sampling of magnitude-bins (constant shear modulus magnitude)
#

# Random scenarios -- importance sampling.
random_scenarios_stage_weighted = ptha18$randomly_sample_scenarios_by_Mw_and_rate(
    event_rates=event_rates,
    event_Mw=event_Mw,
    event_importance_weighted_sampling_probs = (event_rates * event_importance),
    samples_per_Mw=function(Mw){ SAMPLES_PER_BIN_BELOW_MW_LIMIT*(Mw < SAMPLING_MW_LIMIT) }
    )

# Compute the max-stage exceedance-rates at the target_site
stage_exrates_rs_stage_weighted = sapply(stage_seq, 
    function(x){
        sum(random_scenarios_stage_weighted$importance_sampling_scenario_rates_basic * 
            (event_peak_stage_ref[random_scenarios_stage_weighted$inds] > x), na.rm=TRUE)
    })

#
# Importance-sampling with non-uniform sampling of magnitude-bins
#
unique_Mws = ptha18$unique_sorted_with_check_for_even_spacing(event_Mw)

# stage_threshold = 1
non_uniform_samples_IS_1 = ptha18$get_optimal_number_of_samples_per_Mw(
    event_Mw, event_rates, event_peak_stage_ref, stage_threshold=stage_threshold_1, 
    event_importance_weighted_sampling_probs = (event_rates*event_importance), # Importance sampling
    total_samples=TOTAL_SAMPLES)

# stage_threshold = 2
non_uniform_samples_IS_2 = ptha18$get_optimal_number_of_samples_per_Mw(
    event_Mw, event_rates, event_peak_stage_ref, stage_threshold=stage_threshold_2, 
    event_importance_weighted_sampling_probs = (event_rates*event_importance), # Importance sampling
    total_samples=TOTAL_SAMPLES)

# Uniform sampling up to max magnitude
uniform_sampling_effort = SAMPLES_PER_BIN_BELOW_MW_LIMIT * (unique_Mws < SAMPLING_MW_LIMIT) 

# Weight the optimal solutions with different stage-thresholds
average_nonuniform_sampling_effort_IS = 
    0.5*(non_uniform_samples_IS_1$Nsamples + non_uniform_samples_IS_2$Nsamples)

# 25% uniform, 75% weighted non-uniform, rounded to integer values
chosen_sampling_effort_IS = round(
    fraction_uniform_sampling*uniform_sampling_effort + 
    (1-fraction_uniform_sampling)*(average_nonuniform_sampling_effort_IS))
samples_per_Mw_stratified_importance = approxfun(unique_Mws, chosen_sampling_effort_IS, method='constant')

write.csv(data.frame(Mws = unique_Mws, sampling_effort=chosen_sampling_effort_IS), 
    file='Non_uniform_sampling_effort_compromise_stratifiedImportance.csv', 
    row.names=FALSE)

#### PLOT
png(paste0('Number_of_samples_', source_zone,'.png'), width=9, height=5, units='in', res=300)
plot(non_uniform_samples_IS_2[,1]+0.025, non_uniform_samples_IS_2[,2], type='h', lend=1, lwd=3, col='red4',
    xlab='Magnitude bin', ylab='Number of scenarios', cex.lab=1.5, cex.axis=1.4, cex.main=1.5,
    main=sampling_effort_plot_title)
points(non_uniform_samples_IS_1[,1]+0.01    , non_uniform_samples_IS_1[,2], type='h', lend=1, lwd=3, col='red3')
points(non_uniform_samples_IS_1[,1]-0.0, uniform_sampling_effort, type='h', lend=1, lwd=2, col='orange')
grid(col='black')
legend('topleft', 
    c('Uniform sampling', paste0('Threshold = ', stage_threshold_1, ' m'), paste0('Threshold = ', stage_threshold_2, ' m')), 
    lwd=c(2, 3, 3), col=c('orange', 'red3', 'red4'), pch=c(NA, NA, NA), cex=1.2)
dev.off()

# Make the random scenarios -- 
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

# Plot importance sampling without focus on particular magnitudes
png(paste0('Stage_vs_exrate_', source_zone, '.png'), width=9, height=8, units='in', res=300)
options(scipen=5)
plot(stage_seq, stage_exrates_ptha18, log='xy', xlab='Stage (m)', ylab='Exceedance rate (events/year)', cex.axis=1.4, cex.lab=1.4)
points(stage_seq, stage_exrates_rs_stage_weighted, t='l', col='red')
points(stage_seq, stage_exrates_rs_stage_mw_weighted, t='l', col='green')
grid()
points(stage_seq, mean_exrate, t='l', col='orange', lty='dashed')
points(stage_seq, mean_exrate + qnorm(0.025)*sqrt(var_exrate), t='l', col='orange', lty='dashed')
points(stage_seq, mean_exrate + qnorm(0.975)*sqrt(var_exrate), t='l', col='orange', lty='dashed')
legend('bottomleft', c('PTHA18', 'Even sampling Mw bins', 'Uneven sampling Mw bins',  '2sd Monte Carlo CI (uneven sampling Mw bins)'),
    col=c('orange', 'red', 'green', 'orange'), lty=c(rep('solid', 3), 'dashed'), bty='n', cex=1.2)
dev.off()
