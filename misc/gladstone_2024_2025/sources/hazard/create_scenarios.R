#
# Sample the source-zones with details specified in sampling_config.R
#

# Configuration for these samples
sc = new.env()
source('sampling_config.R', local=sc)

# For getting PTHA18 results
ptha18 = new.env()
source(sc$get_ptha_results_script_file, local=ptha18, chdir=TRUE)

# For getting detailed PTHA18 results
#ptha18_detailed = new.env()
#source(sc$get_detailed_ptha18_source_zone_info_script_file, local=ptha18_detailed, chdir=TRUE)

# Importance sampling calculations
isu = new.env()
source('../importance_sampling_utilities.R', local=isu)

# Names of the source zones
source_zones = names(sc$source_zones_and_samples)

# Function to sample scenarios on a given source_zone
sample_source<-function(source_zone){

    # Reproducible randomness
    set.seed(sc$REPRODUCIBLE_SEED)

    # Number of scenarios to be sampled
    N_MC = sc$source_zones_and_samples[[source_zone]]

    print(source_zone)

    # For this source zone, get peak stage at the target point
    event_peak_stage_at_targetpt = ptha18$get_peak_stage_at_point_for_each_event(
        target_point = sc$target_pt, # Known location of PTHA18 hazard point
        slip_type='stochastic', # Constant rigidity HS
        all_source_names=source_zone)
    peak_stage_target_pt = event_peak_stage_at_targetpt[[source_zone]]$max_stage

    # For this source zone, get the event table
    source_zone_scenarios = ptha18$get_source_zone_events_data(source_zone,
        slip_type='stochastic')

    # Use the logic-tree-mean rate and importance-factor I(e) to determine
    # the chance of sampling each scenario
    event_rate_logic_tree_mean = source_zone_scenarios$events$rate_annual
    # I(e) 
    importance_factor = sc$importance_function(
        peak_stage_target_pt,
        event_rate_logic_tree_mean,
        sc$stage_vs_rate_all_source_zones_at_target_pt)

    # Chance of sampling each scenario -- this is designed to well represent
    # events of interest near sc$target_pt
    sampling_prob = event_rate_logic_tree_mean * importance_factor
    sampling_prob = sampling_prob/sum(sampling_prob)

    # Compute the analytical exceedance-rate and its MC variance for the logic-tree mean
    stage_thresholds = sc$stage_thresholds
    analytical_mc_results = lapply(stage_thresholds, function(x){
        isu$analytical_exrate_and_MC_variance(
            event_rate=event_rate_logic_tree_mean, 
            event_stage=peak_stage_target_pt, 
            stage_threshold=x, 
            sampling_prob=sampling_prob, 
            N_MC=N_MC)})
    analytical_exrate = unlist(lapply(analytical_mc_results, function(x) x[1]))     
    mc_exrate_variance = unlist(lapply(analytical_mc_results, function(x) x[2]))

    # Get a couple of MC samples
    test_sample = isu$sample_scenarios_monte_carlo(
        event_rate=event_rate_logic_tree_mean, 
        event_stage=peak_stage_target_pt, 
        stage_thresholds=stage_thresholds, 
        sampling_prob=sampling_prob, 
        N_MC=N_MC)
    alt_test_sample = isu$sample_scenarios_monte_carlo(
        event_rate=event_rate_logic_tree_mean, 
        event_stage=peak_stage_target_pt, 
        stage_thresholds=stage_thresholds, 
        sampling_prob=sampling_prob, 
        N_MC=N_MC)

    # Get some stage data at other points
    peak_stage_other_pts = lapply(sc$other_test_points, function(pt){
       ptha18$get_peak_stage_at_point_for_each_event(
        target_point = pt, # Known location of PTHA18 hazard point
        slip_type='stochastic', # Constant rigidity HS
        all_source_names=source_zone)}
       )

    ##
    ## Prepare data for unsegmented and segmented source representations
    ## We don't need this here
    ##
    ## Unsegmented index
    #unseg_ind = which(names(ptha18_detailed$crs_data$source_envs) == source_zone)
    ## Segmented index -- they start with "sourcezonename_"
    #seg_ind = grep(paste0(source_zone, '_'), names(ptha18_detailed$crs_data$source_envs))
    #stopifnot( (length(unseg_ind) == 1) ) # Only 1 segmented source
    #if(length(seg_ind) > 0){
    #    # Segmented sources come straight after unsegmented sources, and in a sequence
    #    stopifnot( (seg_ind[1] == unseg_ind[1] + 1) & all(range(diff(seg_ind)) == 1) )
    #}
    #all_source_names = names(ptha18_detailed$crs_data$source_envs)[c(unseg_ind, seg_ind)]
    #all_source_segment_names = gsub('^_', '', gsub(paste0(source_zone), '', all_source_names))
    #all_source_unseg_or_seg = c('unsegmented', rep('segment', length(all_source_names)-1))


    # Save the sampled scenarios for later usage
    # Compared to our old importance sampling method, with this method the calculations of epistemic
    # uncertainties will need to follow a different approach, and we don't setup those calculations as much here.
    output_sample = data.frame(
        inds = test_sample$si, # Row indices
        event_Mw = source_zone_scenarios$events$Mw[test_sample$si], # Magnitude
        sampling_prob=sampling_prob[test_sample$si], # sampling probabilities applied to the sampled scenarios
        importance_factor=importance_factor[test_sample$si], # I(e) for sampled scenarios
        event_rate_logic_tree_mean = event_rate_logic_tree_mean[test_sample$si], # Logic-tree-mean rate for sampled scenarios
        basic_importance_sampling_weights_logic_tree_mean = test_sample$sampled_bis_wts, # Basic importance sampling weights for samlped scenarios
        TOTAL_RATE_logic_tree_mean = rep(test_sample$TOTAL_RATE, length=length(test_sample$si)),
        basic_importance_sampling_rates_logic_tree_mean = 
            test_sample$TOTAL_RATE * test_sample$sampled_bis_wts / length(test_sample$sampled_bis_wts))

    output_dir = paste0('random_', source_zone)
    dir.create(output_dir, showWarnings=FALSE)

    write.csv(output_sample, 
        file=paste0(output_dir, '/random_scenarios_', source_zone, '_logic_tree_mean_HS.csv'), 
        row.names=FALSE)

    return(environment())
}

source_zone_data = lapply(source_zones, sample_source)
names(source_zone_data) = source_zones

saveRDS(source_zone_data, file='source_zone_data_backup.RDS')

# Plot
pdf('source_zone_performance.pdf', width=8, height=6)

for(i in 1:length(source_zone_data)){
    source_zone = names(source_zone_data)[i]
    stage_thresholds = source_zone_data[[i]]$stage_thresholds
    analytical_exrate = source_zone_data[[i]]$analytical_exrate
    mc_exrate_variance = source_zone_data[[i]]$mc_exrate_variance
    test_sample_exrates = source_zone_data[[i]]$test_sample$exrates
    alt_test_sample_exrates = source_zone_data[[i]]$alt_test_sample$exrates

    plot(stage_thresholds, analytical_exrate, ylim=c(1e-06, 1), log='xy', t='l')
    grid(); add_log_axis_ticks(1); add_log_axis_ticks(2)
    ci_offset = 1.96*sqrt(mc_exrate_variance)
    points(stage_thresholds, analytical_exrate + ci_offset, t='l', col='red')
    points(stage_thresholds, analytical_exrate - ci_offset, t='l', col='red')
    points(stage_thresholds, test_sample_exrates, t='l', col='green', lty='dashed')
    points(stage_thresholds, alt_test_sample_exrates, t='l', col='purple', lty='dashed')
    title(paste0(source_zone, ': N_MC=', sc$source_zones_and_samples[[source_zone]]))
    legend('topright', legend=c('Analytical', '95% CI', 'MC sample', 'MC sample alternative'), 
           col=c('black', 'red', 'green', 'purple'), lty=c('solid', 'solid', 'dashed', 'dashed'))
}

combined_exrate = source_zone_data[[1]]$analytical_exrate*0
combined_mc_exrate_variance = source_zone_data[[1]]$mc_exrate_variance*0
combined_sample = source_zone_data[[1]]$test_sample$exrates*0
combined_sample2 = source_zone_data[[1]]$alt_test_sample$exrates*0
for(i in 1:length(source_zone_data)){
    combined_exrate = combined_exrate + source_zone_data[[i]]$analytical_exrate
    combined_mc_exrate_variance = combined_mc_exrate_variance + source_zone_data[[i]]$mc_exrate_variance
    combined_sample = combined_sample + source_zone_data[[i]]$test_sample$exrates
    combined_sample2 = combined_sample2 + source_zone_data[[i]]$alt_test_sample$exrates
}
plot(stage_thresholds, combined_exrate, ylim=c(1e-06, 1), log='xy', t='l')
grid(); add_log_axis_ticks(1); add_log_axis_ticks(2)
ci_offset = 1.96*sqrt(combined_mc_exrate_variance)
points(stage_thresholds, combined_exrate - ci_offset, col='red', t='l')
points(stage_thresholds, combined_exrate + ci_offset, col='red', t='l')
points(stage_thresholds, combined_sample, t='l', col='green', lty='dashed')
points(stage_thresholds, combined_sample2, t='l', col='purple', lty='dashed')
title('Sum over source-zones')

#
# PTHA18 all-source-zones result vs the above sum
#
plot(stage_thresholds, combined_exrate, ylim=c(1e-06, 1), log='xy', t='l')
grid(); add_log_axis_ticks(1); add_log_axis_ticks(2)
points(sc$stage_vs_rate_all_source_zones_at_target_pt$stage, 
       sc$stage_vs_rate_all_source_zones_at_target_pt$stochastic_slip_rate,
       t='l', col='blue')
points(stage_thresholds, combined_sample, t='l', col='green', lty='dashed')
points(stage_thresholds, combined_sample2, t='l', col='purple', lty='dashed')
title('Sum over source-zones: Blue curve gives PTHA18 all source-zones result \n (interpolated from different absicca)')

dev.off()

# Summary of sampled stages (for the single MC sample)
stage_sam = lapply(source_zone_data, function(x){
    sort(x$test_sample$sampled_stages) })

# Summary of sampled Mw (for the single MC sample)
Mw_sam = lapply(source_zone_data, function(x){ 
    sort(x$source_zone_scenarios$events$Mw[x$test_sample$si]) })

unique_scenarios = lapply(source_zone_data, function(x) length(unique(x$test_sample$si)))

#
# Plot performance at other spatial locations
#
other_sites_store = list()
pdf('source_zone_performance_other_points.pdf', width=8, height=6)
for(i in 1:length(source_zone_data)){

    source_zone = names(source_zone_data)[i]
    szdi = source_zone_data[[i]]

    other_sites_store[[i]] = vector(mode='list', length=length(sc$other_test_points))

    par(mfrow=c(2, ceiling(length(sc$other_test_points)/2)))
    for(j in 1:length(sc$other_test_points)){

        # Get the PTHA18 results at the 'other test point'
        site_name = names(sc$other_test_points)[j]
        event_rate = szdi$event_rate_logic_tree_mean
        peak_stage_pt_j = szdi$peak_stage_other_pts[[j]][[source_zone]]$max_stage
        sampling_prob = szdi$sampling_prob
        N_MC = szdi$N_MC
        stage_thresholds = szdi$stage_thresholds

        analytical_mc_results = lapply(stage_thresholds, function(stg){
            isu$analytical_exrate_and_MC_variance(
                event_rate = event_rate,
                event_stage = peak_stage_pt_j,
                stage_threshold=stg,
                sampling_prob=sampling_prob,
                N_MC=N_MC) })
        analytical_exrate = unlist(lapply(analytical_mc_results, function(x) x[1]))
        mc_exrate_variance = unlist(lapply(analytical_mc_results, function(x) x[2]))

        # Get the test_sample exceedance-rates at the test point
        test_sample = isu$sample_scenarios_monte_carlo(
             event_rate = event_rate,
             event_stage = peak_stage_pt_j,
             stage_thresholds=stage_thresholds,
             sampling_prob=sampling_prob,
             N_MC=N_MC,
             myseed = szdi$test_sample$myseed) # Ensure same random sample as previous -- only stage changes

        plot(stage_thresholds, analytical_exrate, ylim=c(1e-06, 1), log='xy', t='l')
        grid(); add_log_axis_ticks(1); add_log_axis_ticks(2)
        ci_offset = 1.96*sqrt(mc_exrate_variance)
        points(stage_thresholds, analytical_exrate + ci_offset, t='l', col='red')
        points(stage_thresholds, analytical_exrate - ci_offset, t='l', col='red')
        points(stage_thresholds, test_sample$exrates, t='l', col='green', lty='dashed')
        title(paste0(source_zone, ': N_MC=', sc$source_zones_and_samples[[source_zone]], ' \n', site_name))

        other_sites_store[[i]][[j]] = list(
            site_name = site_name,
            source_zone = source_zone,
            stage_thresholds = stage_thresholds,
            analytical_exrate = analytical_exrate,
            mc_exrate_variance = mc_exrate_variance,
            test_sample_exrates = test_sample$exrates)
    }
}

# Combined source plot
par(mfrow=c(2, ceiling(length(sc$other_test_points)/2)))

stage_thresholds = other_sites_store[[1]][[1]]$stage_thresholds
for(j in 1:length(sc$other_test_points)){
    site_name = names(sc$other_test_points)[j]
    analytical_exrate = colSums(do.call(rbind, lapply(other_sites_store, function(x) x[[j]]$analytical_exrate)))
    mc_exrate_variance = colSums(do.call(rbind, lapply(other_sites_store, function(x) x[[j]]$mc_exrate_variance)))
    test_sample_exrates = colSums(do.call(rbind, lapply(other_sites_store, function(x) x[[j]]$test_sample_exrates)))

    plot(stage_thresholds, analytical_exrate, ylim=c(1e-06, 1), log='xy', t='l')
    grid(); add_log_axis_ticks(1); add_log_axis_ticks(2)
    ci_offset = 1.96*sqrt(mc_exrate_variance)
    points(stage_thresholds, analytical_exrate + ci_offset, t='l', col='red')
    points(stage_thresholds, analytical_exrate - ci_offset, t='l', col='red')
    points(stage_thresholds, test_sample_exrates, t='l', col='green', lty='dashed')
    # Add "all source zones" solution from PTHA18 in blue
    points(sc$stage_vs_rate_all_source_zones_other_test_pts[[j]]$stage,
           sc$stage_vs_rate_all_source_zones_other_test_pts[[j]]$stochastic_slip_rate,
           t='l', col='blue')
    title(paste0('Combined sources (& PTHA18 in blue): \n', site_name))
}

dev.off()

