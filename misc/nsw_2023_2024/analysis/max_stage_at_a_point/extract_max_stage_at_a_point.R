#
# Plots comparing max-stage and exceedance-rates in PTHA18 (frictionless, linear solver, 36hrs)
# with nonlinear model (high-res, 24 hours) using importance sampling.
#

# SWALS routines (used to find the domain containing a point)
swals = new.env()
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R', local=swals)

# Importance sampling utilities
isu = new.env()
source('../../sources/hazard/importance_sampling_utilities.R', local=isu)

# PTHA18 reader
ptha18 = new.env()
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '../../../../../../../AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

# Misc functions that were previously inside this script
source('extract_max_stage_utilities.R')

#
# INPUTS
#

# Use a command-line integer argument to choose which point to plot
point_code = as.numeric(commandArgs(trailingOnly=TRUE)[1])

# Define the background sea-level in the nonlinear model (so we can compare nonlinear model
# wave-heights to PTHA18)
nonlinear_model_MSL = 1.1

# Working with tarred raster files
tarred_raster_files =
    Sys.glob('../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/*/*/raster_output_files.tar')
# Each raster_output_files.tar contains a set of rasters, including 'raster_filename'

# A multidomain directory with design that matches the one used for the
# importance-sampling model runs (in terms of the domain layout).
# Will be used to associate the target_point with a particular
# max_stage_domain_XX.tif raster.
reference_multidomain_dir =
    '../../swals/OUTPUTS/run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0/RUN_20231129_102535028/'

# How many cores (cannot trust auto-detect on NCI).
# The script isn't very parallel, but assuming we need 16GB memory, then on the Cascade Lake nodes we may as well use 
# the corresponding number of cores { = 16/(192/48) cores }
MC_CORES = 4

# Files with scenarios used for importance-sampling -- logic_tree_mean curve results
importance_sampling_scenarios_logic_tree_mean = list(
    alaskaaleutians = read.csv('../../sources/hazard/scenarios_ID710.5/random_alaskaaleutians/random_scenarios_alaskaaleutians_logic_tree_mean_HS.csv'),
    kermadectonga2 = read.csv('../../sources/hazard/scenarios_ID710.5/random_kermadectonga2/random_scenarios_kermadectonga2_logic_tree_mean_HS.csv'),
    newhebrides2 = read.csv('../../sources/hazard/scenarios_ID710.5/random_newhebrides2/random_scenarios_newhebrides2_logic_tree_mean_HS.csv'),
    outerrise_kermadectonga = read.csv('../../sources/hazard/scenarios_ID710.5/random_outerrise_kermadectonga/random_scenarios_outerrise_kermadectonga_logic_tree_mean_HS.csv'),
    outerrisenewhebrides = read.csv('../../sources/hazard/scenarios_ID710.5/random_outerrisenewhebrides/random_scenarios_outerrisenewhebrides_logic_tree_mean_HS.csv'),
    outerrise_puysegur = read.csv('../../sources/hazard/scenarios_ID710.5/random_outerrise_puysegur/random_scenarios_outerrise_puysegur_logic_tree_mean_HS.csv'),
    puysegur2 = read.csv('../../sources/hazard/scenarios_ID710.5/random_puysegur2/random_scenarios_puysegur2_logic_tree_mean_HS.csv'),
    solomon2 = read.csv('../../sources/hazard/scenarios_ID710.5/random_solomon2/random_scenarios_solomon2_logic_tree_mean_HS.csv'),
    southamerica = read.csv('../../sources/hazard/scenarios_ID710.5/random_southamerica/random_scenarios_southamerica_logic_tree_mean_HS.csv'))

# Points at which we compute the exceedance-rate curve for the logic-tree-mean calculations
STAGE_POINTS_FOR_EXRATE_CURVE = seq(0.01, 5.01, len=1001)

# Y-range for plots of exceedance-rate curves for the logic-tree-mean calculations
EXRATE_PLOT_YLIM = c(1.0e-05, 0.1)

# Turn on/off the epistemic uncertainty calculations.
DO_EPISTEMIC_UNCERTAINTY_CALCULATIONS = TRUE # FALSE

# Random sampling is used to compute percentiles from combined probability
# distributions -- how many random samples should be used?
epistemic_uncertainty_Nsamples = 1e+04
epistemic_uncertainty_random_seed = 123 # Reproducible randomness
# Epistemic uncertainty calculations are a bit expensive -- compute at fewer max-stage values
epistemic_uncertainty_threshold_stage_values = nonlinear_model_MSL +
    c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)

#
# Choose site
#
# FIXME: Might need to avoid points exactly on domain bin boundaries -- currently
# confusing the raster lookup
if(point_code == 1){
    # Well offshore of Ulladulla, 4850m water depth
    target_point =c(151.3333, -35.6666) # gauge ID 710.5 
}else if(point_code == 2){
    # Well offshore of Tweed river entrance, 4500 m water depth.
    target_point = c(154.6666, -28.33333)
}else if(point_code == 3){
    # Well offshore of Eden, ~ 4600m water depth
    target_point = c(151.33333, -37.66666)
}else if(point_code == 4){
    # Well offshore of Port Stephens, ~ 4700m water depth
    target_point = c(154, -32.66666)
}else{
    stop('unknown point code')
}


#
# END INPUTS
#

# Get the name of the max-stage raster file that contains the point
tmp = swals$find_domain_containing_point(target_point, multidomain_dir=reference_multidomain_dir)
raster_filename = paste0('max_stage_domain_', swals$domain_index_from_folder(basename(tmp$domain_dir)), '.tif'); rm(tmp)

# Get the cell indices on a raster file that correspond to the target_point
tmp = get_indices_of_target_point_on_raster_file(tarred_raster_files[1], raster_filename, target_point)
XIND = tmp$XIND
YIND=tmp$YIND
coordinate_errtol = tmp$coordinate_errtol
rm(tmp)

# Get all max-stage values from the nonlinear model runs.
# This only includes each unique scenario. Some scenarios may have been sampled twice,
# but the model is run only once.
library(parallel)
all_max_stage = mclapply(tarred_raster_files,
    get_max_stage_at_target_point,
    raster_filename=raster_filename,
    XIND=XIND, YIND=YIND,
    target_point=target_point,
    coordinate_errtol = coordinate_errtol,
    mc.cores=MC_CORES)
# Convert to a data.frame
run_stage_at_target_point = data.frame(file=tarred_raster_files, max_stage=unlist(all_max_stage))

if(any(is.na(run_stage_at_target_point$max_stage))){
    stop(paste0("Finding NA max-stage values. ",
        "This can happen if target_point is right on the boundary of two domains, ",
        "given how the raster lookup works here. The simplest workaround is to use a nearby ",
        "point that isn't right on the boundary. A more complex fix would be to make a ",
        "vrt file with multiple neighbouring tifs, and extract from that."))
}

rm(all_max_stage)
# Save for easier use later
output_RDS = paste0('run_stage_at_target_point_', target_point[1], '_', target_point[2], '.RDS')
saveRDS(run_stage_at_target_point, output_RDS)

#
# Extract relevant scenario information from the folder names
# This involves a few potentially fragile hacks (depending on how the folder names are defined)
#
bdr = basename(dirname(run_stage_at_target_point$file)) # Can extract info from this piece of file path
split_bdr = strsplit(bdr, '_') # Intermediate step for extracting important information
source_name = unlist(lapply(split_bdr, function(x){
    k0 = which(x == "scenarios")
    k1 = which(x == 'row')
    return(paste(x[(k0+1):(k1-1)], collapse="_"))
}))
scenario_row = unlist(lapply(split_bdr, function(x){
    k0 = which(x == 'row')
    return(as.numeric(x[k0+1]))
}))
scenario_Mw = unlist(lapply(split_bdr, function(x){
    k0 = which(x == 'Mw')
    return(as.numeric(x[k0+1])/10)
}))
# Do a few sanity checks (since a change in the folder name structure could break the above code)
stopifnot(all(is.finite(scenario_Mw) & scenario_Mw > 7 & scenario_Mw < 10))
stopifnot(all(unique(source_name) %in% ptha18$config_env$source_names_all))
stopifnot(all(is.finite(scenario_row) & scenario_row > 0 & scenario_row == round(scenario_row)))

#
# Get the PTHA18 max stage values and other metadata
#
all_max_stage_ptha18 = ptha18$get_peak_stage_at_point_for_each_event(
    target_point=target_point,
    all_source_names = unique(source_name))
# Make a convenient data structure with PTHA18 max-stage and other info.
ptha18_max_stage = rep(NA, nrow(run_stage_at_target_point))
ptha18_Mw = rep(NA, nrow(run_stage_at_target_point))
for(i in 1:nrow(run_stage_at_target_point)){
    ptha18_max_stage[i] = all_max_stage_ptha18[[source_name[i]]]$max_stage[scenario_row[i]]
    ptha18_Mw[i] = all_max_stage_ptha18[[source_name[i]]]$Mw[scenario_row[i]]
}
scenario_max_stages = cbind(run_stage_at_target_point,
    data.frame('ptha18_Mw' = ptha18_Mw, 'ptha18_max_stage' = ptha18_max_stage,
               'nonlinear_max_stage_minus_MSL' = run_stage_at_target_point$max_stage - nonlinear_model_MSL,
               'nonlinear_Mw' = scenario_Mw,
               'nonlinear_source_name' = source_name,
               'nonlinear_scenario_row' = scenario_row))
out_file = paste0('scenario_max_stages_', target_point[1], '_', target_point[2], '.RDS')
saveRDS(scenario_max_stages, out_file)

#
# Make some plots of the tsunami maxima in PTHA18 and the nonlinear model.
# This can help us interpret differences in the stage vs exceedance-rate for each
# (which are due to a mixture of model differences, and monte carlo sampling).
#
plot_tsunami_maxima_in_nonlinear_model_and_PTHA18(target_point, scenario_max_stages)

#
# Get all the events data on all the source zones (avoid a heavy download if possible)
#
unique_sources = names(all_max_stage_ptha18)
events_file_local = paste0('all_events_', paste0(unique_sources, collapse="_"), '.RDS')
if(!file.exists(events_file_local)){
    all_events = lapply(unique_sources, ptha18$get_source_zone_events_data)
    names(all_events) = unique_sources
    saveRDS(all_events, events_file_local)
}else{
    all_events = readRDS(events_file_local)
}


#
# Exceedance-rates: Mean curve from PTHA18
#
ptha18_curves = vector(mode='list', length=length(all_max_stage_ptha18))
names(ptha18_curves) = names(all_max_stage_ptha18)
for(nm in names(ptha18_curves)){
    exrates = sapply(STAGE_POINTS_FOR_EXRATE_CURVE,
        function(x){ sum(all_events[[nm]]$events$rate_annual *
                         (all_max_stage_ptha18[[nm]]$max_stage > x)) })
    ptha18_curves[[nm]] = data.frame(max_stage=STAGE_POINTS_FOR_EXRATE_CURVE, exrate=exrates)
}

#
# Exceedance-rates: Mean curves from importance-sampling scenarios.
#
stopifnot(all(names(importance_sampling_scenarios_logic_tree_mean) %in% names(ptha18_curves)))
# Compared to the PTHA18 exceedance-rate curves, it helps to adjust the max-stage increments
# to account for the nonlinear model MSL (although not strictly needed)
output_max_stages_IS = STAGE_POINTS_FOR_EXRATE_CURVE + nonlinear_model_MSL
# Get the nonlinear-model max-stage for each importance_sampling_scenario
nonlinear_model_curves = list()
for(nm in names(importance_sampling_scenarios_logic_tree_mean)){

    k = which(scenario_max_stages$nonlinear_source_name == nm)
    # Index from the randomly sampled scenarios (which can include repeated scenarios)
    # to the scenario_max_stages (which only include each scenario once)
    mtch = match(importance_sampling_scenarios_logic_tree_mean[[nm]]$inds,
        scenario_max_stages$nonlinear_scenario_row[k])

    # Key terms used to calculate monte carlo exceedance rates with importance sampling.
    max_stages_sampled_scenarios = scenario_max_stages$max_stage[k[mtch]]
    max_stages_PTHA18_sampled_scenarios = scenario_max_stages$ptha18_max_stage[k[mtch]]
    rates_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean[[nm]]$event_rate_logic_tree_mean
    sampling_prob_sampled_scenarios = importance_sampling_scenarios_logic_tree_mean[[nm]]$sampling_prob

    # Double check we didn't make a mistake indexing into the data.frame's above
    r1 = scenario_max_stages$nonlinear_scenario_row[k[mtch]]
    r2 = importance_sampling_scenarios_logic_tree_mean[[nm]]$inds
    stopifnot(all(r1 == r2) & (length(r1) == length(r2)))
    
    # Get the importance-sampling-based exceedance rates, and an approximate variance
    exrates_with_uncertainty = isu$estimate_exrate_and_variance_sampled_scenarios(
            rates_sampled_scenarios=rates_sampled_scenarios,
            sampling_prob_sampled_scenarios=sampling_prob_sampled_scenarios,
            stage_sampled_scenarios = max_stages_sampled_scenarios,
            desired_stage_thresholds = output_max_stages_IS)

    # Pack the data into a data-frame with useful names
    output_df = data.frame(max_stage = output_max_stages_IS,
        exrate=exrates_with_uncertainty[,1],
        exrate_variance = exrates_with_uncertainty[,2])


    # As above, but using the PTHA18 max-stages.
    # This lets us see how good the approximation would be IF the nonlinear
    # model tsunami agreed exactly with PTHA18 in terms of max-stage. This is
    # useful because, with that assumption, we know all the error is due to
    # Monte Carlo sampling.
    exrates_with_uncertainty_using_PTHA18_max_stage = isu$estimate_exrate_and_variance_sampled_scenarios(
            rates_sampled_scenarios=rates_sampled_scenarios,
            sampling_prob_sampled_scenarios=sampling_prob_sampled_scenarios,
            # Pass the PTHA18 max-stage, but add in the nonlinear model MSL so it's comparable to the real
            # nonlinear model calculations
            stage_sampled_scenarios = (max_stages_PTHA18_sampled_scenarios + nonlinear_model_MSL),
            desired_stage_thresholds = output_max_stages_IS)

    # Store results assuming the PTHA18 agrees exactly with the nonlinear model, for testing.
    output_df = cbind(output_df, data.frame(
        exrate_using_ptha18_stages_plus_nonlinear_model_MSL=exrates_with_uncertainty_using_PTHA18_max_stage[,1],
        exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL=exrates_with_uncertainty_using_PTHA18_max_stage[,2]))

    nonlinear_model_curves[[nm]] = output_df
}

# Useful to store for later cross-checks in probabilistic raster calculations.
combined_output = list(nonlinear_model_curves=nonlinear_model_curves, ptha18_curves=ptha18_curves)
rds_file = paste0('nonlinear_model_curves_and_ptha18_curves_', target_point[1], '_', target_point[2], '.RDS')
saveRDS(combined_output, rds_file)
rm(combined_output); gc()

#
# Sum the logic-tree-mean exceedance-rate curves from PTHA18 and nonlinear model
# over the source-zones
#
ptha18_combined_curve = ptha18_curves[[1]]
ptha18_combined_curve[,2] = 0
nonlinear_model_combined_curve = nonlinear_model_curves[[1]]
nonlinear_model_combined_curve$exrate = 0
nonlinear_model_combined_curve$exrate_variance=0
nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL = 0
nonlinear_model_combined_curve$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL = 0
for(nm in names(nonlinear_model_curves)){
    # Sum of PTHA18 curves
    stopifnot(all(ptha18_curves[[nm]][,1] == ptha18_combined_curve[,1]))
    ptha18_combined_curve[,2] = ptha18_combined_curve[,2] + ptha18_curves[[nm]][,2]

    # Sum of nonlinear model curves (including the Monte-Carlo variance)
    stopifnot(all(nonlinear_model_curves[[nm]]$max_stage == nonlinear_model_combined_curve$max_stage))
    nonlinear_model_combined_curve$exrate =
        nonlinear_model_combined_curve$exrate + nonlinear_model_curves[[nm]]$exrate
    nonlinear_model_combined_curve$exrate_variance =
        nonlinear_model_combined_curve$exrate_variance + nonlinear_model_curves[[nm]]$exrate_variance

    # For testing, we also sum nonlinear_model_curves using the PTHA18 max-stages.
    # This is useful for testing, since any difference purely reflects Monte-Carlo sampling, whereas
    # for the "real" results the differences are a mixture of Monte Carlo sampling and differences
    # in the linear/nonlinear models.
    nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL = 
        nonlinear_model_combined_curve$exrate_using_ptha18_stages_plus_nonlinear_model_MSL +
          nonlinear_model_curves[[nm]]$exrate_using_ptha18_stages_plus_nonlinear_model_MSL
    nonlinear_model_combined_curve$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL = 
        nonlinear_model_combined_curve$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL +
          nonlinear_model_curves[[nm]]$exrate_variance_using_ptha18_stages_plus_nonlinear_model_MSL

}

# Make plots
plot_logic_tree_mean_exrate_curves_from_nonlinear_model_and_PTHA18(
    target_point, nonlinear_model_curves, nonlinear_model_combined_curve, nonlinear_model_MSL,
    ptha18_curves, ptha18_combined_curve, EXRATE_PLOT_YLIM)


if(DO_EPISTEMIC_UNCERTAINTY_CALCULATIONS){
    #
    # Epistemic uncertainty
    #

    # PTHA18 detailed reader
    ptha18_detailed = new.env()
    file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_detailed_PTHA18_source_zone_info.R'
    file_home = '../../../../../../../AustPTHA/CODE/ptha/ptha_access/get_detailed_PTHA18_source_zone_info.R'
    source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18_detailed, chdir=TRUE)

    # We do the calculations using the regular nonlinear model results, and separately, by
    # replacing the nonlinear model stages by the PTHA18 value. The latter is useful for testing
    # because it eliminates differences due to the different hydrodynamic approximations in PTHA18
    # vs the nonlinear model. Thus we can consider differences due only to sampling.
    for(max_stage_type in c('nonlinear_model', 'ptha18')){



        # For each source zone we have "unsegmented" and (possibly) "union-of-segments" models.
        # For the unsegmented model and each segment, the PTHA18 scenario conditional probability (given Mw)
        # is available from e.g.
        #   x = ptha18_detailed$get_PTHA18_scenario_conditional_probability_and_rates_on_segment(
        #      'kermadectonga2', 'hikurangi')
        # where the second argument is "" for the unsegmented source. 
        # The scenario probabilities (conditional on Mw) are in 
        #   x$HS_prob_given_Mw
        # and we will use the values for the sampled scenarios only. 
        # The logic tree branches and probabilities are under, e.g., 
        # tmp = ptha18_detailed$crs_data$source_envs$kermadectonga2_hikurangi$mw_rate_function(NA, 
        #    return_all_logic_tree_branches=TRUE)
        # From this we only need the Mw bin rate for each logic-tree branch, and so can calculate
        # r_i(e) from that and the conditional probability 
        #   [r_i(e) = scenario-e-mwbin-rate_model_i * scenario-e-cond-prob-given-Mw]
        # Then for each logic-tree branch we can directly compute the Monte Carlo exceedance-rate -vs- stage-threshold curves. 
        # Combined with the weights for each rate model we have a distribution of exceedance rates for any stage threshold.
        # - The distribution F^{-1} (percentile --> exceedance-rate) can be evaluated using rptha::weighted_percentile
        # - We can then sample from F^{-1} using interpolation at random percentiles (random_percentile --> random_exceedance-rate)
        # - On source zones where there is only an unsegmented model, the 84th percentile of F^{-1} is of interest
        #   + We could compute it directly, or by randomly sampling from F^{-1} and taking the 84th percentile of that. 
        #   + The latter approach generalises better in the more complex case below.
        # - Consider a source-zone that combines unsegmented and 'union of segments' source representations. 
        #   For each stage threshold, the union of segments result can be obtained by (repeatedly) sampling F^{-1} at
        #   random_percentiles on each segment and summing over segments. To get co-monotonic
        #   results, the sum involves the same random percentile on each segment (whereas
        #   for independent results, the random percentile would be different).
        # - When there is a 50:50 mixture of unsegmented and union-of-segments models, we use
        #   the random approach with a 50:50 chance of each model.
        
        # Loop over each source zone
        percentile_uncertainty_results = vector(mode='list', length=length(importance_sampling_scenarios_logic_tree_mean))
        names(percentile_uncertainty_results) = names(importance_sampling_scenarios_logic_tree_mean)
        for(source_zone in names(importance_sampling_scenarios_logic_tree_mean)){

            # Get names of the segmented sub-sources (which could be empty on sources without segmentation)
            source_segmentation_info = ptha18_detailed$get_unsegmented_and_segmented_source_names_on_source_zone(
                source_zone)
            segmented_source_names = source_segmentation_info$segments
            segment_names = gsub(paste0(source_zone, '_'), '', segmented_source_names)
            unsegmented_source_name = source_segmentation_info$unsegmented_name
            stopifnot(unsegmented_source_name == source_zone)

            # Get the max stage for each sampled scenario on the current source_zone, ordered
            # in the same was as importance_sampling_scenarios_logic_tree_mean[[source_zone]]
            k = which(scenario_max_stages$nonlinear_source_name == source_zone)
            scenario_max_stages_SZ = scenario_max_stages[k,]
            mtch = match(importance_sampling_scenarios_logic_tree_mean[[source_zone]]$inds,
                scenario_max_stages_SZ$nonlinear_scenario_row)
                # Note that mtch might contain repeated values because scenario_max_stages_SZ will only contain each scenario once,
                # whereas importance_sampling_scenarios_logic_tree_mean might contain scenarios that have been sampled more than once.
            if(max_stage_type == 'nonlinear_model'){
                # Use the high-res nonlinear model
                max_stages_sampled_scenarios = scenario_max_stages_SZ$max_stage[mtch]
            }else if(max_stage_type == 'ptha18'){
                # Use PTHA18 -- this enables us to isolate the effect of scenario sampling (i.e. no differences in
                # the hydrodynamic models).
                max_stages_sampled_scenarios = scenario_max_stages_SZ$ptha18_max_stage[mtch] + nonlinear_model_MSL
            }else{
                stop(paste0('unrecognized max_stage_type: ', max_stage_type))
            }

            # For every unsegmented rate model, compute a monte carlo stage-vs-exrate curve
            unsegmented_stage_vs_exrate_curves = estimate_stage_exrate_curves_with_IS_for_all_logic_tree_branches(
                ptha18_detailed,
                source_zone,
                source_zone_segment="",
                importance_sampling_scenarios_logic_tree_mean_on_source_zone=
                    importance_sampling_scenarios_logic_tree_mean[[source_zone]],
                sampled_scenarios_max_stage=max_stages_sampled_scenarios,
                threshold_stage_values = epistemic_uncertainty_threshold_stage_values)

            # As above, segmented models
            segmented_stage_vs_exrate_curves = vector(mode='list', length=length(segment_names))
            if(length(segment_names) > 0){
                names(segmented_stage_vs_exrate_curves) = segment_names
                for(nm_i in segment_names){
                    segmented_stage_vs_exrate_curves[[nm_i]] =
                        estimate_stage_exrate_curves_with_IS_for_all_logic_tree_branches(
                            ptha18_detailed,
                            source_zone,
                            source_zone_segment=nm_i,
                            importance_sampling_scenarios_logic_tree_mean_on_source_zone=
                                importance_sampling_scenarios_logic_tree_mean[[source_zone]],
                            sampled_scenarios_max_stage=max_stages_sampled_scenarios,
                            threshold_stage_values = epistemic_uncertainty_threshold_stage_values)
                }
            }

            # Get exceedance-rate percentiles for the combination of sources above.
            set.seed(epistemic_uncertainty_random_seed) # Reproducible randomness
            percentile_uncertainty_results[[source_zone]] =
                ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
                    unsegmented_stage_vs_exrate_curves,
                    segmented_stage_vs_exrate_curves,
                    N = epistemic_uncertainty_Nsamples,
                    unsegmented_wt=0.5 + 0.5*(length(segment_names) == 0),
                    union_of_segments_wt=0.5 - 0.5*(length(segment_names) == 0),
                    segments_copula_type = 'comonotonic',
                    print_progress_every_nth_threshold=99999)

        }

        if(max_stage_type == 'nonlinear_model'){
            # Normal case
            out_file = paste0('epistemic_percentile_uncertainties_nonlinear_model_',
                target_point[1], '_', target_point[2], '.RDS')
        }else if(max_stage_type == 'ptha18'){
            # Test case
            out_file = paste0('epistemic_percentile_uncertainties_replacing_nonlinear_model_with_ptha18_',
                target_point[1], '_', target_point[2], '.RDS')
        }
        saveRDS(percentile_uncertainty_results, out_file)

        # Get the PTHA18 results with epistemic uncertainty
        ptha18_percentile_uncertainty_results = vector(mode='list', length=length(percentile_uncertainty_results))
        names(ptha18_percentile_uncertainty_results) = names(percentile_uncertainty_results)
        for(source_zone in names(ptha18_percentile_uncertainty_results)){

            ptha18_percentile_uncertainty_results[[source_zone]] = 
                ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
                    target_point =target_point, 
                    source_name = source_zone,
                    make_plot = FALSE,
                    non_stochastic_slip_sources=FALSE,
                    only_mean_rate_curve=FALSE)
        }

        # Plot the results by source zone
        plot_epistemic_uncertainties_in_PTHA18_and_nonlinear_model_by_source_zone(
            target_point,
            percentile_uncertainty_results,
            nonlinear_model_MSL,
            ptha18_percentile_uncertainty_results,
            max_stage_type)

        #
        # Compute epistemic uncertainties in PTHA18 and nonlinear model, summed over source-zones
        #

        # Get sum of nonlinear model percentile exrates
        nonlinear_model_percentile_uncertainty = 0 * percentile_uncertainty_results[[source_zone]]$percentile_exrate
        # Get sum of PTHA18 percentile exrates
        ptha18_value_dummy = 0 * ptha18_percentile_uncertainty_results[[source_zone]]$stochastic_slip_rate 
        ptha18_vals = list(
            'stochastic_slip_rate_lower_ci' = ptha18_value_dummy,
            'stochastic_slip_rate_16pc' = ptha18_value_dummy,
            'stochastic_slip_rate_median' = ptha18_value_dummy,
            'stochastic_slip_rate_84pc' = ptha18_value_dummy,
            'stochastic_slip_rate_upper_ci' = ptha18_value_dummy)
        for(source_zone in names(percentile_uncertainty_results)){
            nonlinear_model_percentile_uncertainty = nonlinear_model_percentile_uncertainty + 
                percentile_uncertainty_results[[source_zone]]$percentile_exrate

            for(nm in names(ptha18_vals)){
                ptha18_vals[[nm]] = ptha18_vals[[nm]] + ptha18_percentile_uncertainty_results[[source_zone]][[nm]]
            }
        }

        plot_summed_exrates_with_epistemic_uncertainty_in_PTHA18_and_nonlinear_model(
            target_point,
            percentile_uncertainty_results,
            nonlinear_model_percentile_uncertainty,
            nonlinear_model_MSL,
            epistemic_uncertainty_threshold_stage_values,
            ptha18_percentile_uncertainty_results,
            ptha18_vals,
            max_stage_type)

    }
}
