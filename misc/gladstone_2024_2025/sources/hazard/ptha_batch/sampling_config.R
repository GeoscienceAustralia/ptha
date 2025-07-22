
# Point used to define importance criteria. This must be very close to a PTHA18 hazard point.
target_pt = c(152.667, -23.333)  # gauge ID1666.6

# Number of samples (sampling with replacement, so likely fewer scenarios to simulate)
NTOT = 360

# For reproducible randomness
REPRODUCIBLE_SEED = 123456

# When we tablulate exceedance-rate curves, use these stage thresholds
stage_thresholds = seq(0.02, 5.5, by = 0.05)

# Source zones to include, and their number of samples
source_zones_and_samples = list(
    'kermadectonga2'         = round(NTOT * 4./17.),
    'solomon2'               = round(NTOT * 3./17.),
    'newhebrides2'           = round(NTOT * 3./17.),
    'southamerica'           = round(NTOT * 3./17.),
    'outerrisenewhebrides'   = round(NTOT * 1./17.),
    'alaskaaleutians'        = round(NTOT * 1./17.),
    'kurilsjapan'            = round(NTOT * 1./17.),
    'outerrisesolomon'       = round(NTOT * 1./17.)
)

# Test point data (for plotting only). These should be very close to PTHA18 hazard points.
other_test_points = list(
    "near_offshore_target" = c(152.667, -23.667), # 350 m deep
    "rosslyn_bay_nearshore" = c(150.818, -23.166),  # 3km east of Rosslyn Bay, 12 m deep
    "gladstone_nearshore" = c(151.523, -23.741), # Near Gladstone, 25 m deep. Way too nearshore for PTHA18
    "Baffle_creek_offshore" = c(152.109, -24.403)  # 20 m deep, 13km NE from Baffle Creek
)

# For parallel initial condition tif creation
MC_CORES = 48

# Get the scripts to access the PTHA18
# Path to the get_PTHA18_results.R script
get_ptha_results_script_file = '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
source(get_ptha_results_script_file, local=ptha18, chdir=TRUE)

# Path to the get_detailed_PTHA18_source_zone_info.R script
get_detailed_ptha18_source_zone_info_script_file = paste0(
    dirname(get_ptha_results_script_file), 
    '/get_detailed_PTHA18_source_zone_info.R')

#' Read the "all source zones" PTHA18 stage vs rate curve at a point "tp", using a cache.
#'
#' The coordinate of tp must be near an actual PTHA18 point.
#'
#' @param tp the coordinate 
#' @return object with various PTHA18 stage-vs-exrate curves at TP (considering
#' all source-zones).
get_stage_vs_rate_all_source_zones<-function(tp){
    stage_vs_rate_all_source_zones_cache_file = paste0(
        'stage_vs_rate_all_source_zones_cache_', round(tp[1], 3), '_', 
        round(tp[2], 3), '.RDS')
    if(file.exists(stage_vs_rate_all_source_zones_cache_file)){
        # Read the cached version
        stage_vs_rate_all_source_zones_at_tp = readRDS(stage_vs_rate_all_source_zones_cache_file)
    }else{
        # Make and cache for later
        stage_vs_rate_all_source_zones_at_tp = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
            target_point=tp)
        saveRDS(stage_vs_rate_all_source_zones_at_tp, 
            file=stage_vs_rate_all_source_zones_cache_file)
    }
    return(stage_vs_rate_all_source_zones_at_tp)
}

# Use the "all source zones" stage vs rate curve at the target point to help
# define the importance criteria
stage_vs_rate_all_source_zones_at_target_pt = get_stage_vs_rate_all_source_zones(target_pt)
# Get curves at test points for plotting
stage_vs_rate_all_source_zones_other_test_pts = lapply(other_test_points, get_stage_vs_rate_all_source_zones)

#' Function I(e) used for importance sampling.
#'
#' The chance of sampling a scenario will be
#'    I(e) * r(e)
#' where r(e) is the logic-tree-mean rate for scenario e.
#' We use basic importance sampling to do this without bias.
#'
#' @param peak_stage_target_pt vector with peak-stage values at target_pt for
#' every scenario on one source-zone
#' @param scenario_rate_logictreemean vector with the logic-tree-mean scenario
#' rates for every scenario on one source-zone
#' @param stage_vs_rate_all_source_zones_at_target_pt Object defined above,
#' holding various PTHA18 stage-vs-rate curves (all source-zones combined).
#' This is used in an ad-hoc way, to help determine how we sample different stages.
#' The idea is to well-resolve stages around return periods of interest, using the
#' 84th percentile curve for better emphasis of large scenarios.
#' @return The importance factor I(e) for every scenario e on the source-zone. This 
#' can later be used to sample the scenarios.
importance_function<-function(peak_stage_target_pt, scenario_rate_logictreemean,
    stage_vs_rate_all_source_zones_at_target_pt){
    # Importance criteria. 

    # We want approximately even sampling of different stage bins, with bin
    # boundaries defined by the offshore PTHA at the target_pt. This is better
    # than evenly sampling stages, because in practice we will want to resolve scenarios
    # near some target exceedance-rates.
    max_possible_stage = max(peak_stage_target_pt * (scenario_rate_logictreemean > 0))
    min_possible_stage = 0
    # Find stage values corresponding to "all source zones" exceedance rates AT
    # THE 84th PERCENTILE. These will be used for bin boundaries. The choice is
    # ad-hoc (others could be made without logical problems). But the idea is to
    # evenly sample between these return periods, and using the 84th percentile curve
    # because we'll probably use the 1/2500 @ 84%. 
    exrate_boundaries = c(1/250, 1/2500, 1/10000)
    middle_stages = approx(
        stage_vs_rate_all_source_zones_at_target_pt$stochastic_slip_rate_84pc, 
        stage_vs_rate_all_source_zones_at_target_pt$stage, 
        xout=exrate_boundaries, 
        ties='min')$y

    # Drop middle stages outside min/max for this source-zone
    middle_stages = middle_stages[middle_stages < max_possible_stage] 
    middle_stages = middle_stages[middle_stages > min_possible_stage]

    # Split the exrate boundaries into multiple bins to enhance the evenness of sampling
    upper_bin_buffer = 0.05 # Make highest bin exceed largest stage
    nsub_bins = 5 # This encourage even more even sampling of stages, within the bins
    stage_threshes = approx( 
        c(min_possible_stage, middle_stages, max_possible_stage + upper_bin_buffer), 
        n=(length(exrate_boundaries)+1)*nsub_bins + 1)$y
    event_peak_stage_invrate = peak_stage_target_pt * 0

    # Importance criteria will lead to equal chance of sampling in each stage bin
    for(i in 1:(length(stage_threshes)-1)){
        ii = which(peak_stage_target_pt > stage_threshes[i] & peak_stage_target_pt <= stage_threshes[i+1])
        event_peak_stage_invrate[ii] = sum(scenario_rate_logictreemean)/sum(scenario_rate_logictreemean[ii])
        if(sum(scenario_rate_logictreemean[ii]) == 0) event_peak_stage_invrate[ii] = 0
    }
    return(event_peak_stage_invrate)
}


