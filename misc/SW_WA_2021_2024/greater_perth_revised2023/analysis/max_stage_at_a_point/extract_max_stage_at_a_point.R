#
# Plots comparing max-stage and exceedance-rates in PTHA18 (frictionless, linear solver, 36hrs)
# with nonlinear model (high-res, 24 hours) using stratified/importance sampling.
#

# SWALS routines (used to find the domain containing a point)
swals = new.env()
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R', local=swals)

# STARS to pick cells from rasters
library(stars)

# PTHA18 reader
ptha18 = new.env()
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '../../../../../../AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

#
# INPUTS
#

# Use a command-line integer argument to choose which point to plot
point_code = as.numeric(commandArgs(trailingOnly=TRUE)[1])

# Define the background sea-level in the nonlinear model (so we can compare nonlinear model
# wave-heights to PTHA18)
nonlinear_model_MSL = 0.6

# Working with tarred raster files
tarred_raster_files =
    Sys.glob('../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/*/*/raster_output_files.tar')
# Each raster_output_files.tar contains a set of rasters, including 'raster_filename'

# A multidomain directory with design that matches the one used for the
# importance-sampling model runs (in terms of the domain layout).
# Will be used to associate the target_point with a particular
# max_stage_domain_XX.tif raster.
reference_multidomain_dir =
    '../../swals/OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20230831_172007573/'

# How many cores (cannot trust auto-detect on NCI).
# The script isn't very parallel, but assuming we need 16GB memory, then on the Cascade Lake nodes we may as well use 
# the corresponding number of cores { = 16/(192/48) cores }
MC_CORES = 4

# Files with scenarios used for importance-sampling -- logic_tree_mean curve results
importance_sampling_scenarios_logic_tree_mean = list()
importance_sampling_scenarios_logic_tree_mean$outerrisesunda =
    read.csv('../../sources/hazard/random_outerrisesunda/random_scenarios_outerrisesunda_logic_tree_mean_curve_HS.csv')
importance_sampling_scenarios_logic_tree_mean$sunda2 =
    read.csv('../../sources/hazard/random_sunda2/random_scenarios_sunda2_logic_tree_mean_curve_HS.csv')

# Points at which we compute the exceedance-rate curve for the logic-tree-mean calculations
STAGE_POINTS_FOR_EXRATE_CURVE = seq(0.01, 5.01, len=1001)

# Y-range for plots of exceedance-rate curves for the logic-tree-mean calculations
EXRATE_PLOT_YLIM = c(1.0e-05, 0.1)

# Turn on/off the epistemic uncertainty calculations.
DO_EPISTEMIC_UNCERTAINTY_CALCULATIONS = TRUE #FALSE

# Files with scenarios used for epistemic uncertainties -- split up by source-zone and unsegmented/segmented
importance_sampling_scenarios_detailed = list(
    # Each entry is a list (one per source_zone) with
    #    - random_scenarios:
    #        List of data.frames with the random_scenarios, with weights based
    #        on the unsegmented source and each segment. The names of this list
    #        must match the source names in PTHA18.
    #    - segment_names:
    #        Character vector with segment name for each entry of random_scenarios,
    #        ordered in the same way. For the unsegmented source, the entry should be "".
    #
    outerrisesunda = list(
        random_scenarios = list(
            outerrisesunda_unsegmented = read.csv(
                '../../sources/hazard/random_outerrisesunda/random_scenarios_outerrisesunda_unsegmented_HS.csv')
            ),
        # segment_names must be ordered as above, with empty string for unsegmented source
        segment_names = ""
        ),
    sunda2 = list(
        random_scenarios = list(
            # In these data.frames, the scenarios are the same in all cases. But
            # the scenario_weights/rates change, depending on the source_zone segmentation
            sunda2_unsegmented = read.csv(
                '../../sources/hazard/random_sunda2/random_scenarios_sunda2_unsegmented_HS.csv'),
            sunda2_arakan_segment = read.csv(
                '../../sources/hazard/random_sunda2/random_scenarios_sunda2_arakan_segment_HS.csv'),
            sunda2_andaman_segment = read.csv(
                '../../sources/hazard/random_sunda2/random_scenarios_sunda2_andaman_segment_HS.csv'),
            sunda2_sumatra_segment = read.csv(
                '../../sources/hazard/random_sunda2/random_scenarios_sunda2_sumatra_segment_HS.csv'),
            sunda2_java_segment = read.csv(
                '../../sources/hazard/random_sunda2/random_scenarios_sunda2_java_segment_HS.csv')
        ),
        # segment_names must be ordered as above, with empty string for unsegmented source
        segment_names = c('', 'arakan', 'andaman', 'sumatra', 'java')
        )
    )

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

if(point_code == 1){
    ## About 100m depth, NE of Rotnest Island
    target_point = c(115.234169, -31.6791019)
}else if(point_code == 2){
    # Quite offshore ~ 1000m depth
    target_point = c(114.7244, -31.62706)
}else if(point_code == 3){
    # Quite offshore ~ 1000m depth
    target_point = c(114.73777, -31.110929)
}else if(point_code == 4){
    # Very nearshore, the comparison should be terrible!
    target_point = c(115.64612, -31.80478)
}else if(point_code == 5){
    # About 100m depth
    target_point = c(115.27728, -31.850109)
}else if(point_code == 6){
    # Well offshore, 4km depth
    target_point = c(114.3333, -31.3333)
}else if(point_code == 7){
    # About 1km depth, further north than point 6
    target_point = c(114.3333, -30.3333)
}else if(point_code == 8){
    # About 1km depth, quite north
    target_point = c(113.0708, -28.5679)
}else if(point_code == 9){
    # Near point 8, in the leap-frog domain
    target_point = c(112.843, -28.482)
}else if(point_code == 10){
    # Fremantle boat harbour
    target_point = c(115.74660, -32.06373)
}else if(point_code == 11){
    # Coast near Penguin Island
    target_point = c(115.6994, -32.3050)
}else if(point_code == 12){
    # Hillarys Harbour
    target_point = c(115.7393, -31.8225)
}else if(point_code == 13){
    # Offshore nearer to Bunbury/Busselton
    target_point = c(114.0, -33.3333)
}else if(point_code == 14){
    # Offshore also near-ish to Bunbury/Busselton
    target_point = c(114.0, -32.6666)

}else{
    stop('unknown point code')
}


#
# END INPUTS
#

# Get the name of the max-stage raster file that contains the point
tmp = swals$find_domain_containing_point(target_point, multidomain_dir=reference_multidomain_dir)
raster_filename = paste0('max_stage_domain_', swals$domain_index_from_folder(basename(tmp$domain_dir)), '.tif'); rm(tmp)

# Get the indices of the target_point in our raster file
get_indices<-function(tarred_raster_file, raster_filename){
    # Get the indices of interest by checking one tif
    # Here I assume all tifs with the same raster_filename have the same extent (but that will be checked)
    r1 = read_stars(paste0('/vsitar/', tarred_raster_file, '/', raster_filename))
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
tmp = get_indices(tarred_raster_files[1], raster_filename)
XIND = tmp$XIND; YIND=tmp$YIND; coordinate_errtol = tmp$coordinate_errtol; rm(tmp)

# Get the max-stage at the desired point from a single file.
get_max_stage_at_point<-function(tarred_raster_file, raster_filename, XIND, YIND, target_point, coordinate_errtol){

    # Read the file
    r1 = read_stars(paste0('/vsitar/', tarred_raster_file, '/', raster_filename))
    coordinate_vals = c(
        st_get_dimension_values(r1, 'x', where='center')[XIND],
        st_get_dimension_values(r1, 'y', where='center')[YIND])

    val = r1[[1]][XIND, YIND]

    rm(r1); gc()

    # Check the coordinates are as expected
    if(any(abs(coordinate_vals - target_point) > coordinate_errtol)){
        return(NA)
    }else{
        # Return the value of interest
        return(val)
    }

}

# Get all max-stage values from the nonlinear model runs
library(parallel)
all_max_stage = mclapply(tarred_raster_files,
    get_max_stage_at_point,
    raster_filename=raster_filename,
    XIND=XIND, YIND=YIND,
    target_point=target_point,
    coordinate_errtol = coordinate_errtol,
    mc.cores=MC_CORES)

run_stage_at_target_point = data.frame(file=tarred_raster_files, max_stage=unlist(all_max_stage))
rm(all_max_stage)

output_RDS = paste0('run_stage_at_target_point_', target_point[1], '_', target_point[2], '.RDS')
saveRDS(run_stage_at_target_point, output_RDS)

#
# Now get the PTHA18 max-stage values for all these runs.
#
bdr = basename(dirname(run_stage_at_target_point$file)) # Can extract info from this piece of file path
source_name  = unlist(lapply(bdr, function(x) strsplit(x, '_')[[1]][4]))
scenario_row = unlist(lapply(bdr, function(x) as.numeric(strsplit(x, '_')[[1]][6])))
scenario_Mw  = unlist(lapply(bdr, function(x) as.numeric(strsplit(x, '_')[[1]][8])/10))

# A few sanity checks (since a change in the folder name structure could break the above code)
stopifnot(all(scenario_Mw > 7 & scenario_Mw < 10))
stopifnot(all(unique(source_name) %in% ptha18$config_env$source_names_all))
stopifnot(all(scenario_row > 0 & scenario_row == round(scenario_row)))

all_max_stage_ptha18 = ptha18$get_peak_stage_at_point_for_each_event(
    target_point=target_point,
    all_source_names = unique(source_name))

#
# Get the PTHA18 max-stage for comparison
#
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
# Compare nonlinear model and PTHA18 (36 hours, frictionless)
# PTHA18 tends to be larger due to lack of friction and longer run-time,
# but differences are smaller very offshore
#
sourceID = as.integer(as.factor(scenario_max_stages$nonlinear_source_name))
COLZ = colorRampPalette(c('green', 'darkgreen'))(max(sourceID)) # For plot

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
dev.off()

# Another way of displaying this
out_file = paste0('Tsunami_maxima_in_PTHA18_and_Nonlinear_Model_relative_', target_point[1], '_', target_point[2], '.png')
mean_max_stage = scenario_max_stages$ptha18_max_stage #0.5 * (scenario_max_stages$ptha18_max_stage + scenario_max_stages$nonlinear_max_stage_minus_MSL)
diff_max_stage = (scenario_max_stages$ptha18_max_stage - scenario_max_stages$nonlinear_max_stage_minus_MSL)
reldiff =  diff_max_stage/mean_max_stage
png(out_file, width=12, height=9, units='in', res=300)
scatter.smooth(mean_max_stage, reldiff, xlab='PTHA18 max-stage',
     ylab='[PTHA18 - Nonlinear]/PTHA18',
     col=COLZ[sourceID],
    cex.lab=1.5, cex.axis=1.5, pch=19)
abline(h=0, col='red', lty='dashed')
title(main='PTHA18 max-stage vs Relative difference in nonlinear model', cex.main=1.7)
dev.off()

# Get the events data (avoid a heavy download if possible)
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

# Exceedance-rates: Mean curves from importance-sampling scenarios.
stopifnot(all(names(importance_sampling_scenarios_logic_tree_mean) %in% names(ptha18_curves)))
# Compared to the PTHA18 exceedance-rate curves, it helps to adjust the max-stage increments
# to account for the nonlinear model MSL (although not strictly needed)
output_max_stages_IS = STAGE_POINTS_FOR_EXRATE_CURVE + nonlinear_model_MSL
# Get the nonlinear-model max-stage for each importance_sampling_scenario
nonlinear_model_curves = list()
for(nm in names(importance_sampling_scenarios_logic_tree_mean)){
    k = which(scenario_max_stages$nonlinear_source_name == nm)
    mtch = match(importance_sampling_scenarios_logic_tree_mean[[nm]]$inds,
        scenario_max_stages$nonlinear_scenario_row[k])

    # Need one max-stage for every event in the PTHA18 -- but non NA values only for modelled events
    max_stages = rep(NA, max(importance_sampling_scenarios_logic_tree_mean[[nm]]$inds))
    max_stages[importance_sampling_scenarios_logic_tree_mean[[nm]]$inds] = scenario_max_stages$max_stage[k[mtch]]

    # Get the importance-sampling-based exceedance rates, and an approximate variance
    exrates_with_uncertainty = matrix( unlist(lapply(output_max_stages_IS, function(x){
        ptha18$estimate_exrate_uncertainty(
            importance_sampling_scenarios_logic_tree_mean[[nm]],
            max_stages, threshold_stage=x)
        })), ncol=2, byrow=TRUE)

    # Pack the data into a data-frame with useful names
    output_df = data.frame(max_stage = output_max_stages_IS,
        exrate=exrates_with_uncertainty[,1],
        exrate_variance = exrates_with_uncertainty[,2])

    nonlinear_model_curves[[nm]] = output_df
}

# Useful to store for later cross-checks in probabilistic raster calculations.
combined_output = list(nonlinear_model_curves=nonlinear_model_curves, ptha18_curves=ptha18_curves)
rds_file = paste0('nonlinear_model_curves_and_ptha18_curves_', target_point[1], '_', target_point[2], '.RDS')
saveRDS(combined_output, rds_file)
rm(combined_output); gc()

#
# Compare logic-tree-mean exceedance-rate curves from PTHA18 and nonlinear model, separately for each source-zone
#
out_file = paste0('Exceedance_Rates_in_PTHA18_and_Nonlinear_Model_', target_point[1], '_', target_point[2], '.png')
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
           lower_CI, col='black', t='l', lty='dashed')
    points(nonlinear_model_curves[[nm]]$max_stage - nonlinear_model_MSL,
           upper_CI, col='black', t='l', lty='dashed')

    grid(col='orange')

    legend('topright', c('PTHA18', 'High-res nonlinear model', 'High-res 95% CI'),
           col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
           pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
}
dev.off()

#
# Compare logic-tree-mean exceedance-rate curves from PTHA18 and nonlinear model, summed over source-zones
#
ptha18_combined_curve = ptha18_curves[[1]]
ptha18_combined_curve[,2] = 0
nonlinear_model_combined_curve = nonlinear_model_curves[[1]]
nonlinear_model_combined_curve$exrate = 0; nonlinear_model_combined_curve$exrate_variance=0
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
}

out_file = paste0('Exceedance_Rates_Summed_in_PTHA18_and_Nonlinear_Model_', target_point[1], '_', target_point[2], '.png')
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
       lower_CI, col='black', t='l', lty='dashed')
points(nonlinear_model_combined_curve$max_stage - nonlinear_model_MSL,
       upper_CI, col='black', t='l', lty='dashed')

grid(col='orange')

legend('topright', c('PTHA18', 'High-res nonlinear model', 'High-res 95% CI'),
       col=c('skyblue', 'black', 'black'), cex=1.7, bty='n',
       pch=c(NA, NA, NA), lty=c('solid', 'solid', 'dashed'), lwd=c(2,2,1) )
dev.off()


if(DO_EPISTEMIC_UNCERTAINTY_CALCULATIONS){
    #
    # Epistemic uncertainty
    #

    # PTHA18 detailed reader
    ptha18_detailed = new.env()
    file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_detailed_PTHA18_source_zone_info.R'
    file_home = '../../../../../../AustPTHA/CODE/ptha/ptha_access/get_detailed_PTHA18_source_zone_info.R'
    source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18_detailed, chdir=TRUE)

    # Loop over each source zone
    percentile_uncertainty_results = vector(mode='list', length=length(importance_sampling_scenarios_detailed))
    names(percentile_uncertainty_results) = names(importance_sampling_scenarios_detailed)
    for(source_zone in names(importance_sampling_scenarios_detailed)){

        # Get the random scenarios on this source-zone, as represented for both the unsegmented model
        # and all of the segments. Note we should be using the same scenarios in all cases, but the weights
        # will change.
        random_scenarios = importance_sampling_scenarios_detailed[[source_zone]]$random_scenarios
        segment_names = importance_sampling_scenarios_detailed[[source_zone]]$segment_names

        # Make a variable like scenario_max_stages, but restricted to the current source
        k = which(scenario_max_stages$nonlinear_source_name == source_zone)
        scenario_max_stages_SZ = scenario_max_stages[k,]
        # Pack values only for the sampled events
        max_stage = rep(NA, max(scenario_max_stages_SZ$nonlinear_scenario_row))
        max_stage[scenario_max_stages_SZ$nonlinear_scenario_row] = scenario_max_stages_SZ$max_stage

        # Compute epistemic uncertainties in the exceedance-rate at these threshold stage values
        threshold_stage_values = epistemic_uncertainty_threshold_stage_values

        # Compute exceedance-rates from random scenarios for all logic-tree branches: Unsegmented model
        unsegmented_ind = which(grepl('unsegmented', names(random_scenarios)))
        stopifnot(length(unsegmented_ind) == 1)
        unsegmented_scenario_exrates_logic_tree =
            ptha18_detailed$random_scenario_exceedance_rates_all_logic_tree_branches(
                source_zone=source_zone,
                segment='',
                random_scenarios = random_scenarios[[unsegmented_ind]],
                all_scenario_stage = max_stage,
                threshold_stages = threshold_stage_values)

        # As above, segmented models
        segment_inds = which(!grepl('unsegmented', names(random_scenarios))) # Indices of segments in the lists above
        segmented_scenario_exrates_logic_tree = vector(mode='list', length=length(segment_inds))
        if(length(segment_inds) > 0){
            names(segmented_scenario_exrates_logic_tree) = names(random_scenarios)[segment_inds]
            for(i in segment_inds){
                nm_i = names(random_scenarios)[i]
                segmented_scenario_exrates_logic_tree[[nm_i]] =
                    ptha18_detailed$random_scenario_exceedance_rates_all_logic_tree_branches(
                        source_zone=source_zone,
                        segment=segment_names[i],
                        random_scenarios = random_scenarios[[nm_i]],
                        all_scenario_stage = max_stage,
                        threshold_stages = threshold_stage_values)
            }
        }

        set.seed(epistemic_uncertainty_random_seed) # Reproducible randomness
        percentile_uncertainty_results[[source_zone]] =
            ptha18$compute_exceedance_rate_percentiles_with_random_sampling(
                unsegmented_scenario_exrates_logic_tree,
                segmented_scenario_exrates_logic_tree,
                N = epistemic_uncertainty_Nsamples,
                unsegmented_wt=0.5 + 0.5*(length(segment_inds) == 0),
                union_of_segments_wt=0.5 - 0.5*(length(segment_inds) == 0),
                segments_copula_type = 'comonotonic',
                print_progress_every_nth_threshold=99999)

    }
    out_file = paste0('epistemic_percentile_uncertainties_nonlinear_model_',
        target_point[1], '_', target_point[2], '.RDS')
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

    #
    # Plot epistemic uncertainties in PTHA18 and nonlinear model, by source-zone
    #
    out_file = paste0('Exceedance_Rates_with_epistemic_uncertainty_in_PTHA18_and_Nonlinear_Model_', 
        target_point[1], '_', target_point[2], '.png')
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

        title(
            paste0('Max-stage exceedance-rate with epistemic uncertainty \n PTHA18 & nonlinear model: ', 
                   source_zone), 
            cex.main=1.6)
        legend('topright', c('2.5%', '16%', '50%', '84%', '97.5%', 'PTHA18'), 
               col=c(1:nr, 1), lty=c(rep('solid', nr), 'dotted'), pch=c(rep(1, nr), NA))
    }
    dev.off()

    #
    # Plot epistemic uncertainties in PTHA18 and nonlinear model, summed over source-zones
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


    out_file = paste0('Exceedance_Rates_Summed_with_epistemic_uncertainty_in_PTHA18_and_Nonlinear_Model_', 
        target_point[1], '_', target_point[2], '.png')
    png(out_file, width=HT*1.2, height=HT, units='in', res=300)
    options(scipen=5)

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
        points(ptha18_percentile_uncertainty_results[[source_zone]]$stage,
               ptha18_vals[[nm]], t='l', lty='dotted')
    }

    title(
        paste0('Max-stage exceedance-rate with epistemic uncertainties \n PTHA18 & nonlinear model, sum over source-zones'), 
        cex.main=1.7)
    legend('topright', c('2.5%', '16%', '50%', '84%', '97.5%', 'PTHA18'), 
           col=c(1:nr, 1), lty=c(rep('solid', nr), 'dotted'), pch=c(rep(1, nr), NA))

    dev.off()

}
