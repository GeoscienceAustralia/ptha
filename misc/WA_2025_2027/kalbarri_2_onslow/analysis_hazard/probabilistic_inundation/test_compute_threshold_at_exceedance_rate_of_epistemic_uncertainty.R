# Test of the epistemic uncertainty percentile calculations
#
# We find the max-stage threshold at a site for a for a given exceedance-rate
# at the 84th percentile logic-tree uncertainty, considering the sum over
# source zones.
#
# This was separately computed in ../max_stage_at_a_point
#
# Below we compare these two calculations.
#
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)
swals = new.env()
source(asfm$swals_plots_script_file, local=swals)

# Test inputs 
# The site / stage / percentiles should correspond to values in the comparison_calculation_RDS
target_point = c(113.65335, -21.50723)
stage_threshold = 2.0 # We will back calculate this
uniroot_tol = 1.0e-03 # Error tolerance for uniroot.
extra_stage_tol = 0.0 # Additional error tolerance (probably not needed anymore since we use get_allowable_max_stage_range)
epistemic_uncertainty_percentile = 0.84
comparison_calculation_RDS = '../max_stage_at_a_point/epistemic_percentile_uncertainties_nonlinear_model_113.65335_-21.50723.RDS'
max_stage_RDS = '../max_stage_at_a_point/run_stage_at_target_point_113.65335_-21.50723.RDS'
tmp_output_dir = 'tmpdir_threshold_at_epistemic_uncertainty'

# Remove any tifs from the tmp_output_dir (since the test will make 1, which should be the only tif in the directory)
existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
if(length(existing_tifs) > 0){
    clear_file = file.remove(existing_tifs)
}

# We have a separate calculation of the solution, from which we can get the exceedance_rate
# to pass to the raster calculation (so we can check if it produces the correct stage_threshold).
get_exrate_for_test<-function(comparison_calculation_RDS, stage_threshold, 
    epistemic_uncertainty_percentile){
    # Extract the "same" value from separate calculations in ../max_stage_at_a_point
    # They include exactly the value of interest
    x = readRDS(comparison_calculation_RDS)

    # Ensure 'stage_threshold' corresponds to one of the tabulated stages
    stopifnot(min(abs(x$sunda2$threshold_stages - stage_threshold)) < 1.0e-06)

    # Find the exceedance-rate at stage threshold (sum over source zones)
    ind_j = which.min(abs(x$sunda2$threshold_stages - stage_threshold))
    ind_i = which.min(abs(x$sunda2$percentile_probs - epistemic_uncertainty_percentile))
    return(x$sunda2$percentile_exrate[ind_i, ind_j] + x$outerrisesunda$percentile_exrate[ind_i, ind_j])
}
exrate_for_test = get_exrate_for_test(comparison_calculation_RDS,
    stage_threshold, epistemic_uncertainty_percentile)

get_allowable_max_stage_range<-function(max_stage_RDS, stage_threshold){
    # Get the values in the scenarios file max_stage either side of the stage_threshold.
    # The stage threshold matching the target exceedance-rate could
    # reasonably be anywhere between these values.

    x = readRDS(max_stage_RDS)
    max_stages = x$max_stage

    # Find the closest value below and above the stage_threshold
    sorted_stages = sort(max_stages)
    max_below = max(sorted_stages[sorted_stages < stage_threshold])
    max_above = min(sorted_stages[sorted_stages > stage_threshold])

    return(c(max_below, max_above))
}
range_allowable_max_stage = get_allowable_max_stage_range(max_stage_RDS, stage_threshold)

domain_containing_point = swals$find_domain_containing_point(matrix(target_point, nrow=1, ncol=2), 
    multidomain_dir=asfm$reference_multidomain_dir)
domain_index = swals$domain_index_from_folder(basename(domain_containing_point$domain_dir))

# Run the epistemic uncertainty calculations on the raster containing the point of interest
run_command = paste0('Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R', ' ',
    'max_stage', ' ', 
    domain_index, ' ', 
    epistemic_uncertainty_percentile, ' ',
    exrate_for_test, ' ',
    'adaptive_minimum', ' ',
    'adaptive_maximum', ' ',
    uniroot_tol, ' ',
    tmp_output_dir)
system(run_command)

# Extract the value at the point of interest
library(raster)
test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
stopifnot(length(test_raster) == 1) # There should only be 1 raster here
x = raster(test_raster[1])
result = extract(x, matrix(target_point, ncol=2))

# We should have got close to stage_threshold
too_low = result < range_allowable_max_stage[1] - (uniroot_tol + extra_stage_tol)
too_high = result > range_allowable_max_stage[2] + (uniroot_tol + extra_stage_tol)
if(too_low | too_high){
    print('FAIL')
}else{
    print('PASS')
}

