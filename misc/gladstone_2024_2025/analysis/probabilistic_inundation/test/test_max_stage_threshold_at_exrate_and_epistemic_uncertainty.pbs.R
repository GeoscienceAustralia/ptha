# Test of the epistemic uncertainty percentile calculations
#
# We find the max-stage threshold at a site for a for a given exceedance-rate
# at the 84th percentile logic-tree uncertainty, considering the sum over
# source zones.
#
# This was separately computed in ../../max_stage_at_a_point
#
# Below we compare these two calculations.
#
library(raster)


asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)
swals = new.env()
source(asfm$swals_plots_script_file, local=swals)

target_point = rbind(c(152.66667175293, -23.333333961162))
# The site / stage / percentiles should correspond to values in the comparison_calculation_RDS
stage_threshold = 0.0 + 0.5  # MUST correspond to a stage threshold in the RDS file
uniroot_tol = 0.001 # Error tolerance for uniroot.
epistemic_uncertainty_percentile = 0.84
comparison_calculation_RDS_epi =
    'max_stage_at_a_point/ptha/data/epistemic_percentile_uncertainties_nonlinear_model_152.66667175293_-23.333333961162.RDS'
comparison_calculation_RDS_scenarios =
    'max_stage_at_a_point/ptha/data/scenario_max_stages_152.66667175293_-23.333333961162.RDS'

tmp_output_dir = 'tmp/threshold_at_epistemic_uncertainty'
dir.create(tmp_output_dir, showWarnings=FALSE)
path_test_from_thresh = '../test'
path_thresh_from_test = '../threshold_given_exrate'

# When comparing the solution with the one computed elsewhere, there could be differences due to uniroot_tol, or these aspects of the main script
# - interpolation via neighbour cells (if SUB_SAMPLE > 1), 
# - the discreteness of the stage-vs-exrate-at-epistemic-uncertainty curve (which makes the inverse ambiguous). 
# We give an extra tolerance to account for that. This might need modification for other sites.
extra_stage_tol = 2.0e-03


# Remove any tifs from the tmp_output_dir (since the test will make 1, which should be the only tif in the directory)
existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
if(length(existing_tifs) > 0){
    clear_file = file.remove(existing_tifs)
}

# We have a separate calculation of the solution, from which we can get the exceedance_rate
# to pass to the raster calculation (so we can check if it produces the correct stage_threshold).
get_exrate_for_test<-function(comparison_calculation_RDS, stage_threshold, 
    epistemic_uncertainty_percentile){
    # Extract the "same" value from separate calculations in ../../max_stage_at_a_point
    # They include exactly the value of interest
    x = readRDS(comparison_calculation_RDS)

    nms = names(x)

    result = 0
    for(nm in nms){
        # Ensure 'stage_threshold' corresponds to one of the tabulated stages
        stopifnot(min(abs(x[[nm]]$threshold_stages - stage_threshold)) < 1.0e-06)

        # Find the exceedance-rate at stage threshold (sum over source zones)
        ind_j = which.min(abs(x[[nm]]$threshold_stages - stage_threshold))
        ind_i = which.min(abs(x[[nm]]$percentile_probs - epistemic_uncertainty_percentile))

        result = result + x[[nm]]$percentile_exrate[ind_i, ind_j]
    }

    return(result)
}

get_allowable_max_stage<-function(comparison_calculation_RDS_scenarios, stage_threshold){
    # Extract the "same" value from separate calculations in ../../max_stage_at_a_point
    # Get the values in the scenarios file $max_stage either side of the stage_threshold

    x = readRDS(comparison_calculation_RDS_scenarios)
    max_stages = x$max_stage
    # ind_below = which.max(max_stages[max_stages <= stage_threshold])
    # ind_above = which.min(max_stages[max_stages >= stage_threshold])

    # return(c(max_stages[ind_below], max_stages[ind_above]))

    # Find the closest value below and above the stage_threshold
    sorted_stages = sort(max_stages)
    max_below = max(sorted_stages[sorted_stages <= stage_threshold])
    max_above = min(sorted_stages[sorted_stages >= stage_threshold])

    return(c(max_below, max_above))
}

exrate_for_test = get_exrate_for_test(comparison_calculation_RDS_epi,
    stage_threshold, epistemic_uncertainty_percentile)


domain_containing_point = swals$find_domain_containing_point(
    target_point, 
    multidomain_dir=asfm$reference_multidomain_dir
)
domain_index = swals$domain_index_from_folder(basename(domain_containing_point$domain_dir))

range_allowable_max_stage = get_allowable_max_stage(comparison_calculation_RDS_scenarios, stage_threshold)

# # TEST 1 
# # only compute at target point. Use a fixed search range
# print('TEST 1: Only compute at target point. Use a fixed search range')
# # The target point is passed as an environment variable and is used in
# # epistemic_uncertainty_percentile.R -> make_all_pixel_data() 
# # to compute the exceedance rate only at the target point.
# TARGET_POINT_GLOBAL = target_point
# Sys.setenv(TARGET_POINT_GLOBAL = paste(TARGET_POINT_GLOBAL, collapse = ","))
# run_command = paste0('cd ', path_thresh_from_test, '; Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R', ' ',
#     'max_stage', ' ', 
#     domain_index, ' ', 
#     epistemic_uncertainty_percentile, ' ',
#     exrate_for_test, ' ',
#     stage_threshold - 1.0, ' ',
#     stage_threshold + 1.0, ' ',
#     uniroot_tol, ' ',
#     paste(path_test_from_thresh, tmp_output_dir, sep="/"))
# print(run_command)
# system(run_command)

# # Extract the value at the point of interest
# test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
# stopifnot(length(test_raster) == 1) # There should only be 1 raster here
# x = raster(test_raster[1])
# result = extract(x, target_point)

# print(paste0('Result: ', result))
# print(paste0('Stage threshold: ', stage_threshold))
# print(paste0('Min allowable max-stage: ', range_allowable_max_stage[1]))
# print(paste0('Max allowable max-stage: ', range_allowable_max_stage[2]))
# print(paste0('Tolerance: ', uniroot_tol + extra_stage_tol))

# too_low = result < range_allowable_max_stage[1] - (uniroot_tol + extra_stage_tol)
# too_high = result > range_allowable_max_stage[2] + (uniroot_tol + extra_stage_tol)
# if(too_low || too_high){
#     print('FAIL')
# }else{
#     print('PASS')
# }

# # TEST 2 
# # Only compute at target point. Use an adaptive search range
# print('TEST 2: Only compute at target point. Use an adaptive search range')
# tmp_output_dir = 'tmpdir_threshold_at_epistemic_uncertainty_adaptive'
# dir.create(tmp_output_dir, showWarnings=FALSE)
# existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
# if(length(existing_tifs) > 0){
#     clear_file = file.remove(existing_tifs)
# }
# # Run the epistemic uncertainty calculations on the raster containing the point of interest
# run_command = paste0('cd ', path_thresh_from_test, '; Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R', ' ',
#     'max_stage', ' ', 
#     domain_index, ' ', 
#     epistemic_uncertainty_percentile, ' ',
#     exrate_for_test, ' ',
#     'adaptive_minimum ',
#     'adaptive_maximum ',
#     uniroot_tol, ' ',
#     paste(path_test_from_thresh, tmp_output_dir, sep="/"))
# print(run_command)
# system(run_command)

# # Extract the value at the point of interest
# test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
# stopifnot(length(test_raster) == 1) # There should only be 1 raster here
# x = raster(test_raster[1])
# result = extract(x, target_point)

# print(paste0('Result: ', result))
# print(paste0('Stage threshold: ', stage_threshold))
# print(paste0('Allowable max stage: ', range_allowable_max_stage))
# print(paste0('Difference: ', result-stage_threshold))
# print(paste0('Tolerance: ', uniroot_tol + extra_stage_tol))

# too_low = result < range_allowable_max_stage[1] - (uniroot_tol + extra_stage_tol)
# too_high = result > range_allowable_max_stage[2] + (uniroot_tol + extra_stage_tol)
# if(too_low || too_high){
#     print('FAIL')
# }else{
#     print('PASS')
# }

# TEST 3
# Compute on whole domain. Use a adaptive search range
print('TEST 3: Compute on whole domain. Use an adaptive search range')
Sys.unsetenv("TARGET_POINT_GLOBAL")
tmp_output_dir = 'tmp/threshold_at_epistemic_uncertainty_adaptive_domain'
dir.create(tmp_output_dir, showWarnings=FALSE)
existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
if(length(existing_tifs) > 0){
    clear_file = file.remove(existing_tifs)
}
run_command = paste0('cd ', path_thresh_from_test, '; Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R', ' ',
    'max_stage', ' ', 
    domain_index, ' ', 
    epistemic_uncertainty_percentile, ' ',
    exrate_for_test, ' ',
    'adaptive_minimum ',
    'adaptive_maximum ',
    uniroot_tol, ' ',
    paste(path_test_from_thresh, tmp_output_dir, sep="/"))
print(run_command)
system(run_command)

# Extract the value at the point of interest
test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
stopifnot(length(test_raster) == 1) # There should only be 1 raster here
x = raster(test_raster[1])
result = extract(x, target_point)

print(paste0('Result: ', result))
print(paste0('Stage threshold: ', stage_threshold))
print(paste0('Allowable max stage: ', range_allowable_max_stage))
print(paste0('Difference: ', result-stage_threshold))
print(paste0('Tolerance: ', uniroot_tol + extra_stage_tol))

too_low = result < range_allowable_max_stage[1] - (uniroot_tol + extra_stage_tol)
too_high = result > range_allowable_max_stage[2] + (uniroot_tol + extra_stage_tol)
if(too_low || too_high){
    print('FAIL')
}else{
    print('PASS')
}

# clear all temp directories
tmp_dirs <- Sys.glob('tmp*')
for(tmp_dir in tmp_dirs){
    clear_dir = unlink(tmp_dir, recursive=TRUE)
}
