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
target_point = c(154.6666, -28.33333)
stage_threshold = 1.1 + 0.3 # MUST correspond to a stage threshold in the RDS file
uniroot_tol = 0.001 # Error tolerance for uniroot.
epistemic_uncertainty_percentile = 0.84
comparison_calculation_RDS_epistemic_uncertainties = '../max_stage_at_a_point/epistemic_percentile_uncertainties_nonlinear_model_154.6666_-28.33333.RDS'
comparison_calculation_RDS_scenarios = gsub('epistemic_percentile_uncertainties_nonlinear_model_', 'scenario_max_stages_', comparison_calculation_RDS_epistemic_uncertainties)

tmp_output_dir = 'tmpdir_threshold_at_epistemic_uncertainty'

# When comparing the solution with the one computed elsewhere, there could be differences due to uniroot_tol, or these aspects of the main script
# - interpolation via neighbour cells (if SUB_SAMPLE > 1), 
# - The discreteness of the stage-vs-exrate-at-epistemic-uncertainty curve (which makes the inverse ambiguous). 
#   - This is now dealt with by looking at the relevant discrete stage values, see get_allowable_max_stage below.
# We give an extra tolerance to account for these factors. It might need modification for other sites.
extra_stage_tol = 0.0e-03

# Remove any tifs from the tmp_output_dir (since the test will make 1, which should be the only tif in the directory)
existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
if(length(existing_tifs) > 0){
    clear_file = file.remove(existing_tifs)
}

# We have a separate calculation of the solution, from which we can get the exceedance_rate
# to pass to the raster calculation (so we can check if it produces the correct stage_threshold).
get_exrate_for_test<-function(comparison_calculation_RDS_epistemic_uncertainties, stage_threshold, 
    epistemic_uncertainty_percentile){
    # Extract the "same" value from separate calculations in ../max_stage_at_a_point
    # They include exactly the value of interest
    x = readRDS(comparison_calculation_RDS_epistemic_uncertainties)

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
exrate_for_test = get_exrate_for_test(comparison_calculation_RDS_epistemic_uncertainties,
    stage_threshold, epistemic_uncertainty_percentile)

# Determine a more realistic error tolerance
get_allowable_max_stage<-function(comparison_calculation_RDS_scenarios, stage_threshold){
    # The calculations in ../../max_stage_at_a_point involve computing the
    # exceedance-rate at a set of a-priori tabulated max-stage values, whereas
    # the raster 'threshold-at-exceedance-rate' calculations use root-finding,
    # without a-prior tabulated max-stage values. This can lead to legitimate
    # differences in the solutions. Say the scenarios actually have max-stages
    # of 
    #     s1, s2, ... s_(j-1), s_j, ... s_N
    # and exceedance-rates at the percentile of interest equal to
    #     e1, e2, ... e_(j-1), e_j, ... e_N
    # Then any "tabulated max-stage value" s_t, in between s_(j-1) and s_j,
    # will be assigned an exceedance-rate at the percentile equal to e_(j-1). 
    # But suppose I want to look up the max-stage value matching the
    # exceedance-rate e_(j-1). Due to the discreteness of the exceedance-rate
    # curve, valid solutions would be 
    #     s_(j-1) <= solution < s_(j).
    # So we need to find this range to evaluate the accuracy (and also consider the
    # uniroot tolerance, which is separate)

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
range_allowable_max_stage = get_allowable_max_stage(comparison_calculation_RDS_scenarios, stage_threshold)

domain_containing_point = swals$find_domain_containing_point(matrix(target_point, nrow=1, ncol=2), 
    multidomain_dir=asfm$reference_multidomain_dir)
domain_index = swals$domain_index_from_folder(basename(domain_containing_point$domain_dir))

# Run the epistemic uncertainty calculations on the raster containing the point of interest
run_command = paste0('Rscript compute_threshold_at_exceedance_rate_of_epistemic_uncertainty_percentile.R', ' ',
    'max_stage', ' ', 
    domain_index, ' ', 
    epistemic_uncertainty_percentile, ' ',
    exrate_for_test, ' ',
    stage_threshold - 1, ' ',
    stage_threshold + 1, ' ',
    uniroot_tol, ' ',
    tmp_output_dir)
system(run_command)

# Extract the value at the point of interest
library(raster)
test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
stopifnot(length(test_raster) == 1) # There should only be 1 raster here
x = raster(test_raster[1])
result = extract(x, matrix(target_point, ncol=2))

print(paste0('Result: ', result))
print(paste0('  error tolerance: ', uniroot_tol + extra_stage_tol))
print(paste0('Stage threshold: ', stage_threshold))
print(paste0('  spacing between nearest scenario stages: ', range_allowable_max_stage[1], ', ', range_allowable_max_stage[2]))

too_low = result < range_allowable_max_stage[1] - (uniroot_tol + extra_stage_tol)
too_high = result > range_allowable_max_stage[2] + (uniroot_tol + extra_stage_tol)
if(too_low | too_high){
    print('FAIL')
}else{
    print('PASS')
}


