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
STAGE_TOL = 0.001 # Error tolerance for uniroot.
epistemic_uncertainty_percentile = 0.84
comparison_calculation_RDS = '../max_stage_at_a_point/epistemic_percentile_uncertainties_nonlinear_model_154.6666_-28.33333.RDS'
tmp_output_dir = 'tmpdir_threshold_at_epistemic_uncertainty'

# When comparing the solution with the one computed elsewhere, there could be differences due to STAGE_TOL, or these aspects of the main script
# - interpolation via neighbour cells (if SUB_SAMPLE > 1), 
# - the discreteness of the stage-vs-exrate-at-epistemic-uncertainty curve (which makes the inverse ambiguous). 
# We give an extra tolerance to account for that. This might need modification for other sites.
EXTRA_STAGE_TOL = 2.0e-03

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
exrate_for_test = get_exrate_for_test(comparison_calculation_RDS,
    stage_threshold, epistemic_uncertainty_percentile)


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
    STAGE_TOL, ' ',
    tmp_output_dir)
system(run_command)

# Extract the value at the point of interest
library(raster)
test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
stopifnot(length(test_raster) == 1) # There should only be 1 raster here
x = raster(test_raster[1])
result = extract(x, matrix(target_point, ncol=2))

if(abs(result - stage_threshold) < (STAGE_TOL + EXTRA_STAGE_TOL)){
    print('PASS')
}else{
    print(c('FAIL', result, stage_threshold, result-stage_threshold))
}

