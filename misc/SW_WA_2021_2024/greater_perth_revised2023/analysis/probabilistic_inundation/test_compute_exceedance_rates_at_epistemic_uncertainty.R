# Test of the epistemic uncertainty percentile calculations
#
# We find the max-stage exceedance rate at site c(113.0708, -28.5679) for a
# stage of 1.6m (noting MSL is already 0.6 m) at the 84th percentile logic-tree
# uncertainty.
#
# This was separately computed in ../max_stage_at_a_point
# The separate calculation yielded 
#   1.639500e-04 
# using the REPRODUCIBLE_SEED=123
#
# Below we compare these two calculations.
#
asfm = new.env()
source('application_specific_file_metadata.R', local=asfm)
swals = new.env()
source(asfm$swals_plots_script_file, local=swals)

# Test inputs 
# The site / stage / percentiles should correspond to values in the comparison_calculation_RDS
target_point = c(113.0708, -28.5679)
stage_threshold = 1.6
epistemic_uncertainty_percentile = 0.84
sourcezone_runs_base_dir = '../../swals/OUTPUTS/ptha18-GreaterPerth2023-sealevel60cm/random_sunda2/'
comparison_calculation_RDS = '../max_stage_at_a_point/epistemic_percentile_uncertainties_nonlinear_model_113.0708_-28.5679.RDS'
tmp_output_dir = 'tmpdir_sunda2'

# Remove any tifs from the tmp_output_dir (since the test will make 1, which should be the only tif in the directory)
existing_tifs = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
if(length(existing_tifs) > 0){
    clear_file = file.remove(existing_tifs)
}

domain_containing_point = swals$find_domain_containing_point(matrix(target_point, nrow=1, ncol=2), 
    multidomain_dir=asfm$reference_multidomain_dir)
domain_index = swals$domain_index_from_folder(basename(domain_containing_point$domain_dir))

# Run the epistemic uncertainty calculations on the raster containing the point of interest
run_command = paste0('Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R max_stage ', 
    sourcezone_runs_base_dir, ' ', 
    domain_index, ' ', 
    epistemic_uncertainty_percentile, ' ',
    stage_threshold, ' ',
    tmp_output_dir)
system(run_command)

# Extract the value at the point of interest
library(raster)
test_raster = Sys.glob(paste0(tmp_output_dir, '/*.tif'))
stopifnot(length(test_raster) == 1) # There should only be 1 raster here
x = raster(test_raster[1])
result = extract(x, matrix(target_point, ncol=2))

get_comparison_value<-function(comparison_calculation_RDS, stage_threshold, 
    epistemic_uncertainty_percentile){
    # Extract the "same" value from separate calculations in ../max_stage_at_a_point
    # They include exactly the value of interest
    x = readRDS(comparison_calculation_RDS)
    ind_j = which.min(abs(x$sunda2$threshold_stages - stage_threshold))
    ind_i = which.min(abs(x$sunda2$percentile_probs - epistemic_uncertainty_percentile))
    return(x$sunda2$percentile_exrate[ind_i, ind_j])
}
comparison_value = get_comparison_value(comparison_calculation_RDS,
    stage_threshold, epistemic_uncertainty_percentile)

# Check they agree
tol = abs(result)*1e-04
if(abs(result - comparison_value) < tol){
    print('PASS')
}else{
    print('FAIL')
}

