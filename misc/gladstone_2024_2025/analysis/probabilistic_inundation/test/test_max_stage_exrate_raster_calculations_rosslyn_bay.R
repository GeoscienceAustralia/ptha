#
# Compare the max-stage exceedance-rate against a separately
# calculated value at a point (from '../../max_stage_at_a_point')
#
asfm = new.env()
source('../application_specific_file_metadata.R', local=asfm)

swals = new.env()
source(asfm$swals_plots_script_file, local=swals)

errc = new.env()
source('../exrate_given_threshold/exceedance_rate_raster_calculations.R', local=errc)

##
## Inputs 
##
target_point = c(150.789448242126, -23.1609999987272)

# Get "max-stage-exceedance' of "MSL + 2.5"
MSL = 0.0
exceedance_threshold = MSL + 2.5
comparison_calculation_RDS = 
    '../../max_stage_at_a_point/ptha/data/nonlinear_model_curves_and_ptha18_curves_150.789448242126_-23.1609999987272.RDS'
test_source_name = 'kermadectonga2'

# Tarred models
tarred_multidomain_dirs = Sys.glob(
    '../../../swals/OUTPUTS/ptha/sea_level_vary/random_kermadectonga2/ptha*/RUN*.tar')

# Logic tree mean scenarios
scenario_databases = list()
scenario_databases$logic_tree_mean_curve_HS = read.csv(
    "../../../sources/hazard/ptha_batch/random_kermadectonga2/random_scenarios_kermadectonga2_logic_tree_mean_HS.csv")

# Domain index
domain_containing_point = swals$find_domain_containing_point(
    matrix(target_point, nrow=1, ncol=2), 
    multidomain_dir=asfm$reference_multidomain_dir)
domain_index = swals$domain_index_from_folder(basename(domain_containing_point$domain_dir))

input_domain_index_and_scenarios_name = list(
    domain_index = domain_index, # Compare against a site on domain 39
    scenarios_name = 'logic_tree_mean_curve_HS') # Use logic-tree-mean scenrios

output_directory = 'tmp/exceedance_rate_raster_calculations_offshore'
raster_name_start = asfm$get_raster_name_stub_from_variable_name('max_stage')

##
## END INPUTS
##

for (i in 1:length(scenario_databases)) {
    scenario_databases[[i]]$md_dir = asfm$find_matching_md_data(
        scenario_databases[[i]]$inds, tarred_multidomain_dirs, test_source_name
    )
}

dir.create(output_directory, showWarnings=FALSE, recursive = TRUE)
# Remove any rasters from old tests
existing_rasters = Sys.glob(paste0(output_directory, '/*.tif'))
if(length(existing_rasters) > 0) tmp = file.remove(existing_rasters)

#
# Test 1 -- basic exceedance-rate calculation
#
make_raster = errc$compute_exceedance_rates_on_tile_with_importance_sampling(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs,
    scenario_databases,
    output_directory,
    exceedance_threshold,
    raster_name_start,
    print_progress_debug=TRUE) # If slow, see below
# Attempts to speedup NCI GADI gdata file IO access when it was running slowly.
# * Touch the files first (in parallel)
#   `mclapply(tarred_rasters, function(x) system(paste0('tar --list max_stage_domain_39.tif -f ', x)), mc.cores=48, mc.preschedule=FALSE)`
#   Then run errc$compute_exceedance_rates_on_tile(....) 
#   It definitely read more quickly, although still slow (vs when gdata is working well)

x = raster(paste0(
        output_directory, '/', 
        input_domain_index_and_scenarios_name$scenarios_name,
        '_domain_', domain_index, 
        '_max_stage_domain__exceedance_rate_with_threshold_', 
        round(exceedance_threshold, 3),'.tif'))

result = extract(x, matrix(target_point, ncol=2))

get_reference_exrate_and_variance<-function(comparison_calculation_RDS, exceedance_threshold, source_zone){
    x = readRDS(comparison_calculation_RDS)
    out_exrate = approx(
        x$nonlinear_model_curves[[source_zone]]$max_stage, 
        x$nonlinear_model_curves[[source_zone]]$exrate, 
        xout=exceedance_threshold)
    out_exrate_variance = approx(
        x$nonlinear_model_curves[[source_zone]]$max_stage, 
        x$nonlinear_model_curves[[source_zone]]$exrate_variance, 
        xout=exceedance_threshold)

    return(list(exrate=out_exrate$y, exrate_variance=out_exrate_variance$y))
}

expected_result = get_reference_exrate_and_variance(comparison_calculation_RDS, exceedance_threshold, test_source_name)
# Previously I computed the exceedance-rates at this site using very different code, see here:
#    ../max_stage_at_point/
# It is  0.0001561856
if(abs(result - expected_result$exrate) < 1e-10){
    print('PASS')
}else{
    print('FAIL')
}


#
# Test 2 -- exceedance-rate with variance. Here the exceedance-rate should
# be the same as the above up to floating-point differnces due to
# reordering of a sum, and we also get the error variance
output_directory = 'tmp/exceedance_rate_raster_calculations_offshore_with_variance'
dir.create(output_directory, showWarnings=FALSE, recursive = TRUE)
# Remove any rasters from old tests
existing_rasters = Sys.glob(paste0(output_directory, '/*.tif'))
if(length(existing_rasters) > 0) tmp = file.remove(existing_rasters)

make_raster = errc$compute_exceedance_rates_and_error_variance_on_tile_with_importance_sampling(
    input_domain_index_and_scenarios_name,
    tarred_multidomain_dirs,
    scenario_databases,
    output_directory,
    exceedance_threshold,
    raster_name_start,
    print_progress_debug=TRUE)

xnew = raster(Sys.glob('tmp/exceedance_rate_raster_calculations_offshore_with_variance/logic_tree_mean_curve_HS_domain_*_max_stage_domain__exceedance_rate_with_threshold_*.tif'))
xnew_var = raster(Sys.glob('tmp/exceedance_rate_raster_calculations_offshore_with_variance/logic_tree_mean_curve_HS_domain_*_max_stage_domain__variance_of__exceedance_rate_with_threshold_*.tif'))
result = extract(xnew, matrix(target_point, ncol=2))


# Previously I computed the exceedance-rates at this site using very different code, see here:
#    ../max_stage_at_point/
expected_result = get_reference_exrate_and_variance(comparison_calculation_RDS, exceedance_threshold, test_source_name)
if(abs(result - expected_result$exrate) < 1e-10){
    print('PASS')
}else{
    print('FAIL')
}

# Check the exceedance-rates are just like before
tol = 1.0e-08
if(all(abs(as.matrix(xnew - x)) <= tol*as.matrix(xnew), na.rm=TRUE)){
    print('PASS')
}else{
    print('FAIL')
}

# From a separate calculation, the variance should be
result = extract(xnew_var, matrix(target_point, ncol=2))
rtol = 1.0e-6
if(abs(result - expected_result$exrate_variance)/expected_result$exrate_variance < rtol){
    print('PASS')
}else{
    print('FAIL')
}
