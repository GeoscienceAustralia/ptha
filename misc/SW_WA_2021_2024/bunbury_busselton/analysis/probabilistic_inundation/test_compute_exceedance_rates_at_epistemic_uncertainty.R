# Find max-stage exceedance rates on domain 39 for 1.6m (which is 1m above MSL of 0.6), 84th percentile.
# This was separately computed from the gauge data in the Greater Perth model 
# `../../../[name of greater perth modelling_directory]/analysis/max_stage_at_a_point`
#
# A separate calculation yielded 
#   1.639500e-04 
# using the REPRODUCIBLE_SEED=123
#
run_command = 'Rscript compute_exceedance_rates_at_epistemic_uncertainty_percentile.R max_stage ../../swals/OUTPUTS/ptha18-BunburyBusseltonRevised-sealevel60cm/random_sunda2/ 39 0.84 1.6 tmpdir_sunda2'

system(run_command)

library(raster)
x = raster('tmpdir_sunda2/random_sunda2_max_stage_rast_threshold_2_percentile_84_subsam_1_Nrand_10000_seed_123_domain_index_39.tif')
target_point = c(113.0708, -28.5679)

result = extract(x, matrix(target_point, ncol=2))

if(abs(result - 1.639500e-04) < 1.0e-06){
    print('PASS')
}else{
    print('FAIL')
}

