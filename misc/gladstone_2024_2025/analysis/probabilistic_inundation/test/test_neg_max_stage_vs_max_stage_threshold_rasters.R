#
# Compare the max-stage thresholds at given exceedance-rate against negative max-stage's. The rasters should be identical but flipped in sign. It requires the hard-coded inputs below to be already generated (requiring compute). This test requires memory to hold and manipulate the rasters, e.g. 31 GB.
#
library(terra)


atol = 1e-6

neg_max_stage = rast('../threshold_given_exrate/ptha_highres_domains_neg_max_stage_at_epistemic_uncertainty_percentile_0.84_exrate_0.0004/sum_of_all_sources_neg_max_stage_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_12.vrt')

max_stage = rast('../threshold_given_exrate/ptha_highres_domains_max_stage_at_epistemic_uncertainty_percentile_0.84_exrate_0.0004/sum_of_all_sources_max_stage_rast_exrate_4e-04_percentile_84_subsam_1_Nrand_10000_seed_123.vrt')

neg_max_stage <- crop(neg_max_stage, max_stage, mask=TRUE)
residual <- abs(max_stage + neg_max_stage)
residual_range <- global(residual, fun = "range", na.rm = TRUE)

if (all(abs(range) < atol)) {
    print("PASS")
} else {
    print("FAIL")
}
