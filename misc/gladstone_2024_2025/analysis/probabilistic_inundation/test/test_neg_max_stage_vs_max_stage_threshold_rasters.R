#
# Compare the max-stage thresholds at given exceedance-rate against negative max-stage's. The rasters should be identical but flipped in sign. It requires the hard-coded inputs below to be already generated (requiring compute). This test requires memory to hold and manipulate the rasters, e.g. 31 GB.
#
library(terra)


atol = 1e-6

# Before running this script you'll need to generate the max stage and negative max stage virtual rasters.
# Follow along in the ../threshold_given_exrate/README.md. Generate the pbs scripts, run them and then tidy the bounds.
neg_max_stage = rast('../threshold_given_exrate/ptha_max_stage_1in2500_84pc/ptha_neg_max_stage_1in2500_84pc.vrt')

max_stage = rast('../analysis/probabilistic_inundation/threshold_given_exrate/ptha_max_stage_1in2500_84pc/ptha_max_stage_1in2500_84pc.vrt')

neg_max_stage <- crop(neg_max_stage, max_stage, mask=TRUE)
residual <- abs(max_stage + neg_max_stage)
residual_range <- global(residual, fun = "range", na.rm = TRUE)

if (all(abs(range) < atol)) {
    print("PASS")
} else {
    print("FAIL")
}
