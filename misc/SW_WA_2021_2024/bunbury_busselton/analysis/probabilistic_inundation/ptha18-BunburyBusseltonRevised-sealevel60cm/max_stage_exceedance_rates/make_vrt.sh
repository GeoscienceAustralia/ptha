# Run this from INSIDE the folders containing the tifs

gdalbuildvrt -resolution highest max_stage_binned_values_with_exrate_exceeding_0.01.vrt *.tif;
#gdalbuildvrt -resolution highest max_stage_binned_values_with_exrate_exceeding_0.002.vrt *.tif;
#gdalbuildvrt -resolution highest max_stage_binned_values_with_exrate_exceeding_0.0004.vrt *.tif;
#gdalbuildvrt -resolution highest max_stage_binned_values_with_exrate_exceeding_0.0001.vrt *.tif;
