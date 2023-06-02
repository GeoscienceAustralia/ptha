# Run this from INSIDE the folders containing the tifs
export RESVALS='0.000102880658436 0.000102880658436'

gdalbuildvrt -tr $RESVALS max_depth_binned_values_with_exrate_exceeding_0.0001.vrt *.tif;
