# Set these to point to the shifted GEBCO file, and the original GA250 dem
GEBCO_SHIFTED_1m='../GEBCO_2014_1m/GEBCO_2014_1minx1min_W-39.9958333-E320.0041667.tif'
GA250_dem='../../../../../DATA/ELEVATION/GA250m/ER_Mapper_ers/ausbath_09_v4_ex_ex.ers'

# Convert GA250 to wgs84, and down-sample to 1-arc-minute (1/60 = 0.0166666..)
# using bilinear interpolation
echo 'Warping....'
gdalwarp -t_srs epsg:4326 -te 92 -60 171 -9 -tr 0.01666667 0.01666667 -r bilinear $GA250_dem 'ga250_wgs84_newextent.tif'

# Merge the two DEMS
echo 'Merging....'
gdal_merge.py -o 'merged_gebco_ga250_dem.tif' -co "COMPRESS=DEFLATE" $GEBCO_SHIFTED_1m 'ga250_wgs84_newextent.tif'

# PATCH some obviously erronious sections of GA250 #
Rscript patch_dem.R
