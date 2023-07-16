# Run this from INSIDE the folders containing the tifs

# The percentile is stored in the name of the directory 2-levels above
export percentile_tag=$(basename $(dirname $(pwd)))
echo $percentile_tag
gdalbuildvrt -resolution highest "exceedance_rate_percentile_"$percentile_tag".vrt" *.tif;
