# Run this from INSIDE the folders containing the tifs
export RESVALS=' -resolution highest '

for namestart in unsegmented_HS sumatra_segment_HS java_segment_HS andaman_segment_HS arakan_segment_HS; do 
    gdalbuildvrt $RESVALS $namestart'_exrate_all.vrt' $namestart'_domain_'*domain__exceedance*.tif;
    gdalbuildvrt $RESVALS $namestart'_exrate_variance_all.vrt' $namestart'_domain_'*domain__variance*.tif;
    done
