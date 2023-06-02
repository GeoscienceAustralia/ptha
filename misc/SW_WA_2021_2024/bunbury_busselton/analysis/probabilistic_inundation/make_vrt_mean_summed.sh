# Run this from INSIDE the folders containing the tifs
export RESVALS='0.000102880658436 0.000102880658436'

for namestart in summed_HS; do 
    gdalbuildvrt -tr $RESVALS $namestart'_exrate_all.vrt' $namestart'_domain_'*domain__exceedance*.tif;
    done
