# Run this from INSIDE the folders containing the tifs

gdalbuildvrt -resolution highest 'summed_HS_exrate_all.vrt' 'summed_HS_domain_'*domain__exceedance*.tif;
gdalbuildvrt -resolution highest 'summed_HS_exrate_variance_all.vrt' 'summed_HS_domain_'*domain__variance*.tif;
gdalbuildvrt -resolution highest 'summed_HS_exrate_upperMonteCarloCI_all.vrt' 'summed_HS_domain_'*domain__Monte*.tif;
