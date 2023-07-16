runname='ptha18-BunburyBusseltonRevised-sealevel60cm'

# Move logic-tree-mean hazard & variance rasters to a clean location
mkdir -p $runname/max_stage_exceedance_rates
mv $runname-random* $runname/max_stage_exceedance_rates    

