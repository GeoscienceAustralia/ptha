#
# Make a clean directory structure for the tifs generated by 
#    run_compute_exceedance_rates_for_threshold_depth_logic_tree_mean_newParallelPartition.sh
# AND 
#    all the *.sh scripts generated by `make_exceedance_rate_jobs.R` (epistemic uncertainty calcs over multiple nodes)
#

# This should match the beginning of the folders that will be moved to better locations
runname='ptha18-NSW2023-ID4186.3-sealevel110cm'

for variable in depth max_stage; do
  # Move logic-tree-mean hazard & variance rasters to a clean location
  mkdir -p ${runname}/highres_${variable}_with_variance ;
  mv ${runname}-${variable}-LogicTreeMean* ${runname}/highres_${variable}_with_variance;
  
  # Move epistemic uncertainty rasters to a clean location
  # 16th percentile
  mkdir -p ${runname}/highres_${variable}_epistemic_uncertainty/16pc; 
  mv ${runname}-${variable}_exrate_*_0.16_* ${runname}/highres_${variable}_epistemic_uncertainty/16pc/;
  # 84th percentile
  mkdir -p ${runname}/highres_${variable}_epistemic_uncertainty/84pc;
  mv ${runname}-${variable}_exrate_*_0.84_* ${runname}/highres_${variable}_epistemic_uncertainty/84pc/;
done
