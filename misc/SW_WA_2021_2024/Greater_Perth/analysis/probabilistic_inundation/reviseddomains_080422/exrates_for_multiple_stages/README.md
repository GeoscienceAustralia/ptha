# Calculation of 'binned tsunami maxima' satisfying a particular exceedance-rate.
---------------------------------------------------------------------------------

We compute the tsunami maxima (above the ambient sea-level) for a specified exceedance-rate (based on the logic-tree mean model).

For computational efficiency, the calculations bin the max-stage values to 
thresholds of "ambient sea level + [0.001, 1, 2, 3, ... 10]" meters. This means
that the tsunami maxima results are rounded down to 0, 1, 2, 3, ... 10m above the ambient sea-level.

## Key files and folders:

* The input datasets are in folders beginning with `ptha18`, and depict exceedance-rate maps for the different max-stage thresholds. They were created with the script in [../../run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh](../../run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh). 
* To combine the unsegmented and segmented model results we use the script [compute_exrates_multiple_stages.R](compute_exrates_multiple_stages.R). 
    * Run with `Rscript compute_exrates_with_multiple_stages.R`. 
    * This will make folders like `summed_results_XX` where `xx` is the maximum stage. These contain exceedance-rate rasters
* Once the above have been run, the script [compute_banded_stages.R](compute_banded_stages.R) can be used to compute the tsunami maxima (above ambient sea-level) corresponding to the given exceedance-rate. 
    * In this script we specify the exceedance-rate of interest. 
    * Run with `Rscript compute_banded_stages.R`
* Once the above has been run, the script [make_vrt.sh](make_vrt.sh) is used to make a VRT. Run from INSIDE THE FOLDER CONTAINING THE TIFS with `source ../make_vrt.sh`. 

