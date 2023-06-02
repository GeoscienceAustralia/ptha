# Calculation of 'binned depths' satisfying a particular exceedance-rate.
---------------------------------------------------------------------------------

We compute the tsunami depth maxima for a specified exceedance-rate (based on the logic-tree mean model).

For computational efficiency, the calculations bin the depth values to 
thresholds of `c(0.001, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10)` meters.
The tsunami depth maxima results are rounded down to these values.

## Key files and folders:

* The input datasets are in folders beginning with `ptha18`, and depict exceedance-rate maps for the different depth thresholds. They were created with the script in [../../run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh](../run_compute_exceedance_rates_for_multiple_threshold_depths_logic_tree_mean_newParallelPartition.sh) and then manually moved to this location.
* To combine the unsegmented and segmented model results we use the script [compute_exrates_multiple_depths.R](compute_exrates_multiple_depths.R). 
    * Run with `Rscript compute_exrates_with_multiple_depths.R`. 
    * This will make folders like `summed_results_XX` where `xx` is the maximum depth. These contain exceedance-rate rasters
* Once the above have been run, the script [compute_banded_depths.R](compute_banded_depths.R) can be used to compute the depth corresponding to the given exceedance-rate (discretized to a set of values by rounding down). 
    * In this script we specify the exceedance-rate of interest. 
    * Run with `Rscript compute_banded_depths.R`
* Once the above has been run, the script [make_vrt.sh](make_vrt.sh) is used to make a VRT. Run from INSIDE THE FOLDER CONTAINING THE TIFS with `source ../make_vrt.sh`. 

