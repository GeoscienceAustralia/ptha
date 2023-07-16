# Calculation of 'binned tsunami maxima' satisfying a particular exceedance-rate.
---------------------------------------------------------------------------------

We compute the tsunami maxima (as height above the ambient sea-level) for a
specified exceedance-rate, based on the logic-tree mean model. 

The main outputs are contained in folders named like
[binned_max_stage_exceeding_4e-04](binned_max_stage_exceeding_4e-04). Here the
number `4e-04 (=1/2500)` is the exceedance-rate. The folder contains rasters
taking the values 0.001, 1, 2, 3, ... 10m, which give the max-stage ABOVE THE
BACKGROUND SEA-LEVEL of 0.6m AHD, at the 1/2500 exceedance-rate, rounded down to the nearest
integer (with an upper limit of 10, and a lower limit of 1mm). 

Here use rounded-down max-stage values for computational efficiency reasons. It
would be possible to provide results with greater resolution of the max-stage
(just less convenient).

We include outputs for exceedance-rates of:
* 1/100 - [binned_max_stage_exceeding_0.01](binned_max_stage_exceeding_0.01)
* 1/500 - [binned_max_stage_exceeding_0.002](binned_max_stage_exceeding_0.002)
* 1/2500 - [binned_max_stage_exceeding_4e-04](binned_max_stage_exceeding_4e-04)
* 1/10000 - [binned_max_stage_exceeding_1e-04](binned_max_stage_exceeding_1e-04)


## Key files and folders:

* The input datasets are in folders beginning with `ptha18`, and depict exceedance-rate maps for the different max-stage thresholds. They were created with the script in [../../run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh](../../run_compute_exceedance_rates_for_threshold_stage_logic_tree_mean_newParallelPartition.sh). 
* To combine the unsegmented and segmented model results we use the script [compute_exrates_multiple_stages.R](compute_exrates_multiple_stages.R). 
    * Run with `Rscript compute_exrates_with_multiple_stages.R`. 
    * This will make folders like `summed_results_XX` where `xx` is the maximum stage. These contain exceedance-rate rasters
* Once the above have been run, the script [compute_banded_stages.R](compute_banded_stages.R) can be used to compute the tsunami maxima (above ambient sea-level) corresponding to the given exceedance-rate. 
    * In this script we specify the exceedance-rate of interest (EDIT ACCORDINGLY)
    * Run with `Rscript compute_banded_stages.R`
* Once the above has been run, the script [make_vrt.sh](make_vrt.sh) is used to make a VRT. You'll need to edit the filename for the output vrt. Run from INSIDE THE FOLDER CONTAINING THE TIFS with `source ../make_vrt.sh`. 

