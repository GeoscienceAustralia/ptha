## Exceedance-rate calculations

Below we show how to compute rasters depicting:
1. A: The logic-tree-mean rate of inundation (depth > 1mm).
2. B: The logic-tree-mean exceedance-rate for a range of max-stage values (e.g. 1.1, 2.1, 3.1, ..., 10.1)
3. The rate of inundation (depth > 1mm) at the 16th and 84th percentile epistemic uncertainty
4. Organise files and sum above results over source zones 
5. Compress outputs


## 0. Modify application_specific_inputs.R for your case

## 1. Inundation exceedance-rate calculations - Logic-tree-mean case 1st step
### A Depth
Modify the main run script to suit your case, then run it
```bash
qsub compute_exceedance_rates_for_threshold_depth_logic_tree_mean.pbs
```

### B A range of max-stage values
Modify the main run script to suit your case, then run it
```bash
qsub compute_exceedance_rates_for_threshold_max_stages_logic_tree_mean.pbs
```
NOTE: The above script may run out of memory when using the CascadeLake nodes.
An alternative version on the SapphireRapids nodes is also provided.
If you use this, be sure to first change MC_CORES=104 in [../application_specific_file_metadata.R](../application_specific_file_metadata.R).
Then change it back after the script has finished (since the other scripts below use CascadeLake).

## 2. Inundation exceedance-rate calculations (epistemic uncertainty case)
### Step 2a
Check and edit make_exceedance_rate_jobs.R and the associated template script.
These will make a set of qsub files to run the calculations
Make the qsub files with:
```bash
Rscript make_exceedance_rate_jobs.R
```

### Step 2b. Check a few of those qsub files (in case you made a mistake in setup)
If they are OK then submit them, moving to a folder after submission.

### Step 2c. Submit jobs to compute exceedance rates at epistemic uncertainty
```bash
# mkdir -p submitted_exceedance_rate_at_epistemic_jobs
for i in compute_exceedance_rates_at_epi*.pbs; do echo $i; qsub $i; mv $i submitted_exceedance_rate_at_epistemic_jobs/; done
```

## 3. Organise files
At this point the outputs are inside new folders within the current directory.
It's a bit messy and should be better organised.
Move them into organised sub-folders by editing "make_directory_structure.sh"
for your case. Then run it.
```bash
source make_directory_structure.sh
```
You should end up with a folder structure like this:
- ptha/
  - highres_depth_with_variance/  <!-- DEPTH, LOGIC TREE MEAN RESULTS -->
    - ptha-depth-LogicTreeMean-outerrisesunda/
    - ptha-depth-LogicTreeMean-sunda2/
    - ... other source zones if present ...
  - highres_max_stage_with_variance/  <!-- MAX_STAGE, LOGIC TREE MEAN RESULTS -->
    - ptha-max_stage-LogicTreeMean-outerrisesunda/
    - ptha-max_stage-LogicTreeMean-sunda2/
    - ... other source zones if present ...
  - highres_depth_epistemic_uncertainty/  <!-- DEPTH, EPISTEMIC UNCERTAINTY RESULTS -->
    - 84pc/
      - ptha-depth_exrate_0.001_0.84_outerrisesunda/
      - ptha-depth_exrate_0.001_0.84_sunda2/
      - ... other source zones if present ...
    - 16pc/
      - ptha-depth_exrate_0.001_0.16_outerrisesunda/
      - ptha-depth_exrate_0.001_0.16_sunda2/
      - ... other source zones if present ...
  - EMPTY FOLDERS FOR MAX STAGE EPISTEMIC UNCERTAINTIES

## 4.  Use calculations above to make the logic-tree-mean inundation exceedance-rate and variance, summed over sources.
The command line argument gives the path to the logic-tree-mean results above,
the variable (depth) and the exceedance-threshold (0.001). **This uses parallel computing. Get compute.**
```bash
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_stage_with_variance/ max_stage 0.001 &
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_speed_with_variance/ max_speed 1.543 &
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_speed_with_variance/ max_speed 3.087 &
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_speed_with_variance/ max_speed 4.630
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_depth_with_variance/ depth 0.001
```
This created a folder inside the 'highres_depth_with_variance' sub-folder above, containing
results summed over source-zones (exceedance-rate, upper 95% CI for true exceedance-rate, variance). 
The folder name is like:
./ptha/sea_level_vary/highres_depth_with_variance/ptha/sea_level_vary-depth-LogicTreeMean-sum_of_source_zones

Use calculations above to get the logic-tree mean exceedance-rates for all the alternative max-stage thresholds. The thresholds must be the same as
specified in the earlier qsub script. **This uses parallel computing. Get compute.**
``` bash
# for stage_threshold in 1.101 2.1 3.1 4.1 5.1 6.1 7.1 8.1 9.1 10.1 ; do
for stage_threshold in 0.001 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 ; do
    echo $stage_threshold ;
    Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_stage_with_variance max_stage $stage_threshold ;
done
```
Rscript compute_mean_exrate_upper_CI.R ptha/sea_level_vary/highres_max_stage_with_variance max_stage 0.001

## 5. Use calculations above to get the inundation exceedance-rates at different epistemic uncertainty percentiles, summed over source-zones.
This assumes co-monotonic epistemic uncertainties between the sources (conservative).
**Runs in parallel. Get compute.**
``` bash
Rscript compute_sum_of_percentiles.R ptha/sea_level_vary/highres_depth_epistemic_uncertainty/ 84 depth 0.001
Rscript compute_sum_of_percentiles.R ptha/sea_level_vary/highres_depth_epistemic_uncertainty/ 16 depth 0.001
```
Rscript compute_sum_of_percentiles.R ptha/sea_level_vary/highres_max_speed_epistemic_uncertainty/ 16 max_speed 1.543 2>&1 | tee compute_sum_of_percentiles_speed_16_1.543.log
Rscript compute_sum_of_percentiles.R ptha/sea_level_vary/highres_max_speed_epistemic_uncertainty/ 84 max_speed 1.543 > compute_sum_of_percentiles_speed_84_1.543.log

This created folders containing sums over source zones, with names like:
sea_level_vary/highres_depth_epistemic_uncertainty/84pc/sea_level_vary-depth_exrate_0.001_0.84_sum_of_source_zones/

## 6. Compress the output folders
At this point (or before) you should compress the output folders as they can contain many many files.
