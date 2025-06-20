# PTHA scenario nearshore testing

The simulations for this study were implemented over a long time period
(starting around 2020). Below I have condensed the process to run them. As I
received additional data during the study, some models were re-run with
extended spatial extents, which made some of the initially run models obsolete
(but the obselete runs are not removed from the scripts below).

The code was run on NCI and makes assumptions specific to that environment
(e.g. hard-coded locations of tide gauge data and local copy of `ptha` package)
that would need updating to run elsewhere.

# Notes on the runs

The study includes a few thousand model runs which were split into batches because
of limits in our supercomputer quota. 
* Initially I ran random ptha scenarios with 30-per-source-model-per-event. This was what I could mangage with the compute quota in that quarter. 
* In the subsequent quarter I had more quota and the opportunity to do another 30 runs for each of the previous events, inserting `batch2` in the filenames to avoid breaking things. 
* Later I added more historical tsunami events by running `batch3`  and `batch4`, using 60 scenarios per source model per event. They are stored in `OUTPUTS_2022_new_events`.
* Then I added yet another with `batch5`. The first versions of these runs were not used since the scenarios were later re-run in larger models with more data.
* In early 2025 I re-ran the events affecting WA with a greatly extended version of the model. I only did this for the Sumatra 2004 and Sumatra 2005 events (the South Sandwich 2021 event is only well resolved at the 1min Hillarys tide gauge, which is already covered in the previous model, so this was not rerun).
* Then I added some other historical events with `batch6`, as by this time I had the required data to model them.

To be clear, no tsunami scenarios were discarded in the above process -- but some simulations were updated to leverage more elevation data and tide gauge observations.

The batches used different base output directories, achieved by manually
changing the value of `output_basedir` (defined in `model_local_routines.f90`)
and recompiling. The following commented-out values correspond to random PTHA
scenarios (ignoring cases that were not used because they were later re-run
with more data).
```
!output_basedir = './OUTPUTS/' // &
!output_basedir = './OUTPUTS_SCRATCH/' // &
!output_basedir = './OUTPUTS_2022_new_events/' // &
!output_basedir = './OUTPUTS_2025_extend_WA/' // &
!output_basedir = './OUTPUTS_2025_NWWA/' // &
```
and (for the source inversions discussed below)
```
output_basedir = './OUTPUTS_new_validation_events/' // &
```

To manage file-counts (due to the inode quota on NCI) it is necessary to 'tar'
the multidomain directories, which we do with `run_tar_dir.sh` that calls a
worker R script (you need to make sure the pattern matching Sys.glob call is
finding the folders that you want to work on).

To work with the tar outputs, the easiest approach is to write a script that
extracts data required for analysis, which is small enough to not need tarring.
I did this with `run_extract_key_outputs_from_tarred_multidomain_3.sh` which
calls a worker R script that contains the details (you need to make sure that
the script parameters correspond with what you want to extract).

# HOW TO RUN THE MODELS

1. Create all random scenarios (for all batches) in the folder [../ptha18_scenarios_random](../ptha18_scenarios_random) by following instructions therein.

2. Do a test run for the SWALS model in this folder, and use the results to create a SWALS load balance file. 
  * For any given model run, a load balance file can be created by running the [load_balance_script.R](load_balance_script.R) from inside the multidomain directory (something like `Rscript ../../../load_balance_script.R` ). That will make a file `load_balance_partition.txt` in the folder where it was run, which is then copied to the `load_balance_files` folder with an appropriate name.
  * The load balance files were updated from time to time, as compiler versions were updated on NCI.
  * For more information on load balance files, see the [README for the global dissipation paper codes](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/nearshore_testing_2020/swals/README.md).

3. Create the PBS scripts that are required to run all jobs. 
  * In practice I used multiple batches, due to both limitations in compute quota, and the relatively long time over which the study was implemented.
  * Beware some of the scripts below run early versions of simulations that use less tide gauge data than the final version (e.g. early versions of the Sumatra simulations did not have such extensive tide gauge data in Western Australia, so were ultimately replaced with models that included more tide gauges). This is because we didn't have all the datasets when beginning the study. These cases are run again by the scripts - but you wouldn't do this if starting from scratch with all the data.
```
Rscript create_random_ptha_qsub_scripts.R
Rscript create_random_ptha_batch2_qsub_scripts.R
Rscript create_random_ptha_batch3_qsub_scripts.R
Rscript create_random_ptha_batch4_qsub_scripts.R
# Rscript create_random_ptha_batch5_qsub_scripts.R # No longer needed, runs included below
Rscript create_random_ptha_qsub_scripts_largerWAmodels.R
Rscript create_random_ptha_batch6_qsub_scripts.R
```

4. Submit all the jobs to the PBS queue, and wait until they finish.
  * In practice they were run in batches, with the code recompiled for each batch to change the `output_basedir` as mentioned above.
  * Once each batch is finished, tar the directories with [run_tar_dir.sh](run_tar_dir.sh) to avoid exceeding the file-count qutota. 
    * To do this you'll need to edit the variable `all_md_dir` in [tar_multidomain_dirs.R](tar_multidomain_dirs.R) to ensure it globs the correct files.

5. Run the test events derived from source inversions. For all these runs the code was compiled with `output_basedir = './OUTPUTS_new_validation_events/' // &`
```
# This script collates all the runs into a single source for convenience, but
# it's an over-simplification of the process in which they were often run in
# small batches along with other tests and sometimes output to different
# directories (but collated later). However this will not change the
# simulations (aside from perhaps names of some dirs in output files).
qsub run_inversions_for_each_event.sh
```

Then use the scripts in `./plots/plot_XXXX.R` to move the gauge data to an RDS file, ready for later ploting.


6. Extract tide gauge signals from the simulations using
```
qsub run_extract_key_outputs_from_tarred_multidomain_3.sh
```
It actually needs to be run twice, editing the following variable in
`extract_key_outputs_from_tarred_multidomain_3.R` to be `TRUE` and `FALSE`
respectively.
```
PROCESS_RANDOM_SCENARIOS = FALSE # TRUE for random scenarios, FALSE for validation events 
```

7. At this point the processing moves to an analysis directory somewhere in `../`
