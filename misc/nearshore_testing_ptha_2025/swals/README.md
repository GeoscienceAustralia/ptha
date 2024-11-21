# PTHA scenario nearshore testing

The following codes were used to create qsub files to run random PTHA18 scenarios with location similar to historic events:

    Rscript create_random_ptha_qsub_scripts.R	
    Rscript create_random_ptha_batch2_qsub_scripts.R	
    Rscript create_random_ptha_batch3_qsub_scripts.R	
    Rscript create_random_ptha_batch4_qsub_scripts.R	
    Rscript create_random_ptha_batch5_qsub_scripts.R	

These were split into batches because limits in the compute quota meant this study was undertaken over a multi-year time period. The batches used different base output directories, achieved by manually changing the value of `output_basedir` (coded around line 156 of `model_local_routines.f90`) and recompiling. The following commented-out values correspond to random PTHA scenarios.
```
!output_basedir = './OUTPUTS/' // &
!output_basedir = './OUTPUTS_SCRATCH/' // &
!output_basedir = './OUTPUTS_2022_new_events/' // &
!output_basedir = './OUTPUTS_2023_new_events/' // &
```

To manage file-counts (due to the inode quota on NCI) it is necessary to 'tar' the multidomain directories, which we do with `run_tar_dir.sh` that calls a worker R script (you need to make sure the pattern matching Sys.glob call is finding the folders that you want to work on).

To work with the tar outputs, the easiest approach is to write a script that extracts data required for analysis, which is small enough to not need tarring. I did this with `run_extract_key_outputs_from_tarred_multidomain_2.sh` which calls a worker R script that contains the details (you need to make sure that the script parameters correspond with what you want to extract).

# Notes on the runs

There were multiple 'batches' done over different quarters.

* Initially I ran random ptha scenarios with 30-per-source-model-per-event. This was what I could mangage with the compute quota in that quarter. 
* In the subsequent quarter I had more quota and the opportunity to run more events. 
* Some care was required to do this without breaking things:
    * I had to be careful to not 'contaminate' the file structure of the previous runs, and also, to not overload the file-system.
    * Recall the model output folder name includes information on the 'scenario count' [i.e. number of times each scenario was sampled -- usually 1, but sampling with replacement can lead to more]. If I run more scenarios, I need to be careful to not end up with multiple model runs in the same folder [which will happen if the names end up matching -- in which case the 'count' field in the folder name would be wrong. Alternatively I could manually modify the existing folder names, or do corrections within my script, but both seem error prone]. 
    * The solution I implemented is to change the names of the input rasters -- make them in a neighbouring folder that includes "-batch2" in the scenario name. Because this name appears in the output folder, this ensures there is no output-folder name clash. 
    * To manage disk space, it is better to write to "scratch" rather than gdata, because we have more space and more inodes. To do this we change `OUTPUTS` to `OUTPUTS_SCRATCH` in the model output folder, and make it a symbolic link to `/scratch/w85/gxd547/...somewhere-good....`
        * Later I had to move `OUTPUTS_SCRATCH` back onto gdata, to avoid automatic deletion of the files.
    * In the first batch of runs I saved grid-outputs 3 times per simulation. This inflates the output file size more than needed (I haven't been using those outputs at all), so for the second batch I `changed the code to only write the initial condition grids`.
* Later I added more historical tsunami events by running `batch3`  and `batch4`. They are stored in `OUTPUTS_2022_new_events`.
* Then I added yet another with `batch5`. This is stored in `OUTPUTS_2023_new_events`.

# HOW TO RUN THE MODELS

1. Create all random scenarios (for all batches) in the folder [../ptha_scenarios_random](../ptha_scenarios_random) by following instructions therein.

2. Do a test run for the SWALS model in this folder, and use the results to create a SWALS load balance file. 
  * In practice, the load balance files were originally created using runs for [the global dissipation paper](https://doi.org/10.3389/feart.2020.598235), which uses the same model setup as herein.
    * The load balance files were updated from time to time, as compiler versions were updated on NCI.
    * For any given model run, a load balance file can be created by running the [load_balance_script.R](load_balance_script.R) from inside the multidomain directory (something like `Rscript ../../../load_balance_script.R` ). That will make a file `load_balance_partition.txt` in the folder where it was run, which is then copied to the `load_balance_files` folder with an appropriate name.
    * For more information on load balance files, see the [README for the global dissipation paper codes](https://github.com/GeoscienceAustralia/ptha/tree/master/misc/nearshore_testing_2020/swals/README.md).

3. Create the PBS scripts that are required to run all jobs. 
  * In practice I used multiple batches, due to both limitations in compute quota, and the relatively long time over which the study was implemented.
```
Rscript create_qsub_commands.R
Rscript create_random_ptha_batch2_qsub_scripts.R
Rscript create_random_ptha_batch3_qsub_scripts.R
Rscript create_random_ptha_batch4_qsub_scripts.R
Rscript create_random_ptha_batch5_qsub_scripts.R
```

4. Submit all the jobs to the PBS queue, and wait until they finish.
  * In practice they should be run in batches.
  * Once each batch is finished, tar the directories with [run_tar_dir.sh](run_tar_dir.sh) to avoid exceeding the file-count qutota. 
    * To do this you'll need to edit the variable `all_md_dir` in [tar_multidomain_dirs.R](tar_multidomain_dirs.R) to ensure it globs the correct files.
