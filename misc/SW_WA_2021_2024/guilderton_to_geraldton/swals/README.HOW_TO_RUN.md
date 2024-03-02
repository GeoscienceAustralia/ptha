# Outline of steps to run on NCI

## COMPILE THE MODEL
Starting grom this (swals) directory
``` bash
# Get the right modules
module clear
source modules_SWALS_ifort_2023_B.sh
# Compile
make -B -f make_model_ifort_sapphirerapids
```

## SANITY CHECKS OF AN INITIAL MODEL
```bash
# Run a short model
qsub run/test_model_sapphirerapids.pbs

# Go to the multidomain output folder (here the folder name is indicative)
cd ./OUTPUTS/short_test_run-test_load_balance-ambient_sea_level_0.6/RUN_20230830_135916498
# and open a log file, e.g., multidomain_log_image_00000000000000000001.log
#   Did the model finish? (should end with printing of timers)
#   Are the printed namelist variables correct? (grep for &MULTIDOMAIN)
#   Are any domains unstable with the default timestepping, on any MPI process? [The following should be empty]
grep UNSTABLE multi*.log
#   Did any model produce NaN, on any MPI process? [The following should be empty]
grep NaN multi*.log
```

## Post-process the simulation with R
```bash
# Here we are still inside the model output directory (as above)
module clear
source ../../../modules_R_431.sh

# Create & check the domain shapefile (load R modules as needed)
Rscript ../../../post_process/make_domains_shapefile.R "."
# The above command created a directory "domains_shapefile" in the working folder.
# scp the files to your computer and check in GIS -- are the domains where you expected them to be?

# Make rasters with the elevation (here assuming 24 cores -- which
# needs an interactive job on NCI. Too many cores might run out of memory).
Rscript ../../../post_process/make_rasters.R '.' 24 elevation0
# The above should make many raster (tifs) in the working directory with names matching
# elevation0*.tif. 
# Download them, then combine into a vrt (to treat as a single file) with
gdalbuildvrt -resolution highest elevation_all.vrt elevation0*.tif
```
Look at elevation_all.vrt in GIS -- does it look right? Is the accurate, high
resolution data where you expect it do be? Did you get the file order right?
Are there unexpected artefacts? (e.g. caused by the interpolation used visually
in GIS not matching the bilinear interpolation in SWALS).

## CONSIDER AN ALTERNAITVE global_dt
Check the allowed timestep on each domain. This is in the logs, near where the model writes
"stable" (if the timestepping seems stable) or UNSTABLE (but that was ruled out above). It
can be found with:
```bash
grep -B1 stable *.log 
# To focus just on the global domains, which have domain index 1 (in this case):
grep -B1 stable *.log | grep "0000001 ,"
```
If the minimum global domain timestep is much larger than the chosen global_dt, you might
consider increasing global_dt. However that will fatten the halo regions, so there is
a trade-off. You can test by comparing "load balanced models" with alternative global_dt.

## RUN A LOAD BALANCED MODEL
If everything looks sensible then make a load balance file from within the
multidomain directory
```bash
Rscript ../../../post_process/load_balance_script.R
# This will print information on the expected improvement.

# If the speed up looks good, then store the file and run a new load balanced model.
mkdir -p ../../../load_balance_files
cp load_balance_partition.txt ../../../load_balance_files/load_balance_16MPI_NNL4_defaultres.txt
# Go back to the "swals" directory
cd ../../../
# Make a new model configuration namelist that is similar to the previous one, but uses the newly 
# created load balance file
cp multidomain_design_control_NNL4_defaultres_defaultloadbalance.nml multidomain_design_control_NNL4_defaultres.nml 
code multidomain_design_control_NNL4_defaultres.nml # (EDIT the load balance file)
# Make a run script that uses the new load balance file.
code run/test_model_with_load_balance_sapphirerapids.pbs # (EDIT to use the new multidomain_design....nml)
# Run a short load balanced model and check that the speedup is OK. 
qsub run/test_model_with_load_balance_sapphirerapids.pbs

# If the speedup isn't good, or a good speedup isn't possible according to the load_balance_script.R, then
# go to the multidomain output directory and check the time spent in "evolve:".
cd name_of_multidomain_folder
grep " evolve:" *.log
# (Ignore the overall times which can be strongly affected by the variable disk speed when reading inputs).
# Compare this to the time in parallel communication -- ideally the latter
# should be nearly equal between processes, and small compared to the evolve
# time.
grep "comms1" *.log
# Compute the average of these times using: 
grep " comms1:" *.log | awk 'NF>1{print $(NF-1)}' | awk '{sum += $1; count++} END {if (count > 0) print sum/count}'
# If the communication is quite unequal between processes, check the model for
# individual domains with long run times (comparable to the evolve time you can
# achieve) which can prevent load balancing.
Rscript ../../../post_process/report_domain_runtimes.R
```

Any problematic domains might be split or coarsened or removed, to reduce their run time and
allow balancing the model.

If the time in "comms1" is quite equal, but seems too large, consider running with fewer MPI
processes. [Or, redesigning the model to reduce the halo fatness].

## Run a longer model
Below here I assume the model design matches the intent, and it is reasonably efficient.
Read the scripts to understand/modify for your case.

Go back to the "swals" directory from the model run folder.
```bash
cd ../../..
```

Run some validation tests, and low-resolution variants to check convergence.
In the case with a time-varying source, the source file was created with
```bash
Rscript pre_process/make_initial_conditions_complex_historical_events.R
# Then you can run the historical events
qsub run/Sumatra2004_time_varying.pbs
qsub run/Sumatra2004_time_varying_lowres.pbs
qsub run/Sumatra2005.pbs
qsub run/Sumatra2005_lowres.pbs

# Run an extreme source, and a low-resolution version to check convergence in the inundation regime.
qsub run/extreme_source.pbs
qsub run/extreme_source_lowres.pbs

# Run a very small source -- this can be useful to discover nesting artefacts.
# (especially max_flux).
qsub run/small_source.pbs
```

## After running the above jobs
You might want to make a folder to hold the *.sh.* files resulting from finished runs.
* I usually call this dot_e_o_files

 Check the load balancing again. 
 * The efficiency of regular resolution models may be improved by making a new load_balance_file based on one of those runs. 
* The efficiency of low resolution runs might also be improved (with a separate load-balance file). 

Check the log files for mass conservation, energy decay, reasonable peak stage/speed etc.
See ../analysis/check_log_files for more info.

Make rasters in the usual way and check the results, e.g. on a node with 48 cores from inside the multidomain_dir
```bash
Rscript ../../../post_process/make_rasters.R '.' 48 max_stage max_flux max_speed time_of_max_stage
gdalbuildvrt -resolution highest max_stage_all.vrt max_stage_domain*.tif
gdalbuildvrt -resolution highest max_speed_all.vrt  max_speed_domain_*.tif
gdalbuildvrt -resolution highest time_of_max_stage_all.vrt  time_of_max_stage_domain_*.tif
gdalbuildvrt -resolution highest max_flux_all.vrt  max_flux_domain_*.tif
# If you saved temporal parameters you could additionally create the last timestep rasters
gdalbuildvrt -resolution highest last_timestep_UH_all.vrt last_time_UH_domain_*.tif
gdalbuildvrt -resolution highest last_timestep_VH_all.vrt last_time_VH_domain_*.tif
#
# Process gauges for the historical validation cases (two steps for each event, figures end up in the ./plots/ folder).
# Example for Sumatra 2004
cd plots
Rscript process_gauges_sumatra2004.R ../OUTPUTS/path_to_multidomain_dir_sumatra2004
cd ..
# plot with varying y-scale
Rscript plots/plot_gauges_perth_sumatra2004.R OUTPUTS/path_to_multidomain_dir_sumatra2004/gauge_RDS_file_name.RDS
# uniform y range of 1m
Rscript plots/plot_gauges_perth_sumatra2004.R OUTPUTS/path_to_multidomain_dir_sumatra2004/gauge_RDS_file_name.RDS 1.0

# As above, Sumatra 2005
cd plots
Rscript plots/process_gauges_sumatra2005.R ../OUTPUTS/path_to_multidomain_dir_sumatra2005
cd ..
# varying y range
Rscript plots/plot_gauges_perth_sumatra2005.R OUTPUTS/path_to_multidomain_dir_sumatra2005/gauge_RDS_file_name.RDS
# uniform y range of 0.3m
Rscript plots/plot_gauges_perth_sumatra2005.R OUTPUTS/path_to_multidomain_dir_sumatra2005/gauge_RDS_file_name.RDS 0.3
```


## Run the random scenarios
```bash
# Set-up a script that can create all the PBS submission scripts, for all
# scenarios you want to run.
cd pre_process
vim create_random_ptha_qsub_scripts_sealevel60cm.R
# Then make the scripts.
Rscript create_random_ptha_qsub_scripts_sealevel60cm.R
# Eyeball a few scripts. 
# If they look OK, then submit one from the swals directory to see if it works.
cd ..
qsub run/midwest_sealevel60cm_1001.pbs
# It's probably best to move the submitted scripts to another directory, to keep
# track of which have been submitted.
mkdir run/midwest_submitted_jobs
mv run/midwest_sealevel60cm_1001.pbs run/midwest_submitted_jobs/

# Once you're convinced that it's running correctly (including tarring the
# output folders) you can submit all jobs with e.g.
for i in run/midwest_sealevel60cm_10*.pbs; do echo $i; qsub $i; mv $i run/midwest_submitted_jobs/; done
```
Be careful as this typically uses lots of compute.
* Avoid overwhelming the queue, and ensure you have enough disk/inodes/KSU.
* In practice you might want to submit in batches of 10 (or whatever) to deal with those limits, and double-check everything is OK without using all your compute quota.
* Some runs might fail, and they're rerun using [run/midwest_redo4.pbs](run/midwest_redo4.pbs) for example. The [workflow notes](../WORKFLOW.md) has more details on this.


## After the random scenarios have finished
```bash
# Check that they all completed, screen the log files for weird results
cd OUTPUTS/ptha18-midwest-sealevel60cm
# Scan every log file for NaN values -- ideally this will find no matches
grep NaN */*/multi*00001.log
# Identify which runs produced a NaN
grep NaN */*/multi*00001.log | cut -d " " -f1 | uniq

# Look at the second last line of each run. If it finished, then it should 
# have a timer printout, something like 
#     "Total WALLCLOCK time:       10555.17474110"
for i in */*/multi*00001.log; do tail -n2 $i | head -n1 ; done

# Run the log file checks, specifing the ambient sea level (0.6m in this case) 
# and the number of hours the model runs for (24 in this case)
# and a string that will match one log file for each run (via Sys.glob)
cd ../../../../analysis/check_log_files
Rscript check_log_files.R 0.6 24 "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/*/*/multi*00001.log"
```

Look at the resulting pdf plot, and the summary statistics (example below)
- The statistics should show that all runs finished, had good mass
  conservation, had energy that did not increase significantly from the initial
  value (although tiny increases in the first few steps are common and
  acceptable for models with a static source, while if the source is
  time-varying then there may be large energy increases over time as the source
  is applied).
- The plot should show the energy decaying over time.
  - Tiny increases in the first few steps are acceptable
  - Once flows interact with the boundary we could "in-principle" have energy
    increases, although in practice don't expect that.
- The kinetic energy should also show a general decay (albeit with some bumps etc) 
- The max-stage and max-speed can vary, but beware crazy oscillations or very
  extreme values.

```
[gxd547@gadi-login-04 check_log_files]$ Rscript check_log_files.R 0.6 24 "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/*/*/multi*00001.log"
[1] "Did the model runs finish?"
   Mode    TRUE 
logical     369 
[1] "Mass conservation errors relative to initial volume"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
1.644e-16 1.874e-16 1.977e-16 1.992e-16 2.093e-16 2.576e-16 
[1] "Mass conservation errors relative to boundary flux integral"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
7.220e-11 1.827e-10 3.658e-10 2.077e-08 1.827e-09 1.216e-06 
[1] "(Maximum energy - initial energy) relative to (2x maximum kinetic energy), BEFORE BOUNDARY FLUXES"
[1] "(typically very small unless the source is time-varying)"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000000 0.0002719 0.0007013 0.0006621 0.0009920 0.0015067 
[1] "Time index with largest kinetic energy (usually near start)"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   3.000   3.477   4.000  11.000 
[1] "Maximum energy increase between timesteps, relative to (2x maximum kinetic energy)"
[1] "(typically very small unless the source is time-varying)"
      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
-0.0004571  0.0003264  0.0007534  0.0007239  0.0010201  0.0038118 
```

At this point you might need to further investigate and fix any problematic runs.

Before proceeding to probabilistic computations, the fixing also needs to: 
* Ensure there is only one model simulation per scenario under the
  ./OUTPUTS/name_of_this_batch_of_runs/ folder.
* Ensure that all models have the same domain structure / resolutions / halo
  size etc
  * Or, edit the subsequent processing scripts so that they tolerate any
    differences you introduce.
    * This might not be easy, so try to avoid it.


## Once you think all the runs are OK, make the rasters

```bash
qsub post_process/create_tarred_rasters_from_tarred_multidomains.sh

# Check the max-stage exceedance-rate curves vs PTHA18
# * If these aren't good, consider doing more runs. Multiple Importance Sampling provides the theory for combining several runs.
cd ../analysis/max_stage_at_a_point/
# Edit extract_max_stage_at_a_point.R so all the input parameters are correct
# Maybe edit run_a_few.sh
qsub run_a_few.sh
```
