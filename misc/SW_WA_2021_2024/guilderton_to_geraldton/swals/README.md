# Tsunami model runs for Midwest region of WA, FY2023-2024
--------------------------------------------------

This folder contains code to run the tsunami models, including:
* Scenarios like historical events [(see here)](../../Greater_Perth/sources/like_historic/)
* Random PTHA18 scenarios [(see here)](../../Greater_Perth/sources/hazard/)

The tsunami model includes the area from Seabird to Geraldton in high resolution.
It is similar to earlier work on the Greater Perth area, but here we were able to use some new elevation datasets, refined breakwalls, measured river depths and include some previously excluded areas.

In addition the code has been substantially updated to be easier to use and modify, the model setup is more efficient, and the scripts are setup to use the new SapphireRapids nodes on Gadi.

An outline of the steps to run on NCI is provided below, followed by a summary of the key codes.
See comments within the code for further documentation.
You might also want to review [../GENERAL_GUIDANCE_ON_MODEL_SETUP.md](../GENERAL_GUIDANCE_ON_MODEL_SETUP.md).

# Outline of steps to run on NCI

``` bash

#
# COMPILE THE MODEL
#

# Get the right modules
module clear
source modules_SWALS_ifort_2023_B.sh
# NOTE: I have seen one model that segfaults with these compilers (NSW model, Sumatra with time-varying source).
#       The issue appears to be a compiler bug. (There are no problems on older ifort or gfortran, and the code seems valid).
#       A workaround for that specific case is to use SWALS_ifort_modules_2023_llvm.sh

# Compile
make -B -f make_model_ifort_sapphirerapids

#
# SANITY CHECKS OF AN INITIAL MODEL
#

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

# Post-process the simulation with R
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
# Look at elevation_all.vrt in GIS -- does it look right? Is the accurate, high
# resolution data where you expect it do be? Did you get the file order right?
# Are there unexpected artefacts? (e.g. caused by the interpolation used visually
# in GIS not matching the bilinear interpolation in SWALS).

#
# CONSIDER AN ALTERNAITVE global_dt
#

# Check the allowed timestep on each domain. This is in the logs, near where the model writes
# "stable" (if the timestepping seems stable) or UNSTABLE (but that was ruled out above). It
# can be found with:
grep -B1 stable *.log 
# To focus just on the global domains, which have domain index 1 (in this case):
grep -B1 stable *.log | grep "0000001 ,"
# If the minimum global domain timestep is much larger than the chosen global_dt, you might
# consider increasing global_dt. However that will fatten the halo regions, so there is
# a trade-off. You can test by comparing "load balanced models" with alternative global_dt.

#
# RUN A LOAD BALANCED MODEL
#

# If everything looks sensible then make a load balance file from within the
# multidomain directory
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
# Any problematic domains might be split or coarsened or removed, to reduce their run time and
# allow balancing the model.

# If the time in "comms1" is quite equal, but seems too large, consider running with fewer MPI
# processes. [Or, redesigning the model to reduce the halo fatness].


#
# Below here I assume the model design matches the intent, and it is reasonably efficient.
# Read the scripts to understand/modify for your case.
#

# Run some validation tests, and low-resolution variants to check convergence.
# In the case with a time-varying source, the source file was created with
Rscript make_initial_conditions_complex_historical_events.R
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

#
# After running the above jobs
# 
# You might want to make a folder to hold the *.sh.* files resulting from finished runs.
# * I usually call this dot_e_o_files
#
# Check the load balancing again. 
# * The efficiency of regular resolution models may be improved by making a new load_balance_file based on one of those runs. 
# * The efficiency of low resolution runs might also be improved (with a separate load-balance file). 
#
# Check the log files for mass conservation, energy decay, reasonable peak stage/speed etc.
# See ../analysis/check_log_files for more info.
#
# Make rasters in the usual way and check the results, e.g. on a node with 48 cores from inside the multidomain_dir
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
Rscript plot_gauges_perth_sumatra2004.R OUTPUTS/path_to_multidomain_dir_sumatra2004/gauge_RDS_file_name.RDS
# uniform y range of 1m
Rscript plot_gauges_perth_sumatra2004.R OUTPUTS/path_to_multidomain_dir_sumatra2004/gauge_RDS_file_name.RDS 1.0

# As above, Sumatra 2005
cd plots
Rscript plots/process_gauges_sumatra2005.R ../OUTPUTS/path_to_multidomain_dir_sumatra2005
cd ..
# varying y range
Rscript plot_gauges_perth_sumatra2005.R OUTPUTS/path_to_multidomain_dir_sumatra2005/gauge_RDS_file_name.RDS
# uniform y range of 0.3m
Rscript plot_gauges_perth_sumatra2005.R OUTPUTS/path_to_multidomain_dir_sumatra2005/gauge_RDS_file_name.RDS 0.3
#


#
# Run the random scenarios
#

# Set-up a script that can create all the PBS submission scripts, for all
# scenarios you want to run. 
vim create_random_ptha_qsub_scripts_sealevel60cm.R
# Then make the scripts.
Rscript create_random_ptha_qsub_scripts_sealevel60cm.R
# Eyeball a few scripts. 
# If they look OK, then submit one to see if it works.
qsub run/midwest_sealevel60cm_1001.pbs
# It's probably best to move the submitted scripts to another directory, to keep
# track of which have been submitted.
mkdir run/midwest_submitted_jobs
mv run/midwest_sealevel60cm_1001.pbs run/midwest_submitted_jobs/

# Once you're convinced that it's running correctly (including tarring the
# output folders) you can submit all jobs with e.g.
for i in run/midwest_sealevel60cm_10*.pbs; do echo $i; qsub $i; mv $i run/midwest_submitted_jobs/; done
# Be careful as this typically uses lots of compute.
# * Avoid overwhelming the queue, and ensure you have enough disk/inodes/KSU.
# * In practice you might want to submit in batches of 10 (or whatever) to deal
#   with those limits, and double-check everything is OK without using all your
#   compute quota.


#
# After the random scenarios have finished
#
 
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
# Look at the resulting pdf plot, and the summary statistics (example below)
# - The statistics should show that all runs finished, had good mass
#   conservation, had energy that did not increase significantly from the initial
#   value (although tiny increases in the first few steps are common and
#   acceptable for models with a static source, while if the source is
#   time-varying then there may be large energy increases over time as the source
#   is applied).
# - The plot should show the energy decaying over time.
#   - Tiny increases in the first few steps are acceptable
#   - Once flows interact with the boundary we could "in-principle" have energy
#     increases, although in practice don't expect that.
# - The kinetic energy should also show a general decay (albeit with some bumps etc) 
# - The max-stage and max-speed can vary, but beware crazy oscillations or very
#   extreme values.
#
## [gxd547@gadi-login-04 check_log_files]$ Rscript check_log_files.R 0.6 24 "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/*/*/multi*00001.log"
## [1] "Did the model runs finish?"
##    Mode    TRUE 
## logical     369 
## [1] "Mass conservation errors relative to initial volume"
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 1.644e-16 1.874e-16 1.977e-16 1.992e-16 2.093e-16 2.576e-16 
## [1] "Mass conservation errors relative to boundary flux integral"
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 7.220e-11 1.827e-10 3.658e-10 2.077e-08 1.827e-09 1.216e-06 
## [1] "(Maximum energy - initial energy) relative to (2x maximum kinetic energy), BEFORE BOUNDARY FLUXES"
## [1] "(typically very small unless the source is time-varying)"
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.0000000 0.0002719 0.0007013 0.0006621 0.0009920 0.0015067 
## [1] "Time index with largest kinetic energy (usually near start)"
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   2.000   2.000   3.000   3.477   4.000  11.000 
## [1] "Maximum energy increase between timesteps, relative to (2x maximum kinetic energy)"
## [1] "(typically very small unless the source is time-varying)"
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
## -0.0004571  0.0003264  0.0007534  0.0007239  0.0010201  0.0038118 


#
# At this point you might need to further investigate and fix any problematic runs.
#
# Before proceeding to probabilistic computations, the fixing also needs to: 
# * Ensure there is only one model simulation per scenario under the
#   ./OUTPUTS/name_of_this_batch_of_runs/ folder.
# * Ensure that all models have the same domain structure / resolutions / halo
#   size etc
#   * Or, edit the subsequent processing scripts so that they tolerate any
#     differences you introduce.
#     * This might not be easy, so try to avoid it.

#
# Once you think all the runs are OK, make the rasters
#
qsub post_process/create_tarred_rasters_from_tarred_multidomains.sh

#
# Check the max-stage exceedance-rate curves vs PTHA18
# * If these aren't good, consider doing more runs. Multiple Importance Sampling provides the theory for combining several runs.
cd ../analysis/max_stage_at_a_point/
# Edit extract_max_stage_at_a_point.R so all the input parameters are correct
# Maybe edit run_a_few.sh
qsub run_a_few.sh

```

# Code summary.

## Hydrodynamic model setup and compilation

[model.f90](model.f90) is the main driver program. This takes a number of input arguments (documented in the file). It uses:
* [model_multidomain_design_mod.f90](model_multidomain_design_mod.f90) to define variables that are most often modified.
  * See the code for explanation of the variables, and their default values.
  * Most default values can be overriden at runtime via input namelists, for example [multidomain_design_control_NNL4_defaultres.nml](multidomain_design_control_NNL4_defaultres.nml)
* [model_initial_conditions_mod.f90](model_initial_conditions_mod.f90) to define initial conditions.
* Other inputs are in [folders above this one](../), see README files therein for information.

The model requires loading some NCI modules.
* `source modules_SWALS_ifort_2023_B.sh`
  * I have seen compiler bugs that caused one SWALS model to segfault using these modules. 
    * The solution was to use another modules file (`modules_SWALS_ifort_modules_2023_B.sh`) that specified a different compiler version..
  * The problematic case was a NSW model, setup to use a time-varying Sumatra source model. 
    * The same model runs fine with gfortran and other versions of ifort. 
    * The same source model works fine with in some high-resolution WA models (using the same compilers, and different compilers).
    * I could not find any problem with the code (even using lots of debugging options in compilation). 
    * Hence why I blame a compiler bug.

Once the modules are loaded, the code can be compiled.
* `make -B -f make_model_ifort_sapphirerapids`

An initial model without any load balancing was executed with [run/test_model_sapphirerapids.sh](run/test_model_sapphirerapids.sh). 
* This used the input namelist [multidomain_design_control_NNL4_defaultres_defaultloadbalance.nml](multidomain_design_control_NNL4_defaultres_defaultloadbalance.nml).
* The result was used to create a load-balance file (to more evenly distribute the work, and thus improve parallel efficiency)
  * See discussion of `load_balance_script.R` below.

That load balance file was referenced in a new input namelist [multidomain_design_control_NNL4_defaultres.nml](multidomain_design_control_NNL4_defaultres.nml)
  * The load balanced model was then run in test mode using [run/test_model_with_load_balance_sapphirerapids.sh](run/test_model_with_load_balance_sapphirerapids.sh), to check it was working OK.
    * i.e. run times are acceptable, and not too much time is spent in `comms1`. 

## Simulating historic tsunamis

We test the model against tide-gauge observations for the 2004 and 2005 Sumatra tsunamis.
* [run/Sumatra2004_time_varying.sh](run/Sumatra2004_time_varying.sh) models the 2004 event with a time-varying source
  * The source was constructed with [make_initial_conditions_complex_historical_events.R](make_initial_conditions_complex_historical_events.R).
* [run/Sumatra2005.sh](run/Sumatra2005.sh) models the 2005 event.

Before running the probabilistic scenarios, we also ran convergence tests of
the above models using
[run/Sumatra2004_time_varying_lowres.sh](run/Sumatra2004_time_varying_lowres.sh)
and [run/Sumatra2005_lowres.sh](run/Sumatra2005_lowres.sh). 

We also ran an extreme scenario (derived by taking a big PTHA18 scenario and multiplying by 5) to ensure that the model
was reasonably stable in the inundation regime:
* See [run/extreme_source.sh](run/extreme_source.sh) and [run/extreme_source_lowres.sh](run/extreme_source_lowres.sh)

We also ran a very small scenario using
[run/small_source.sh](run/small_source.sh). From experience with earlier
models, the max-flux results from small scenarios can help to identify
any instabilities at nesting boundaries. Nowadays this is less common thanks
to the SWALS subroutine to smooth domain elevations at fine-to-coarse nesting boundaries.

## Plotting historic tsunamis

Code in [./plots](./plots) is used to process gauges for plotting. Then the actual plots are made with:
* [plots/plot_gauges_perth_sumatra2004.R](plots/plot_gauges_perth_sumatra2004.R) 
* [plots/plot_gauges_perth_sumatra2005.R](plots/plot_gauges_perth_sumatra2005.R)

In addition it is a good idea to inspect raster outputs; see [make_rasters.R](make_rasters.R) for one-off raster creation.

## Creation of qsub scripts to run the hydrodynamic model for random PTHA scenarios

[pre_process/create_random_ptha_qsub_scripts_sealevel60cm.R](pre_process/create_random_ptha_qsub_scripts_sealevel60cm.R) is used to make qsub scripts which run the hydrodynamic model for all the random scenarios (for the full-resolution runs), and also `tar` the resulting output folders (to prevent creation of too many files on NCI).
* The script can be run with
  * `Rscript pre_process/create_random_ptha_qsub_scripts_sealevel60cm.R`.
* It produces approximately 50 qsub scripts, that can be manually submitted later.

## Other miscellaneous code

* [post_process/create_plots_from_tarred_multidomain_dirs.R](post_process/create_plots_from_tarred_multidomain_dirs.R) can make basic png images of various flow maxima (stage, speed, flux) from the tarred multidomain directories.
* [post_process/create_tarred_rasters_from_tarred_multidomains.R](post_process/create_tarred_rasters_from_tarred_multidomains.R) makes rasters for an existing tarred multidomain folder, saving them to another tar archive. This approach of using tarfiles is important to avoid making too many files (violating our iinode quota on NCI).
* [post_process/make_rasters.R](post_process/make_rasters.R) is useful for making rasters from a single model run.
* [post_process/load_balance_script.R](post_process/load_balance_script.R) can be run from inside a multidomain directory, and will produce a file `load_balance_file.txt` which tries to distribute domains among MPI images to equidistribute the work from that run.
    * The files in [./load_balance_files](./load_balance_files) were produced this way (but apply to different model setups and core-counts).
* [post_process/make_domains_shapefile.R](post_process/make_domains_shapefile.R) can make a shapefile depicting the multidomain layout.
* [post_process/report_domain_runtimes.R](post_process/report_domain_runtimes.R) can summarise the time required by each domain, which is useful for optimization.

