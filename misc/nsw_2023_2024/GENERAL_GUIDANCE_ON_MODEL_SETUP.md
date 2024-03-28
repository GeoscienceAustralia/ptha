# Ideas for model setup 

The suggestions below are arranged according to subdirectories where the inputs are created. The order is roughly sequential, although when troubleshooting the model it is often necessary to revise earlier steps.

Each subdirectory should contain a README.md file that gives a more detailed "how-to" explanation.

## `./elevation`: 

Collate elevation data for the model in raster format and make a file that lists the rasters in order of preference.
* The model works with a set of one or more rasters, overlapped with a preference order. 
* They should have the same projection but usually have different resolutions to better resolve important areas. 
* In combination the rasters should completely cover (and ideally exceed) the model domain.
* When the rasters are overlapped in order of preference, they should give a good representation of the elevation for modelling. 
  * Check this in GIS.
* The model elevation is set from these rasters with bilinear interpolation, using lower preference rasters at sites where higher preference rasters were missing data.

## `./point_gauges`: 

Make a point gauges csv file
* This defines locations where we store flow time-series
* Each gauge is assigned an ID. It should be unique when stored as a single precision real.

## `./multidomain_design`: 

Define the model extent and create polygons showing where different nesting levels apply. 
* For each nesting level, you make polygonal regions to be covered by rectangular domains (with a chosen size).
* Some domains can also be flagged for coarsening. This is useful if a few domains have excessively fine resolution and are slowing the model (typical for domains touching deep areas). 
* For efficiency, domains should have several hundred cells east-west and north-south, e.g. 300x400 cells. 
  * Too few cells in the y-direction will not run efficiently if we use many OPENMP threads.
    * In practice, on Sapphire-Rapids nodes we should use `OMP_NUM_THREADS=13`
    * The Cascade Lake nodes are more flexible (e.g. `OMP_NUM_THREADS=2 - 12`). 
  * Too many cells is probably OK, although there may be some cache benefit to not having too many cells ??
    * Hypothesis: If there's an effect, it might matter most in the x-direction.
  * It's fine to have some domains with few cells, so long as most computational work isn't spent there. 

## `./sources`: 

Optionally put files that define tsunami sources here

## `./breakwalls`: 

Optionally put files that define breakwalls here.
* These are 3d lines that will be burned into the model elevation. They are defined in CSV format.
* They are useful to ensure that flow barriers are represented in the model's elevation.
* We burn the maxima of the breakwall elevation and the raster-derived elevation into the model.

## `./inverts`:

Just like `./breakwalls` but for elevation minima, rather than elevation maxima. 
* Inverts are typically used to ensure minor flow paths are represented in the model.


## `./swals`: 

Contains code to run the tsunami model for hazard scenarios and testing, and produce raster outputs and basic plots.

### Setup the hydrodynamic model
* Start by modifying the `model_multidomain_design_mod.f90` to mostly match what you want. 
  * Variables can be set at run time via two namelists. These are defined in an input file that is passed to the main program.
  * For more control you may modify `model_initial_conditions_mod.f90` (initial conditions) or `model.f90` (everything else).
* Get a model to compile correctly with your design and input data.
* Try to run it in `test` mode (for a short but non-negligible time). 
  * If it fails, figure out why and fix it.

### Once it runs:
* Check that the domains are in the right place.
  * The script `make_domains_shapefile.R` in the `swals` folder can make a single shapefile, showing all the domains, with attributes having information on output folder names etc. 
  * Run with `Rscript make_domains_shapefile.R path_to_the_multidomain_output_folder`. It will make a folder `domains_shapefile` in the multidomain output folder.
* Check the elevation produced by the model, and update the model inputs as required.
  * The script `make_rasters.R` in the `swals` folder can be used to make rasters (for each domain) in the model output folder. The rasters can be treated collectively by making a VRT file (e.g. `gdalbuildvrt -resolution highest combined_elevation.vrt elevation0_domain*.tif`).
* Grep the log files for `NaN` values, or `UNSTABLE` (which indicates the timestep is too large). Edit the model inputs and/or the multidomain design as needed to ensure there are no unstable domains.
* Check the namelist print outs in a log file (find them by searching for `&MULTIDOMAIN` in a log file). Do the variables match what you think you input?  
* Check the other model outputs
  * e.g., Is the correct timestepping method used? Are gauges in the right place? Does manning friction match expectations?
* Consider adjusting the domain placement and/or coarsening some domains. 
  * Eyeball the model and think what you need.
  * Have a look at the time required by each domain using `report_domain_runtimes.R`. 
    * Domains taking much more time than others are candidates for coarsening/splitting/removing, as they might interfere running the model efficiently (load balancing, see below). 

### Compare the `stationary_timestep` information in the log files with the timestep you have specified
And consider an alternative `global_dt`.
* Most often the `global_dt` is limited by the shortest timestep on the coarsest domains. 
  * The nested domains can usually do local timestepping (finite volume schemes)
* If lots of time is being spent in the global domain, increasing `global_dt` will speed it up.
* If the nested domains have very fat halos, it may be more efficient to reduce `global_dt` (which will reduce the halo thickness).

### At this point, you should have addressed any major problems with the input data or multidomain design. 

Now try to run a `test` model for a non-trivial number of timesteps (e.g. `final_time = 600.0_dp` or similar).
* If your original model meets this criteria, then just use that.
* Check that the initial condition is appearing correctly.
* Check that the model is still stable.

### Ensure it is running efficiently 
by checking and improving the parallel load balancing.
* For the model that you ran for a non-trivial number of timesteps, grep the log files for `comms1:`. 
  * This will show the fraction of time being spent on parallel communication. Ideally it will be a few percent on all processes.
  * If many images are spending a high fraction of time here (and one or more are not) you probably have a load balance problem, i.e., unequal work on the MPI processes. 
    * Images with less work end up wasting time in `MPI_Waitall` while they await data from other processes. 
* To improve the load balancing, make a new `load_balance_file` by running the `load_balance_script.R` from inside a multidomain output directory 
  * Go to the model output folder (e.g. `./swals/OUTPUTS/model_name_as_chosen/RUN_1234556677/`)
  * Start `R` and run `source("../../../load_balance_script.R")`.
    * Alternatively use `Rscript "../../../load_balance_script.R"`
  * This will make a new `load_balance_file.txt` in the working directory that in principle will lead to more equal work.
    * It also prints statistics on the expected improvement (they are usually optimistic).
  * Copy the `load_balance_file.txt` to a new location, update the input namelist (or `multidomain_design_mod.f90`) to point to this file, and try running the short model again. If it works, you should find the time in `comms1:` is much smaller. If not, try again with the model you just ran - sometimes it takes a few tries.
* Beware it may be impossible to achieve great load balance on models for which:
  * The timestep required on each domain changes substantially during the simulation.
    * In practice I haven't found big problems. 
  * A single domain takes lots of time (e.g. more than the average model time). 
    * To check if this is happening use `Rscript ../../../report_domain_runtimes.R` from inside the multidomain directory .
    * If a single domain is taking too much time, consider coarsening it, or removing it, or splitting it into smaller pieces.
      * Coarsening or removal can be implemented in `../multidomain_design`
      * Splitting can be implemented via the load_balance_file (by using more images on the offending domain).

In applications I am usually happy with `comms1` being around 10% of the `evolve:` time. 
* It can be much smaller for models with large domains and narrow halos.
* If the time in communication (`comms1`) is fairly equal among processes but remains a high fraction of the `evolve:` time, consider whether it is more efficient to:
  * Use fewer MPI processes
  * Redesign the multidomain structure so there are fewer, larger domains.
    * Beware this reduces opportunities for local timestepping, so might slow down the model.
  * Wait and run again. 
    * More than once I've found Gadi running slowly, with lots of time in `comms1`, yet the same model ran fast later without changes.
  * Other ideas (not successful to date)
    * Use an alternative MPI routine for the main communication
      * Try compiling with `-DCOARRAY_MPI_USE_ALLTOALLV`
      * Historically, this was similar or less efficient than the default approach.
    * Try coarray communication
      * Remove all preprocessor variables starting with `-DCOARRAY_`, but keep `-DCOARRAY`.
      * Historically this was less efficient than the default approach, especially with Intel compilers, because the coarray implementations were less mature than MPI.

### Repeat the above steps 

until you have a model that is stable, reasonably efficient, and seems to have reasonable resolution in the important areas.

### Run the model with some kind of validation test 
(e.g. a historical tsunami source, or a known reference result).
* Check if it behaves as expected. If not, do you need to improve the elevation/resolution/... ? 
* Check if the load balancing is still OK. If it isn't, try making a new `load_balance_file` using this model (but run some short tests to check it works ).
* Look at the results near nesting boundaries. Are there any unexpected flow features that might be nesting artefacts? 
  * This is more likely where fine-to-coarse transitions occur over steep topography. 
  * If it's a problem, consider moving the nesting-level regions.
* `analysis/check_log_files` may help detect more obvious instabilities.
  * Check the logs show good mass conservation (`Volume statistics (m^3) integrated over all domains and images:`)
  * Check the logs show sensible energy decay over time (`Global energy-total / rho (over all domains and images):`)

### Run the model with an extreme scenario (relative to the ones you intend to simulate). 
* For instance, take the initial condition of a big scenario and multiply it by 5.
    * Not too big, but a bit bigger than the worst scenarios.
* `analysis/check_log_files` may help detect more obvious instabilities.
    * Does it remain stable and have reasonable mass conservation / energy decay? Do the results appear sensible?
* If not, we might expect problems when running many alternative scenarios. Adjust the setup to make the extreme scenario stable.

### Run the model for a very small scenario. Does it show artefacts? 
* You could take a small earthquake, or scale down another initial condition.
* It can be easier to notice nesting artefacts in simulations with small flows.
* If found, they might be fixed by adjusting the nesting locations, or locally smoothing the elevation, or similar.
* `analysis/check_log_files` may help detect more obvious instabilities.
  * I have seen very small scenarios where energy increased over time once flow was interacting with the model boundaries. While we expect energy conservation in closed domains (and energy decay with friction), if the domain is open then the energy may increase or decrease with flows in and out of the domain.

### Check the convergence of the model
Re-run one or more of the above tests above using a finer grid, or a coarser grid if you think that is sufficient
  * i.e. halving or doubling `global_dx_arcmin` and `global_dt`. 
* Are the results the same in all important respects? 
  * If yes, the default grid is probably fine enough. 
  * If no, you need to improve the resolution in areas that show important differences, then repeat the steps above
    * If you use scripts to implement your workflows, this shouldn't be too difficult.

### Repeat the previous steps until you have a well-behaved model setup.
It will reduce the chance of difficult problems later on.

### At this point you can run all the hazard scenarios

Don't skip the earlier steps!
* They are designed to detect problems that could be much more troublesome when you run all scenarios.
* It is much easier to detect and fix these issues early on.

Before running all scenarios at full resolution, consider doing runs on a 2x coarser model.
* See discussion of re-running on coarser grid below.
* This offers another opportunity to notice problems, without blowing all your compute quota.

For models with many runs I usually write a script to create multiple PBS files. Each does a "batch" of runs.
* Consider running just one batch to start - check that it's working as planned.
* Later you can submit many at once, but ensure you don't overwhelm the NCI queue and have enough CPU/disk/inodes to finish the jobs.
  * Info on the latter is in `lquota -v` and `nci_account -P w85`. 
  * Compare with the model resources.
    * SWALS currently creates many output files, so you probably need to `tar` the multidomain folder immediately after the run completes, and before starting the next run.
    * In the longer term we could consider [reducing the number of output files, or combining them after running](https://github.com/GeoscienceAustralia/ptha/issues/10). Main challenges are backwards compatibility, and ensuring it doesn't kill the performance.

### Once all the hazard scenarios are run
* It's a good idea to scan the multidomain log files to ensure
  * The runs finished.
  * There are not NaN values in the results
* Screen the log files (`analysis/check_log_files`) to detect more obvious instabilities.
* When simulating a few hundred scenarios, it isn't unusual for one or two to go unstable, even using the approach above. 
  * I have seen cases where the NCI node failed
    * These instances often led to the creation of various MPI debug files (not from SWALS).
    * They could be fixed by simply re-running the same models.
    * I don't think these cases reflect problems with the code.
      * From experience in many cases, SWALS is bitwise reproducible when 2 runs use the same compiler options and prescribed multidomain partition (irrespective of the number of MPI/openmp processes).
  * Other problems have manifest as a sudden model blow-up (NaN) near a coarse-2-fine nesting boundary, without any precursor that something was wrong. 
    * They were practically solvable by changing minor details of the nesting interpolation (below)

#### Fixing unstable models
If the blow-up is sudden, without prior problematic behaviour, then it might help to change the nesting interpolation.

* To change the nesting interpolation, recompile the model with `-DOLD_PROCESS_DATA_TO_SEND_B4FEB22`. This uses a slightly different nesting algorithm.
  * Most models will not show differences. 
    * An exception is near chaotic processes (e.g. long-term eddy transport) where any minor change to the model will create differences.
  * The change in nesting is often enough to work around the "sudden" blow-ups I have encountered (which do not show precursors of instability).
    * Presumably they reflect some corner case in the nesting algorithms (be good to investigate this further).
      * The interpolation involves stage/elev/UH/VH and does not necessarily bound the velocity. That might cause blowups?
* Then re-run the unstable scenarios, which hopefully will now be stable and well-behaved.
  * Then fix the output folder structure. 
* Another approach to debugging models is to take the problematic case and
  * Make outputs frequent enough for visual inspection
  * Recompile with `-DTRACK_MULTIDOMAIN_STABILITY`, which should print information on the location of regions judged to be going unstable. 
  * The above information is often enough to identify the source of the instability, which can be fixed by editing the domain structure. 
    * Beware that complications will arise if you want to combine results from models with different multidomain setups / grid sizes etc. These might be worked-around, but it's better to avoid this by identifying problems with the models before doing large runs.

### Additional quality control measures

#### Consider re-running all the scenarios on a coarser grid
* Doing this before the full resolution runs could help to identify non-obvious problems with the model setup, before spending all your computational quota.
* Running on a 2x coarser model will take approximately `1/8` of the computational time, `1/4` of the disk space, and the same number of inodes.
* This will allow the convergence of the full analysis to be checked. 
  * Major differences could be a sign of model instabilities, or a poorly resolved model.
  * I did this in the original Greater Perth model, and it worked well. 
    * The noticeable differences were rare, and occurred in areas where the coarser grid wasn't fine enough to resolve important aspects of the topography (as expected). 

#### Consider re-running the analysis with another random sample of scenarios
* Assuming you've got the computational resources.
* I did this for the Tongatapu study, and it worked well.

## `./analysis`:

Contains code to check model results and do hazard calculations.

### Check the log files in `analysis/check_log_files`
* Mass conservation, energy decay, time-series of global maximum stage, speed, potential energy, kinetic energy.
* This can highlight problematic model runs.

### Compare the high-resolution hazard model results with the offshore PTHA in  `analysis/max_stage_at_a_point`
* They should agree quite well in deep water far from the coast
  * i.e. Sites at which differences between the offshore PTHA solver and the high resolution model are not important.
* This can help to show that you've simulated enough scenarios at high resolution, and haven't made other calculation blunders.
* Some differences are expected
  * The Monte Carlo method implies some variability in the high-resolution solution.
  * There are fundamental differences in the hydrodynamics
    * The offshore PTHA solver uses the linear shallow water equations, without any friction, on a 1 arcmin grid.
    * The nonlinear model is expected to deviate from this at sites where:
      * Nonlinearity matters (typically shallow sites, or sites that are affected by waves from shallow sites). 
      * The improved elevation data or resolution have an important effect.
      * At late times (at any site) due to the cumulative effects of friction.

### Run probabilistic inundation calculations in `analysis/probabilistic_inundation`
* Inundation probabilities for the logic-tree mean model, and various percentiles.
* The results are required to create the JATWC-style zones

### Compute preliminary JATWC zones in `analysis/jatwc_to_inundation`
* The results will later be manually modified to make actionable evacuation maps
