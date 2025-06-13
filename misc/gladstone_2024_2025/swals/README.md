# Tsunami model solver for Gladstone
--------------------------------------------------

This folder contains code to run the tsunami models, including:
* Scenarios like [historical events](../sources/like_historic/)
* [Random PTHA18 scenarios](../sources/hazard/)
You might also want to review [../doc/general_guidance_on_model_setup.md](../doc/general_guideance_on_model_setup.md) for further guidance.

The key folders are:
- [test](test) -- test cases
- [test-full](test-full) -- full duration tests and validation cases like historic events
- [run_ptha](run_ptha) -- run all scenarios for PTHA
- [post_process](post_process) -- post-processing of model outputs
- [multidomain_design_control](multidomain_design_control) -- model configuration files
- [OUTPUTS](OUTPUTS) -- model outputs

The model is specified with the following fortran files:
- [model.f90](model.f90) -- main model code which makes use of the following two included files
- [model_initial_conditions_mod.f90](model_initial_conditions_mod.f90) -- initial conditions
- [model_multidomain_design_mod.f90](model_multidomain_design_mod.f90) -- domain design
as well as a namelist file to control commonly modified inputs (see [multidomain_design_control](multidomain_design_control)).

## Outline of steps to run on NCI

## COMPILE THE MODEL
I've updated the compiler to the latest version, i.e. Intel's llvm compiler.
``` bash
# Get the right modules
module clear
source ../modules_SWALS_ifx_2024.sh
# Compile with a shell script that sets parameters then calls make
make build
make debug
make clean  # Optional
```
This will make an executable called [model_build](model_build) and then [model_debug](model_debug). The latter will be used in stress testing below.  If you're not going to compile again, or don't mind waiting a few minutes next time you recompile, you can `make clean` to declutter the swals directory of object files. If you also want to delete the binaries use `make clean_all`.

You can use different compile options `make old_nesting` to have the code use an alternative nesting strategy. See the makefile. Each of these should make a different executable.

## SANITY CHECKS OF AN INITIAL MODEL

### 1. Run a short model
```bash
qsub test/model_debug.pbs
```

### 2. Inspect the log files
Go to the multidomain output folder (here the folder name is indicative) and open a log file, e.g., multidomain_log_image_00000000000000000001.log.
- [x] Did the model finish? (should end with printing of timers)
- [x] Are the printed namelist variables correct? (grep for &MULTIDOMAIN)
- [x] Are any domains unstable with the default timestepping, on any MPI process? 
- [x] Did any model produce NaN, on any MPI process?

``` bash
cd  OUTPUTS/model_debug-test-ambient_sea_level_2.42/RUN_20240708_154120216
grep UNSTABLE multi*.log
grep NaN multi*.log
```
The above two lines should produce no output. If they do produce output, check the global timestep - see changing the global timestep below.

### 3. Post-process the simulation with R
Run the following from inside the current folder.
```bash
module clear
source ../modules_R_431.sh
```

#### 3a. Create & check the domain shapefile (load R modules as needed)
Run the following from inside the multidomain output folder.
```bash
Rscript ../../../post_process/make_domains_shapefile.R "."
```
- [x] The above command created a directory "domains_shapefile" in the working folder.
scp the files to your computer and check in GIS -- are the domains where you expected them to be?

#### 3b. Make and inspect the model elevation
Run the following from inside the multidomain output folder.

Make rasters with the elevation (here assuming 24 cores -- which needs an interactive job on NCI. Too many cores might run out of memory). The below should make many raster (tifs) in the working directory with names matching elevation0*.tif. Combine into a vrt (to treat as a single file) with
```bash
Rscript ../../../post_process/make_rasters.R '.' 48 max_stage max_speed elevation0
gdalbuildvrt -resolution highest elevation0.vrt elevation0*.tif
gdalbuildvrt -resolution highest max_stage.vrt max_stage*.tif
gdalbuildvrt -resolution highest max_speed.vrt max_speed*.tif
```
- [x] Look at elevation_all.vrt in GIS -- does it look right? Is the accurate, high resolution data where you expect it do be? Did you get the file order right? Are there unexpected artefacts? (e.g. caused by the interpolation used visually in GIS not matching the bilinear interpolation in SWALS).
- [x] Look at the max_speed at all the polygons where the initial stage has been modified. If the initial stage polygon doesn't fully cover a low lying region then it can lead to an uneven initial sea level  around the polygon, akin to a dam break.

### 4. Consider an alternative global_dt
Run the following from inside the multidomain output folder.

Check the allowed timestep on each domain. This is in the logs, near where the model writes "stable" (if the timestepping seems stable) or UNSTABLE (but that was ruled out above). It can be found with focus just on the global domains, which have domain index 1 (in this case):
```bash
grep -B1 stable *.log | grep "0000001 ,"
```
If the minimum global domain timestep is much larger than the chosen global_dt, you might consider increasing global_dt. However that will fatten the halo regions, so there is a trade-off. You can test by comparing "load balanced models" with alternative global_dt.
- [x] Here the global_dt is `1.6s` and the smallest global domain timestep is `1.639s` so we are good to go.

## Validation with Historic events
Use source inversions to simulate [Tohoku 2011](validation/tohoku2011.pbs) and [Solomon Islands 2007](validation/solomon2007.pbs). The gauge records for the event are plotted in `validation/gauge_plots` using their respective `compute_gauges_*.R` and then `plot_gauges_*.R`.


## Performance optimisation
Using v5 of the model, the model was load balanced using the Solomon islands event. Cross checking against the Tohoku event (using [post_process/load_balance_script.R](post_process/load_balance_script.R)) leads to the following results, which suggest no further load balancing is required for events like these two.
```text
Range of partition total times: 29700-to-29900s 
                  (difference): 197.6s 
             (as a percentage): 0.6629%
Previous time range: 29520-to-30230s 
      (difference) : 710.4
Potential run-time reduction (difference in max times):327.7s
      Potential run-time reduction (% of old max time):1.084%
```

The three slowest domains are reported from [post_process/report_domain_runtimes.R](post_process/report_domain_runtimes.R) for the Solomon and Tohoku events are the 10m domains at Lady Elliot Island (time = 1.38 hours), One Tree Island (1.38 hours) and Heron Island (1.16 hours). This is expected, since they are in 10 m grids and have deep areas of up to 60 m deep.

The load imbalance due to parallel communication in the Solomon 2007 event is 4.42%, averaged over the 8 images, and is computed using:
```bash
grep " comms1:" *.log | awk 'NF>1{print $(NF-1)}' | awk '{sum += $1; count++} END {if (count > 0) print sum/count}'
```

## Stress Tests
1. Run the [extreme source](test-full/extreme_source.pbs) to test the model with a large source. In this case the initial sea surface deformations are 5 times larger than a magnitude 9.4 scenario from `kermadectonga2`. This extreme source was also used in the NSW study. Based on the results simulating the extreme source, the preliminary model multidomain was adjusted to remove global (leapfrog linear) to first level nesting (finite volume rk2) communication in the south of the model. Additionally, the preliminary model's elevation was smoothed at a nesting boundary at Kolan River. This unusual event spent 8.9% of time on comms1. Using grep, we ensured that there were no "UNSTABLE" or "NaN" found in the multidomain logs (see above).

2. Run the [low res version of the extreme source](test-full/extreme_source_lowres.pbs) to see if the results change with a coarser resolution and larger timestep, by a factor of 2. Version 5 gave small changes in results (acceptably similar).

3. Run the [small source](test-full/small_source.pbs) which can detect instabilities in the nesting boundaries. Make the last_timestep_VH and last_timestep_UH rasters at the last timestep to ensure nothing at the nesting boundaries is indicative of instability. 

## Debug any failed runs
Also `grep NaN *.log` can help. An early version of a model showed instability and was re-run in debug mode using 2 nodes with 16 images (so I could get the result faster). The run time came down from 9 hours to 6 hours and wouldn't be significantly improved by load balancing (15 mins). This check revealed the instability was triggered by a river crossing the first level nesting boundary. So I patched it at 100 m AHD. After fixing this issue the [post_process/check_log_files.R](post_process/check_log_files.R) check printed:
```
[1] "Did the model runs finish?"
   Mode    TRUE 
logical     313 
[1] "Mass conservation errors relative to initial volume"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
1.470e-16 1.766e-16 1.907e-16 2.657e-16 2.114e-16 2.022e-15 
[1] "Mass conservation errors relative to boundary flux integral"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
1.480e-10 8.650e-10 2.856e-09 7.049e-08 7.080e-09 5.024e-06 
[1] "(Maximum energy - initial energy) relative to (2x maximum kinetic energy), BEFORE BOUNDARY FLUXES"
[1] "(typically very small unless the source is time-varying)"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 7.479e-05 2.801e-04 3.233e-04 5.497e-04 9.263e-04 
[1] "Time index with largest kinetic energy (usually near start)"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   3.000   3.709   5.000  15.000 
[1] "Maximum energy increase between timesteps, relative to (2x maximum kinetic energy)"
[1] ", normalised for volume change but not wetting and drying."
[1] "(typically very small unless the source is time-varying)"
      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
-0.0002932  0.0001223  0.0003070  0.0003331  0.0005616  0.0016091 
```

### Once you think all the runs are OK, make the rasters
The script will need editing to your output folders
```bash
cd post_process
qsub create_tarred_rasters_from_tarred_multidomains.pbs
```

### Now go to the Analysis folder
[../analysis](../analysis/)


# Quick checklist:
Once the model compiles and runs, use the following checklist to aid model development.

1. **test/**
- model_load_balance.pbs
    - [x] Create & check the domain shapefile
    - [x] Look at elevation_all.vrt in GIS
    - [x] Check that max_speed is zero near places where the initial stage has been modified
    - [x] Are the printed namelist variables correct? (NNL4_defaultRes.nml) Find '  ' replace '' to make more readable.

2. **test-full/**
- small_source.pbs
    - [x] Inspect energy decay
    - [x] Inspect max_speed, last_timestep_UH and last_timestep_VH on nesting boundaries.
- extreme_source.pbs
    - [x] Inspect energy decay. Check for instabilities in log files  
    - [x] Inspect max_speed, last_timestep_UH and last_timestep_VH

3. **Validation/**
- tohoku2011.pbs
    - [x] Inspect energy decay
    - [x] Inspect gauge plots
- chile2010.pbs
    - [x] Inspect energy decay
    - [x] Inspect gauge plots
- Solomon2007_1_19.pbs  
    - [x] Inspect energy decay
    - [x] Inspect gauge plots
    - [x] Report domain runtimes
    - [x] Can improve load balance?

4. **Convergence of max_stage**
    - [x] validation/tohoku2011_lowRes.pbs
    - [x] validation/solomon_lowRes.pbs
    - [x] test-full/extreme_source_lowRes.pbs 

5. **All**
    - [x] Are any domains unstable with the default timestepping, on any MPI process? `grep UNSTABLE */*/*.log`
    - [x] Did any model produce NaN, on any MPI process? `grep NaN */*/*.log`
    - [x] How much KSU to expect.

Once complete, check the logs:
```
cd ../../OUTPUTS/ptha/sea_level_vary`
-[x] `grep NaN */*/*.log`
-[x] `grep UNSTABLE */*/*.log`
-[x] Use the `post_process/check_log_files.R`
-[x] Inspect the energy decay pdf it made
```
