# Tsunami modelling for the Perth to Geraldton area in 2023-2024

## Part of the WA Tsunami Inundation Modelling Project with DFES

This folder is for placing high-resolution grids throughout the region from perth to Geraldton.
The main coastline in this stretch is the midwest, although it includes other areas further south down to Guilderton.

An outline of how to setup the model is provided in [GENERAL_GUIDANCE_ON_MODEL_SETUP.md](GENERAL_GUIDANCE_ON_MODEL_SETUP.md).
A summary of software that might help is in [USEFUL_TOOLS.md](USEFUL_TOOLS.md).
They can be converted to pdf format with [make_pdf.sh](make_pdf.sh).

Additional useful content may be found in the [Greater Perth Modelling files](../Greater_Perth).
For instance, it includes checks of sea-levels and initial model developments that are also useful here.

Key folders are:

* [./analysis](./analysis) - Computations that use multiple simulated scenarios to produce outputs. Such as hazard products, and inundation maps that correspond to JATWC categories.
* [./breakwalls](./breakwalls) - Create 3d lines that are burned into the elevation model, to enforce breakwalls and other local high-points irrespective of model resolution.
* [./elevation](./elevation) - Elevation data + preference order for the SWALS model
* [./inverts](./inverts) - Create 3d lines that are burned into the elevation model, to enforce channels and other local low-points irrespective of model resolution.
* [./gauges](./gauges) - Locations of point-gauge output used in the SWALS model
* [./multidomain_design](./multidomain_design) - Create the boxes defining domains in the multidomain
* [./sources](./sources) - Tsunami source models
* [./swals](./swals) - Hydrodynamic model and outputs

# Workflow
Description of workflow followed for the PTHA, after model testing and validation.

## Testing and validation
Completed model test runs.
Model run on Sumatra2005 and Sumatra 2004 for comparison with tide gauges.
Compared well in high resolution zones.
Mostly reasonable agreement in medium resolution for initial timing and peak, but wave train and late stage less well represented.
Convergence test: Sumatra2005 and Sumatra 2004 in low res.
Extreme source test for model blow ups.

## Job submission
- 1 submitted are finished
- submitted the rest
- `grep NaN */*/multi*00001.log | wc`
  14529   38951 2611479

`for i in */*/multi*00001.log; do tail -n2 $i | head -n1 ; done`
Mostly around 9400 to 9500 seconds.

**Slowest timing**:
ptha18_random_scenarios_sunda2_row_0035771_Mw_77_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log:Total WALLCLOCK time:       10005.65672157

### One at 149.48093546 seconds. Too short?
`grep "Total WALLCLOCK time:         149.48093546" */*/multi*0001.log`
random_sunda2/ptha18_random_scenarios_sunda2_row_0110918_Mw_96_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log:Total WALLCLOCK time:         149.48093546

### Inspect energy curves:
- High max stage for long time. Floods a high basin.
- Max speed increases in when boundary condition comes in for many scenarios.
- All completed runs had believable decay of total energy. 

### 7 runs failed, had NaN:
`grep NaN */*/multi*00001.log | cut -d " " -f1 | uniq`
- random_sunda2/ptha18_random_scenarios_sunda2_row_0103658_Mw_91_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0108115_Mw_94_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0108515_Mw_94_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0109366_Mw_95_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0109559_Mw_95_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0110775_Mw_96_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log
- random_sunda2/ptha18_random_scenarios_sunda2_row_0110918_Mw_96_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log

### Re-ran the first one using the -DOLD_PROCESS_DATA_TO_SEND_B4FEB22 flag
It still blew up. Searching for the first NaN, it looks like this could be the offending tile:
i:            1  local_ti:           22  local_ni:           48  ti:
           4  ni:           16  lower_left_nx:         5251        2161
  upper_right_nx:         6300        3240  lower-left:    103.500000000000    
  -39.0000000000000       dx:   1.666666666666667E-002  1.666666666666667E-002
  lw:    17.5000000000000        18.0000000000000
### Try re-running with -TRACK_MULTIDOMAIN_STABILITY flag
Confirmed it's the above point.

### Re-ran with SWALS smoothing
Using `domain%smooth_elevation_near_point` in `model_initial_conditions_mod.f90` at this point.
Used script `midwest_redo1.pbs` to run the first one.
It worked.
Reran the remaining 6 that failed with `midwest_redo2.pbs`.
It also worked (checked for NaNs in all the log files after untarring).

### Use the analysis/check_log_files.R to inspect results
The following runs had a finish time more than 10 seconds earlier than 24 hours:
- [1]  "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0108115_Mw_94_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log"
    - inner at  x:    131.480555555556      ; y:   -12.2861111111111    
- [2] "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0108515_Mw_94_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log"
    - inner at  x:    131.489814814815      ; y:   -12.2953703703704     
- [3] "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0109366_Mw_95_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log"
    - inner at  x:    115.000308641975      ; y:   -33.6688271604938     
- [4] "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0109559_Mw_95_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log"
    - inner at  x:    131.475000000000      ; y:   -12.2842592592593     
- [5] "../../swals/OUTPUTS/ptha18-midwest-sealevel60cm/random_sunda2/ptha18_random_scenarios_sunda2_row_0110775_Mw_96_HS-full-ambient_sea_level_0.6/multidomain_log_image_00000000000000000001.log"
    - inner at  x:    131.484259259259      ; y:   -12.2842592592593     

### Fixed the above by adding a smoothing region near (131, -12) and expanding the one near (115, -33).
Rerun using `midwest_redo3.pbs`.
Only one failed `sunda2_row_0109366_Mw_95_HS`.
It had an innerB fail at

In multidomain_log_image_00000000000000000002.log
```
innerB
  (NB: overflow just after mpi-comms can result from error-stop on other images)
 error-count:           93
 md%domains(          22 )%U(          13          17           2 )
 -9.095186590220378E+041
 x:    114.997067901235      ; y:   -33.6683641975309     
 time:    15307.1600000019     
 stg:elv,   5.394341630932504E+032 -9.095186590220378E+041
 -3.131641885406360E+044   3.68805329633864 
```
### Added aggressive smoothing near the problematic region
It's a coarse-to-fine boundary in the second level of nesting.
It finished.

### Ran energy analysis
 "Did the model runs finish?"
 ```
   Mode    TRUE 
logical     369 
[1] "Mass conservation errors relative to initial volume"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
4.490e-17 5.640e-17 6.000e-17 3.054e-15 6.560e-17 1.104e-12 
[1] "Mass conservation errors relative to boundary flux integral"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
2.400e-11 5.550e-11 1.179e-10 7.781e-09 5.612e-10 4.997e-07 
[1] "(Maximum energy - initial energy) relative to (2x maximum kinetic energy), BEFORE BOUNDARY FLUXES"
[1] "(typically very small unless the source is time-varying)"
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000000 0.0002491 0.0007014 0.0006604 0.0009999 0.0015052 
[1] "Time index with largest kinetic energy (usually near start)"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   3.000   3.507   4.000  11.000 
[1] "Maximum energy increase between timesteps, relative to (2x maximum kinetic energy)"
[1] "(typically very small unless the source is time-varying)"
      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
-0.0004518  0.0003264  0.0007587  0.0007250  0.0010227  0.0039387
```
### Post-processing
Create rasters using R script.
```
Error: All the multidomain tar files are completed
Execution halted
```
## Analysis

### Max stage at a point
Ran the max_stage comparison to PTHA18 at several points.
One offshore point at 114.73777, -31.305 looks higher for out-rise but not sunda2.
Otherwise, all looks in general agreement.

### Probabilistic inundation 
Ran the probabilistic inundation exceedance calculations for max depth and max stage.

#### Comments on results
General agreement for rate of inundation between greater Perth revised 2023 run and the midwest run.
Discrepancies:
- It looks higher resolution
- Greenough river. Midwest run does not predict inundation at the caravan park
- Some inland low lying land artefacts are no longer there around Jurain Bay and Leeman.
- Low lying areas near Sandy Cape (south of Green Head) are filled
- Generally further inundation into patched rivers
- Geraldton sheds not submerged
- 1/100 zone (and marine warning) more extensive on land at:
  - Geraldton fisherman's wharf
  - Geraldton Francis street Jetty and train tracks nearby
- Abrohols Islands have significantly better resolution:
 - East Wallabi Island now has no inundation on Flag hill by the airport
- Shell beach (south of Leeman) has 1/100 better defined inundation onto Indian Ocean drive


### JATWC zones
  Note that the folders `ATWS_ZONES`, `elevation_contours` and `Inundation_zones` were missing and needed to be copied in.
  Computed scenario statistics for `Perth Coast`, `Lancelin Coast` and `Geraldton Coast` zones.
  Steps 2 and 3 only done for `Perth Coast`: map threat and convert to polygon file.
