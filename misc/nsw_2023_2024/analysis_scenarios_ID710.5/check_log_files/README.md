The script [check_log_files.R](check_log_files.R) runs a number of basic QC checks on the tsunami inundation simulations by scanning the log file. See notes in the script header for details of the theory.

Run like
```
Rscript check_log_files.R ambient_MSL expected_end_time_h "quoted_string_globbing_the_log_files_of_interest"
```
e.g.
```
Rscript check_log_files.R 0.6 24 "../../swals/OUTPUTS/Sum*0.6/RUN*/multi*0001.log"
```
Here
* `ambient_MSL` --> The mean sea level assumed by the model 
  * This affects the potential energy calculations. If it's wrong you might appear to see poor energy conservation, see notes in the script.
* `expected_end_time_h` --> The time (hours) at which the model should finish.
  * Used to check if the model finished.
* `"quoted_string_globbing_the_log_files_of_interest"` --> A string that matches the log files to study 
  * Only 1 log file per model run is needed, because the script works with the Global properties only.

The checks confirm that:
* The runs finished
* Mass conservation errors are negligible. SWALS computes the mass error as:
  * `spatially_integrated_flow_depth - time_and_space_integrated_boundary_fluxes - spatially_integrated_flow_depth_at_time_0`
  * Theoretically it is zero; numerically, floating point imperfections in the integrals will lead to non-zero values, but we check it is negligible compared to the other terms.
* The energy behaves as expected. We check this because energy misbehaviour can be an indicator of model instability. 
  * The energy should be non-increasing before the tsunami reaches the boundary, UNLESS the source is applied over time. 
   * Tiny numerical increase may occur in the leap-frog solver (even with an instantaneous source) but
     * The script reports their relative size, and they should be very small. 
     * They shouldn't be visually evident in the plots
  * Once the tsunami reaches the boundary, the model's volume becomes time-varying. 
    * But we tend to expect that the "kinetic-energy + available-potential-energy" will be decreasing due to boundary radiation (noting that wetting/drying should affect a small part of the model domain). 
      * We can compute this quantity (up to some unknown constant) using the theory explained in the check_logs.R script comments. 
      * Note that in-principle boundary inflows could cause this quantity to increase.
  * We also compare the kinetic energy over time, which should show a general decrease due to friction and boundary radiation, although in-principle boundary inflows could cause it to increase. 
