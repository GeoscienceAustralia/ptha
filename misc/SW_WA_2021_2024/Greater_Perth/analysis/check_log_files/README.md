The script [check_log_files.R](check_log_files.R) runs a number of basic QC checks on the tsunami inundation simulations by scanning the log file. The checks confirm that:
* The runs finished
* Mass conservation errors are negligible. 
  * SWALS computes the mass error as the `spatially_integrated_flow_depth - time_and_space_integrated_boundary_fluxes - spatially_integrated_flow_depth_at_time_0`. 
  * Theoretically it is zero; numerically, floating point imperfections in the integrals will lead to non-zero values, but we check it is negligible compared to the other terms.
* The energy behaves as expected. We check this because energy misbehaviour can be an indicator of model instability. 
  * In all cases the energy should be non-increasing until the tsunami reaches the boundary, as expected in a closed domain with friction. 
    * Tiny numerical increase are not uncommon in the leap-frog solver, but the script checks for this and we find any increases are tiny. 
  * Once the tsunami reaches the boundary, the model's volume becomes time-varying. 
  * But we tend to expect that the "kinetic-energy + available-potential-energy" will be decreasing due to boundary radiation (noting that wetting/drying should affect a small part of the model domain). 
    * We can compute this quantity (up to some unknown constant) using the theory explained in the check_logs.R script comments. Note that in-principle boundary inflows could cause this quantity to increase.
  * We also compare the kinetic energy over time, which should show a general decrease due to friction and boundary radiation, although in-principle boundary inflows could cause it to increase. 

