# Radial dam break regression test

The code to setup the radial dam-break problem is in [radial_dam_break.f90](./radial_dam_break.f90).

Run the model and compare against the reference solution using:

    source run_model.sh

As this problem does not have an analytical solution, the reference solution was created using a model with grid-size about 8x finer than the default case. See comments and code in [make_highres_outputs.R](make_highres_outputs.R) for details on how this was done. The [test code](plot_results.R) checks that the coarse-grid solution is close to the fine grid solution, with variable-specific tolerances.

The figures below compare the coarse/fine grid solutions, and show excellent agreement.

![Figure 1: Comparison of flow arrival times on coarse and fine grids](radial_dam_break_reference_solution_arrival_time.png)

![Figure 2: Comparison of stage at the last time-step on coarse and fine grids](radial_dam_break_reference_solution_stage_last_time_step.png)

![Figure 3: Comparison of the x-directed flux (along y==0) at the last time-step on coarse and fine grids](radial_dam_break_reference_solution_uh_last_time_step.png)

![Figure 4: Comparison of the velocity at the last time-step (along y==0) on coarse and fine grids](radial_dam_break_reference_solution_vel_last_time_step.png)

![Figure 5: Comparison of the VH (along y==0) on coarse and fine grids; due to the symmetry of the problem this should be close to zero.](radial_dam_break_reference_solution_vh_last_time_step.png)

![Figure 6: Comparison of the flux maxima on coarse and fine grids](radial_dam_break_reference_solution_max_flux.png)

![Figure 7: Comparison of the speed-maxima (along y==0) on coarse and fine grids](radial_dam_break_reference_solution_max_speed.png)

![Figure 8: Comparison of the max-stage (along y==0) on coarse and fine grids](radial_dam_break_reference_solution_max_stage.png)

