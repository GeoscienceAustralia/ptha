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
