Figures and statistics for the paper
----------------------------------------------

Codes here do most of the analysis for the paper (and make most figures).

Before running any code in this directory we first had to run all the models and post-processing scripts in [../swals](../swals) as discussed therein. That created a set of post-processed outputs in [./gauge_RDS_files](./gauge_RDS_files), and we provide a [direct download of those files here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/OUTPUTS.zip). 

In addition you can find time-series plots for all model types and nearshore gauges [here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/all_models_sites_vs_data.zip).

Here we explain the location of code used to make various figures in the paper. 
* Figure 1 - The figure with the five tsunami models was made in powerpoint by combining figures created with [../swals/plot_max_stage_and_elevation.R](../swals/plot_max_stage_and_elevation.R). These figures are also included in the data to download in [./gauge_RDS_files/](./gauge_RDS_files).
* Figure 2 - The figure with all the source models was made in [../sources/figures/](../sources/figures) after all the source models had been created.
* Figure 3 - The figure with the model setup was made by powerpoint collation of several model snapshots. The model snapshots were made with ../swals/domain_plot.R
* Figure 4 - The figure with convergence test results, and associated summary statistics, was made in [./convergence_checks/](./convergence_checks), using the script compare_gauges_yamakazi.R. Another script in that folder similarly treats another source. 
* Figure 5 - The figure comparing "fully nonlinear offshore" with "linear + manning offshore" was made in [./offshore_nonlinear_vs_linear_with_manning/compare_gauges.R](./offshore_nonlinear_vs_linear_with_manning/compare_gauges.R). Note that after most models were run I revised the nesting algorithm for the nonlinear-leapfrog solver -- the results for the paper were re-run with the updated algorithm (although it turned out to have no visible effect on the results for this particular problem).
* Figure 6 - The figure illustrating the "approximate linear-friction solution" was made in [./linear_friction_and_energy](./linear_friction_and_energy) with the script plot_linear_frictionless_and_constant_friction.R. In addition some tsunami energy statistics [e.g. change in energy] are assessed and reported in the script energy_tracking.R
* Figure 7,8,9 -- These figures compare modelled/observed gauge time-series at several sites, and are made in [./time_series_model_types/](./time_series_model_types) with the script gauge_plots_sites.R. There is another script in that folder which can plot the model and data at all gauges.
* Figure 10, 11 -- These figures were made in [./tsunami_size/](./tsunami_size) with analysis.R, specifically the function "combined_plot_model_vs_data". The analysis.R script also reports many statistics that are used in this part of the paper (e.g. comparison of our statistics vs other studies).
* Figure 12 -- This was made in [./tsunami_size/](./tsunami_size) with analysis.R, by a call to the function 'boxplot_relative_errors'
* Figure 13 -- This was made in [./tsunami_size/](./tsunami_size) with analysis.R, by a call to the function model_vs_model_maxima_change()

