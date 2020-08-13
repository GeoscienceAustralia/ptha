Figures and statistics for the paper
----------------------------------------------

Codes here do most of the analysis for the paper (and make most figures).

Before running any code in this directory we first had to run all the models and post-processing scripts in [../swals](../swals) as discussed therein. That created a set of post-processed outputs in [./gauge_RDS_files](./gauge_RDS_files), and we provide a [direct download of those files here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/OUTPUTS.zip). 

In addition you can find time-series plots for all model types and nearshore gauges [here](http://dapds00.nci.org.au/thredds/fileServer/fj6/PTHA/Nearshore_testing_2020/all_models_sites_vs_data.zip).

Here we explain the location of code used to make various figures in the paper. 
* Figure 1 - The figure with all the source models is made in ../sources/figures/
* Figure 2 - The figure with the model setup was made by powerpoint collation of several model snapshots. The model snapshots were made with ../swals/domain_plot.R
* Figure 3 - The figure with convergence test results [+ associated summary statistics] was made in ./convergence_checks/, using the script compare_gauges_yamakazi.R. Another script in that folder similarly treats another source. 
* Figure 4 - The figure comparing "fully nonlinear offshore" with "linear + manning offshore" was made in ./offshore_nonlinear_vs_linear_with_manning/compare_gauges.R. Note that after most models wre run I revised the nesting algorithm for the nonlinear-leapfrog solver -- the results for the paper were re-run with the updated algorithm (although it turned out to have no visible effect on the results for this problem).
* Figure 5 - The figure illustrating the "approximate linear-friction solution" was made in ./linear_friction_and_energy with the script plot_linear_frictionless_and_constant_friction.R. In addition some tsunami energy statistics [e.g. change in energy] are assessed and reported in ./linear_friction_and_energy/energy_tracking.R
* Figure 6,7,8 -- These figures compare modelled/observed gauge time-series at several sites, and are made in ./time_series_model_types/gauge_plots_sites.R. There is also a script in this folder to plot the model at all gauges.
* Figure 9, 10 -- These figures were made in ./tsunami_size/analysis.R, by the function "combined_plot_model_vs_data". The analysis.R script also reports many statistics that are used in this part of the paper (e.g. comparison of our statistics vs other studies).
* Figure 11 -- This was made in ./tsunami_size/analysis.R, by one of the calls to the function 'boxplot_relative_errors'
* Figure 12 -- This was made in ./tsunami_size/analysis.R, by the function model_vs_model_maxima_change()

