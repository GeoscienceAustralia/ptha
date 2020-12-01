Probabilistic inundation calculations
-------------------------------------

This folder contains various scripts to do probabilistic inundation computations.

## Key files

* [make_probabilistic_inundation.sh](make_probabilistic_inundation.sh) creates rasters with max-depth exceedance-rates for various threshold depths. The rasters are created for domains 3,4,5,6,7 around Tongatapu, with a separate raster being created for the unsegmented kermadectonga2 model and each of the three segmented models (as this is how the kermadectonga2 source-zone is represented in PTHA18).
* [probabilistic_inundation.R](probabilistic_inundation.R) - This computes exceedance-rate rasters for a range of inundation-depth-thresholds. It needs to be given a directory name [containing many SWALS model runs] and an integer defining the domain of interest. Only one domain is done at a time. We use the shell script to run multiple domains and multiple simulation-sets [e.g. using different MSL, or different resolution]. Other options are controlled via hard-coded arguments inside the script. 
* [raster_plots.R](raster_plots.R) creates exceedance-rate rasters for the combined `segmented + unsegmented` source representation (this is important!), and makes a basic plot of the depths with a 2\% chance of exceedance in 50 years.
* [depth_vs_exrate_at_gauge.R](depth_vs_exrate_at_gauge.R) is used to compute the exceedance-rates of maximum depth (and maximum-stage) at two sites; the Parliament site (for application), and a PTHA18 point (to compare results with the original PTHA18 study). This is useful both for testing, and for the outputs themselves, see further discussion in the script comments. It is run with the following commandline calls: `Rscript depth_vs_exrate_at_gauge.R parliament` and `Rscript depth_vs_exrate_at_gauge.R 'ptha18_point_3458.3' `. 
* [plot_depth_vs_exrate_at_parliament.R](plot_depth_vs_exrate_at_parliament.R) can be run after the previous script has been run to generate RDS files with outputs at Parliament. It makes a plot of the results.
* [plot_stage_vs_exrate_at_gauge_3458.R](plot_stage_vs_exrate_at_gauge_3458.R) makes a plot comparing the max-stage exceedance-rates from this study and the original PTHA18 at a PTHA18 output point. Some differences are expected due to the use of different hydrodynamic models (which are run for different durations on different datasets with different resolutions and solver options), as well as the use of a random subset of scenarios here rather than all the PTHA18 scenarios. Neverthess we expect reasonable agreement if all is working well, and do obtain that in practice.

## Example of running the codes

Assuming all the hydrodynamic model results have been run, and all the dependencies are installed, one can do:

```
# Create the exceedance-rate rasters, and depth-vs-exceedance-rate curve information at a couple of key points.
source make_probabilistic_inundation.R

# Create the combined segmented/unsegmented exceedance-rate rasters and make a basic plot
Rscript raster_plots.R

# Make depth-vs-exceedance-rate plots at parliament
Rscript plot_depth_vs_exrate_at_parliament.R

# Make stage-vs-exceedance-rate plots at parliament
Rscript plot_stage_vs_exrate_at_gauge_3458.R

```
