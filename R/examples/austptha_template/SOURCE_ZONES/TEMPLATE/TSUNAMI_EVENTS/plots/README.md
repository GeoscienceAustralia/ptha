# Code to compare observed events with model scenarios on a single source-zone.

These notes provide some details on the code used to compare model scenarios and observations on a single source-zone. 

# Usage

Codes in this folder can be run after executing codes in the parent directory (as described [here](../README.md)). They can only be run on source-zones for which DART test data exists. 

Before running these scripts, the source-zone specific script corresponding to [../check_dart_example.R](../check_dart_example.R) should have been run in the parent directory. For instance, on the `southamerica` source-zone that script is called `check_dart_southamerica.R`. See [here](../../../dart_check_codes) for source-zone-specific DART-buoy scripts and further explanation. 

Once the `../check_dart_SOURCE_ZONE_NAME_HERE.R` code has been run, the script [gauge_summary_statistics.R](./gauge_summary_statistics.R) is used to process the resulting set of PTHA18 scenarios that have similar earthquake location and magnitude as the observations. Please note the script does not account for any earthquake-rate information, and it does not attempt to exclude scenarios that are impossible according to the PTHA18 (e.g. based on peak-slip limits). 

The script can be run from the commandline like:

    Rscript gauge_summary_statistics.R

and may optionally be followed by a plotting script:

    Rscript event_plot.R 5 7.5

*External users would have to partially recreate our PTHA18 folder structure to successfully run the gauge_summary_statistics.R script. Because this is non-trivial, links to key output files are provided below*.

## Obtaining the files produced by gauge_summary_statistics.R, without running the code.

The files produced by running [gauge_summary_statistics.R](gauge_summary_statistics.R) are often used in subsequent analysis because they contain observed and modelled tsunami time-series, and useful statistics for comparing those. 

For example, the corresponding file paths from our original analysis are stored [here](../../../../EVENT_RATES/config_DART_test_files.R) and are used to examine the statistical properties of random tsunamis [here](../../../../EVENT_RATES/stage_range_summary.R) and [here](../../../../EVENT_RATES/event_properties_and_GOF.R) and [here](../../../../EVENT_RATES/event_dart_coverage_vs_distance.R) ). 

For simplicity, these output files can be directly downloaded from the NCI THREDDS server at the following locations: [kermadectonga2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/catalog.html),
[kurilsjapan](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/catalog.html), 
[newhebrides2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/catalog.html), 
[puysegur2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/plots/catalog.html), 
[solomon2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/catalog.html), 
[southamerica](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/catalog.html), and
[sunda2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/catalog.html).


## Details of gauge_summary_statistics.R

The [gauge_summary_statistics.R](gauge_summary_statistics.R) script reads all the tsunami scenarios having "similar earthquake location and magnitude" as the observed event, which were actually selected via the `../check_dart_SOURCE_ZONE_NAME_HERE.R` script [discussed here](../../../dart_check-codes). It then extracts the time-series in a convenient form, performs some analyses, and makes some plots. The script saves its own workspace, separately for each tsunami event and each rigidity model. This permits access to all the variables defined by [gauge_summary_statistics.R](./gauge_summary_statistics.R) for each case, by loading the relevant file. PTHA18 analyses repeatedly make use of this.

Beware [gauge_summary_statistics.R](gauge_summary_statistics.R) does not give any consideration of the earthquake rates (which are computed later in [../../../../EVENT_RATES](../../../../EVENT_RATES)). Also, it does not exclude scenarios that are "impossible" according to the peak-slip limits in PTHA18. Such scenarios are excluded in later processing (example - the script [stage_range_summary.R](../../../../EVENT_RATES/stage_range_summary.R) does this around lines 101-116 -- and you will see similar exclusions in other scripts). 

### Structure of output files 

Next we give a more concrete example of the output file structure using the `kurilsjapan` source-zone as an example. On the `kurilsjapan` source-zone we had two test events (both defined in [check_dart_kurilsjapan.R](../../../dart_check_codes/check_dart_kurilsjapan.R)). This means [gauge_summary_statistics.R](gauge_summary_statistics.R) produces 4 different R-workspace files that are [available here](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/catalog.html) - two files per event, with constant and depth-varying rigidity respectively. The latter are distinguished by having `varyMu` in the filname. On other source-zones that have DART test data, there are analogous files, including on [kermadectonga2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kermadectonga2/TSUNAMI_EVENTS/plots/catalog.html),
[kurilsjapan](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/catalog.html), 
[newhebrides2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/newhebrides2/TSUNAMI_EVENTS/plots/catalog.html), 
[puysegur2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/plots/catalog.html), 
[solomon2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/solomon2/TSUNAMI_EVENTS/plots/catalog.html), 
[southamerica](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/plots/catalog.html), and
[sunda2](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/plots/catalog.html).

Going back to the `kurilsjapan` example: if we download any of the `*.Rdata` files [from this directory](https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/plots/catalog.html) (say `gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata`), then it can be loaded from within R using:
    
    # Here we just pick one file as an example
    load('./gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata')

and you can see that many variables have been defined (use the `ls()` command to show all variables in the workspace). The variables correspond to those created by [gauge_summary_statistics.R](./gauge_summary_statistics.R). 

For comparing the modelled and observed tsunami, the uniform-slip model results are in `uniform_slip_stats`, the variable-area-uniform-slip model results are in `variable_uniform_slip_stats`, and the heterogeneous-slip model results are in `stochastic_slip_stats`. These variables are two-dimensional lists, with the first dimension corresponding to the DART buoy (check e.g. `names(stochastic_slip_stats)` to see this), and the second dimension corresponding to the model scenario. For example, to get the 10th scenario at the 2nd DART buoy for the heterogeneous-slip model, we would need to look inside `stochastic_slip_stats[[2]][[10]]`. 

The latter is itself a list (output of the function `plot_model_gauge_vs_data_gauge`). It contains multiple variables - to see their names, use the R command `names(stochastic_slip_stats[[2]][[10]])`. They include the observed time-series (`data_t` and `data_s` giving the times and stage-residuals respectively, over a time-period which focusses on the first few hours of tsunami when high-frequency measurements exist), the modelled time-series (`model_t` and `model_s`, limited to similar times as the data), and some other statistics.  

Please look carefully at the function `plot_model_gauge_vs_data_gauge` inside [gauge_summary_statistics.R](./gauge_summary_statistics.R) to understand what everything is. In particular, the time-ranges used for various statistics are rather involved, and there are some documentation errors on the approach in our reports. 

* Regarding the stored time-series and the stage-range statistic, the start time (seconds after the earthquake) is initially taken as MAX(0, MIN("time at which model exceeds 0.5% of its maxima", "10 minutes before observed maxima")), but is also limited to the time at which the DART buoy measurements have 1 min or faster frequency. The end-time is not more than 12 hours after the start-time (incorrectly stated as 3 hours in our reports), but is likewise limited to the high-frequency DART measurements. In practice, the high-frequency period generally limits the time range to a few hours. 

* For the goodness-of-fit type statistic, we use a possibly shorter time-interval -- the end-time is at most 3 hours after the observed maxima in `data_s` (incorrectly described as "at most 3 hours" in our reports, rather than "at most 3 hours after the observed maxima"). Note it is sensible to use a limited time-period for this goodness-of-fit statistic because it is sensitive to any later-time phase-errors in the model, which are common even if the overall wave size is well represented (i.e. due to limitations of the model physics and input data, and differences between the random earthquakes and the real earthquake). This is one reason that tsunami-based finite-fault inversions often only consider the first hour or two of the tsunami. In contrast the stage-range statistic is robust to such phase-errors, so there is less motivation to restrict the time interval.

## Details of event_plot.R

The [event_plot.R](event_plot.R) script creates model-vs-data time-series plots for scenarios processed by [gauge_summary_statistics.R](gauge_summary_statistics.R). The first numeric argument (e.g. 5) gives the number of hours after tsunami arrival to plot at each DART. The second numeric argument (e.g. 7.5) will exclude scenarios that have peak-slip greater than 7.5 times the mean-scaling-relation-slip inferred from the magnitude. See Section 3.2.3 in the [PTHA18 Report]() for discussion of peak-slip limits, which explains why PTHA18 uses the 7.5 factor. In reality there is much uncertainty around this limit, because slip-maxima are a poorly resolved aspect of earthquake-slip inversions. 

