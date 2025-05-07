# Plots and statistics to compare stochastic earthquake tsunami models with data

Run with:
```
# 1. Compute statistics, takes hours and uses lots of memory
Rscript parse_gauge_outputs.R > parse_gauge_outputs_log.log 

# 2. Plot time-series of model-vs-data at selected gauges (based on the tsunami size)
Rscript plot_best_scenarios_large_and_small_waves.R > plot_best_scenarios_log.log 

# 3. Make plots of "good" scenarios according to various GOF criteria (collapsed over tide gauges, e.g. by taking the median GOF over tide gauges), along with some GOF stats
Rscript analysis_good_scenarios_plot.R

# 4. Rank the different model types based on their GOF collapsed over tide gauges
Rscript summary_statistics_median_GOF_over_gauges.R > summary_statistics_median_GOF_over_gauges_log.log

# 5. Make boxplots showing random models vs observations for various tsunami size metrics
Rscript boxplot_stats_and_scatterplots.R

# 6. Compute tables of statistics
Rscript run_stats_under_null_hypothesis.R
```

## Details

* The parsing script may report errors related to tide-gauges which are included in the `gauge_data_links.R` database but not used in the analysis. That's OK and more convenient than modifying the database.
* The CPU time required for all models can be computed with `sum_walltime_all_tsunami_models.R`
* The time window used for comparison of models and data was computed using the function `.make_times_to_start_comparison()` in the script `create_times_to_start_tide_gauge_comparison_with_models.R`. 
  * This requires preliminary versions of files from `parse_gauge_outputs.R` to exist first 
  * In practice the latter was run once, then the time window was created, and then the parsing was redone.
  * Here I've provided the time window output file directly, so this isn't necessary.
* The gauges that are used was counted with `get_coordinates_of_good_gauges.R`.
* Some information on durations of gauges that were truncated is in `gauge_truncation_info.R`
