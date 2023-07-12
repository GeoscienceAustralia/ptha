# Make a file with mean sea level pressure gauge locations & plot them, and export the metadata to csv
Rscript plot_pressure_gauge_locations_and_write_metadata_to_csv.R

# Make a plot of the mean sea level pressure for a few days before/after the eruption,
# including panels with the short-period component (< 2h) and comparison with a simplified
# theory of Lamb-wave arrivals.
Rscript plot_pressure_time_series_and_isolate_short_period_waves.R

# Make a file with tide-gauge locations & plot them, and export the metadata to csv.
Rscript plot_tide_gauge_locations_and_write_metadata_to_csv.R

# Plot and de-tide the tide-gauge data, and export to csv. 
Rscript plot_tide_gauge_time_series_and_isolate_short_period_waves.R

# Compare gauges that are close to each other
Rscript plot_gauges_near_each_other.R

# Get some information on the gauge duration, temporal interval, etc.
Rscript -e "source('compute_gauge_temporal_interval.R', echo=TRUE)" > GAUGE_TEMPORAL_INTERVAL_SUMMARY.txt

# Make some plots
Rscript plots_for_paper.R
# Make a few 'single-file' images from multiple images
Rscript merge_some_figures.R 
