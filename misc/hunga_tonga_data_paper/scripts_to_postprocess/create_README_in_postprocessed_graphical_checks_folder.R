#
# Create the metadata in the OUTPUT_DIR folder.
# Easiest to do this programatically so I don't need to manually create it every time something is changed.
#
make_OUTPUT_GRAPHICS_DIR_README<-function(){
    source('global_variables.R')

    README_FILENAME = paste0(OUTPUT_GRAPHICS_DIR, 'README.md')

    # Text string containing the README, and a bunch of __FLAGS__ denoting folder/file names, which will be replaced below
    README_TEMPLATE = 'This folder contains figures that were created by the post-processing scripts in the data repository.

Most of these figures are included in the manuscript (standalone or after combining with other figures). Some additional figures are not included, but were useful for quality control, or were of general interest. They are described below.
* [nearby_tide_gauges_comparison.pdf](nearby_tide_gauges_comparison.pdf) compares pairs of tide gauges that were less than 350 m from each other. For each site the plot contains two panels. The top panel shows both time series. The bottom panel shows both high-pass filtered time series near the HTHH explosion, with a small vertical offset to separate them.
* [highpass_filtered_plot_pressure_gauges.pdf](highpass_filtered_plot_pressure_gauges.pdf) includes a 3 panel plot for every MSLP station. The top panel shows the MSLP data. The middle panel shows the high-pass filtered MSLP with vertical lines (brown/blue) corresponding to the theoretical Lamb wave passage times, and a vertical red line showing the HTHH explosion time (with calculation details as described in the manuscript). The bottom panel is a zoom of the middle panel around the time of the initial Lamb wave.
* [highpass_filtered_plot_tide_gauges_before_editing_out_artefacts.pdf](highpass_filtered_plot_tide_gauges_before_editing_out_artefacts.pdf) includes a 3 panel plot for every tide gauge. The top panel shows the tide gauge data before fixing artefacts as described in the manuscript. The middle panel shows the high-pass filtered sea level, with vertical lines showing the explosion time and theoretical Lamb wave arrival times (as described for the previous figure). The bottom panel is a zoom of the high-pass filtered sea level near the start of the tsunami, following the HTHH explosion.
* [tide_gauges_comparison_multiple_detiding_techniques.pdf](tide_gauges_comparison_multiple_detiding_techniques.pdf) includes a 4 panel plot for every tide gauge. The top panel compares the observed data with an empirical "tide" obtained by subtracting the high-pass filtered sea level from the original data (so it includes tidal and non-tidal signals with period longer than 3 hours). The second panel shows the high-pass filtered sea level. The third panel compares the observed stage with an astronomical tidal prediction derived from TPXO9v5a. The fourth panel overplots the high-pass filtered sea level derived with two different detiding techniques as described in the paper (in practice it is difficult to see differences between them).
* [pressure_arrival_time_peak_and_trough.png](pressure_arrival_time_peak_and_trough.png) shows the time between the HTHH explosion and the peak (black) and trough (blue) of the high-pass MSLP residual corresponding to the initial passage of the Lamb wave.
* [distance_tonga_and_pressure_maxima_arrival.png](distance_tonga_and_pressure_maxima_arrival.png) plots the distance from the HTHH volcano against the time of arrival of the initial Lamb wave pressure maxima.
* [distance_tonga_and_pressure_maxima_size.png](distance_tonga_and_pressure_maxima_size.png) plots the distance from the HTHH volcano against the size of the initial Lamb wave pressure maxima.
* [lamb_wave_time_from_peak_to_trough.png](lamb_wave_time_from_peak_to_trough.png) plots the time between the peak and trough of the initial Lamb wave.

'
    writeLines(README_TEMPLATE, README_FILENAME)
    return(invisible(0))
}
