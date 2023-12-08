#
# Create the metadata in the OUTPUT_DIR/
# Easiest to do this programatically so I don't need to manually create it every time something is changed.
#
make_OUTPUT_DIR_README<-function(){
    source('global_variables.R')

    README_FILENAME = paste0(OUTPUT_DIR, 'README.md')

    # Text string containing the README, and a bunch of __FLAGS__ denoting folder/file names, which will be replaced below
    README_TEMPLATE = 'This folder contains the post-processed data and plots derived from the data in the "original" folder, as well as metadata tables describing the gauge locations and linking the filenames of post-processed and original files. 

Subfolders:
* __TIDE_GAUGE_OUTPUT_FOLDER__ - Tide gauge data
* __MSLP_SENSOR_OUTPUT_FOLDER__ - Mean sea level pressure data
* __GRAPHICAL_CHECKS_FOLDER__ - Figures generated for quality control purposes, and for the paper.

Metadata tables:
* __TIDEGAUGE_METADATA_TABLE_FILE__
* __MSLP_METADATA_TABLE_FILE__

List of stations that are included in the original metadata files, but were skipped in post-processing (due to missing or noisy time series):
* __IGNORED_MSLP_FILE__
* __IGNORED_TIDEGAUGE_FILE__

To see details on manual cleaning applied to some tide-gauges, see the file:
    "../scripts_to_postprocess/plot_tide_gauge_time_series_and_isolate_short_period_waves.R"
under the comment: 
    "# CLEANING OF TIDE GAUGES"
'

    README_TEMPLATE = gsub('__TIDE_GAUGE_OUTPUT_FOLDER__', basename(OUTPUT_TIDE_DIR), README_TEMPLATE)
    README_TEMPLATE = gsub('__MSLP_SENSOR_OUTPUT_FOLDER__', basename(OUTPUT_MSLP_DIR), README_TEMPLATE)
    README_TEMPLATE = gsub('__GRAPHICAL_CHECKS_FOLDER__', basename(OUTPUT_GRAPHICS_DIR), README_TEMPLATE)
    README_TEMPLATE = gsub('__TIDEGAUGE_METADATA_TABLE_FILE__', basename(TIDEGAUGE_METADATA_TABLE_FILE), README_TEMPLATE)
    README_TEMPLATE = gsub('__MSLP_METADATA_TABLE_FILE__', basename(MSLP_METADATA_TABLE_FILE), README_TEMPLATE)
    README_TEMPLATE = gsub('__IGNORED_MSLP_FILE__', basename(IGNORED_MSLP_FILE), README_TEMPLATE)
    README_TEMPLATE = gsub('__IGNORED_TIDEGAUGE_FILE__', basename(IGNORED_TIDEGAUGE_FILE), README_TEMPLATE)

    writeLines(README_TEMPLATE, README_FILENAME)
    return(invisible(0))
}

