#
# Inputs for both the 'earthquake_types_plot.R', and 'plot_stage_time_series.R'
#
library(rptha)
# Get the ptha-access routines
ptha=new.env()
source('../../../../CODE/ptha/ptha_access/get_PTHA_results.R', chdir=TRUE, local=ptha)
# Local routine that is helpful
source('./find_unit_sources_near_hypocentre.R')

#
# Specify the source-zone, and provide links to files that describe it. These must be set before running any scripts.
#

# The source-zone name -- must match the name used in PTHA18
source_zone = 'kermadectonga2'
# The PTHA18 unit-sources shapefile. This can be downloaded from the NCI
# THREDDS SERVER. The location is like this (using the example of alaskaaleutians -- other sources are analogous):
# https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/EQ_SOURCE/unit_source_grid/catalog.xml
unit_source_polygon_shapefile = 
    '../../../../CODE/ptha/ptha_access/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/unit_source_grid/kermadectonga2.shp'
# Files with earthquake events for stochastic/variable-uniform/uniform slip, and the unit-source statistics
# Download these from the NCI THREDDS SERVER and place in the current
# directory (or edit the file name to include its full path). The download
# location is like this (for alaskaaleutians -- other sources are analogous):
# https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/catalog.xml
stochastic_slip_events = 'all_stochastic_slip_earthquake_events_kermadectonga2.nc' 
uniform_slip_events = 'all_uniform_slip_earthquake_events_kermadectonga2.nc' 
variable_uniform_slip_events = 'all_variable_uniform_slip_earthquake_events_kermadectonga2.nc' 
unit_source_statistics = 'unit_source_statistics_kermadectonga2.nc'

# A background DEM in lon/lat coordinates that covers the "plot_region" defined below. Any DEM will do
dem = raster('../../../../MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')
# A background zero-contour shapefile -- just make it from your DEM (e.g. using GDAL, something like "gdal_contour -fl 0.0 my_dem_name.tif zero_contour")
zero_contour = '../../../../MODELS/AustPTHA_c/DATA/ELEV/merged_dem/zero_contour/zero_contour.shp'

# Output location for plots
out_dir = 'plots/source_animation'

#
# Parameters affecting the location of the selected events on the source-zone. These must be set before running any scripts.
#

# Reference earthquake location and magnitude
# The user must make sure this point touches the nearest unit-source (move if needed)
event_hypocentre = c(187.4, -15.49 ) 
event_mw = 8.1
# Accept events within this magnitude difference of event_mw. Avoid 0.0 because
# of possible floating-point-storage imperfections.
magnitude_difference = 0.05 

#
# Plot parameters for the earthquake/uplift plots. These must be set before running any scripts.
#

# Region of plot, or NULL for selection based on unit-sources
plot_region = extent(180, 195, -25, -10) # NULL

# Max/min of surface deformation scale
SURFACE_DEF_MIN = -3
SURFACE_DEF_MAX = 5
# Max slip scale
MAX_SLIP_VAL = 20
# Z-range for DEM
DEM_ZLIM = c(-10000, 0)

#
# Plot parameters for the stage time-series plots. These can be changed after running earthquake_types_plot.R, but must
# be set before running plot_stage_time_series.R
#

# Coordinates (lon/lat) of sites at which we plot stage time-series. There
# should be exactly 2 sites, or none (denoted with NULL). The points better be close
# to the PTHA18 hazard points.
#wave_sites = NULL
wave_sites = rbind( c(183.62, -9.46), 
                    c(191.59, -23.07))

# Plot parameters for stage time-series
site1_yrange = 0.2 # Y-limits of stage time-series plot in m will be [-site1_yrange, site1_yrange] 
site1_tstart = 0*3600 # Start time of stage time-series plot in seconds  (ie.. min x-axis limit)
site1_tend   = 6*3600 # End time of stage time-series plot in seconds (i.e. max x-axis limit)
site2_yrange = 0.2 # Like site1_yrange for site2
site2_tstart = 0*3600 # Like site1_tstart for site2
site2_tend   = 6*3600 # Like site1_tend for site2

# Rows in the event table to plot. If NULL, instead use the 'stage_plot_desired_count'
stage_plot_desired_rows = NULL # 49321:49770 # Set this if you know the exact rows to use. The row-numbers correspond to the full event table (i.e. with all magnitudes), so this is quite inconvenient.
stage_plot_desired_count = 100 # Plot at most this many events -- thin the events for which we plotted earthquakes if required. This is usually what you want to do. To plot all events with an earthquake plot, just set this to a large number (e.g. 1e+06). 

# DART buoy files to overplot -- or NULL if there is no data.
#dart_data_site1 = 'samoa_2009_09_29_Mw8.1_51425.csv' # NULL
#dart_data_site2 = 'samoa_2009_09_29_Mw8.1_51426.csv' # NULL
#event_time_GMT = '2009-09-29 17:48:11'
dart_data_site1 = NULL
dart_data_site2 = NULL
event_time_GMT = NULL

#
# Plot in parallel, using this many cores
#
library(parallel)
MC_CORES = detectCores()


#
# Check
#

# 2 wave sites only
if(!all(is.null(wave_sites))){
    stopifnot(nrow(wave_sites) == 2)
}
