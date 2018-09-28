#
# Inputs for both the 'earthquake_types_plot.R', and 'plot_stage_time_series.R'
#
library(rptha)

# Get the ptha-access routines
ptha=new.env()
source('../../../CODE/ptha/ptha_access/get_PTHA_results.R', chdir=TRUE, local=ptha)
# Local routine that is helpful
source('./find_unit_sources_near_hypocentre.R')

# The source-zone name
source_zone = 'kermadectonga2'
# The PTHA18 unit-sources shapefile
unit_source_polygon_shapefile = 
    '../../../CODE/ptha/ptha_access/SOURCE_ZONES/kermadectonga2/EQ_SOURCE/unit_source_grid/kermadectonga2.shp'
# A background DEM
dem = raster('../../../MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')
# A background zero-contour shapefile
zero_contour = '../../../MODELS/AustPTHA_c/DATA/ELEV/merged_dem/zero_contour/zero_contour.shp'

# Files with earthquake events for stochastic/variable-uniform, and uniform
# Download these from the NCI THREDDS SERVER
stochastic_slip_events = 'all_stochastic_slip_earthquake_events_kermadectonga2.nc' 
uniform_slip_events = 'all_uniform_slip_earthquake_events_kermadectonga2.nc' 
variable_uniform_slip_events = 'all_variable_uniform_slip_earthquake_events_kermadectonga2.nc' 
# Unit-source-statistics file for the source zone
# Download it from NCI THREDDS SERVER
unit_source_statistics = 'unit_source_statistics_kermadectonga2.nc'

# Output location for plots
out_dir = 'plots/source_animation'

# Reference earthquake location and magnitude
# Make sure it touches the nearest unit-source (move if needed)
event_hypocentre = c(187.4, -15.49 ) 
event_mw = 8.1
# Get all events within this magnitude difference
magnitude_difference = 0.05 

# Plot in parallel
MC_CORES=12

#
# Plot parameters
#

# Max/min of surface deformation scale
SURFACE_DEF_MIN = -3
SURFACE_DEF_MAX = 5
# Max slip scale
MAX_SLIP_VAL = 20
# Z-range for DEM
DEM_ZLIM = c(-10000, 0)

# Region of plot, or NULL for selection based on unit-sources
plot_region = extent(180, 195, -25, -10) # NULL

# Sites at which we plot stage time-series
#wave_sites = NULL
wave_sites = rbind( c(183.62, -9.46), c(191.59, -23.07))
# Stage time-series 
site1_yrange = 0.2
site2_yrange = 0.2
site1_tstart = 0
site2_tstart = 0

# Rows in the stochastic event table to plot. 
stage_plot_desired_rows = 26536:26685

# DART buoy files to overplot
dart_data_site1 = 'samoa_2009_09_29_Mw8.1_51425.csv' # NULL
dart_data_site2 = 'samoa_2009_09_29_Mw8.1_51426.csv' # NULL
event_time_GMT = '2009-09-29 17:48:11'

#
# Check
#

# 2 wave sites only
if(!all(is.null(wave_sites))){
    stopifnot(nrow(wave_sites) == 2)
}
