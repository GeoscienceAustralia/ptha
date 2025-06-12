#
# Load data used to make various plots -- we source this directly in other
# scripts, so be careful about changing names.
#

library(terra)

all_events = read.csv('target_scenarios_data_frame.csv')
all_events_focal_mech = read.csv('target_scenarios_earthquake_catalogue_mechanisms_for_plotting_only.csv', 
    comment.char='#')

# Shapefile for every PTHA18 unit-source
all_unit_source_grids = Sys.glob(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/SOURCE_ZONES/*/EQ_SOURCE/unit_source_grid/*.shp'
    )
# Read unit sources
uss = lapply(all_unit_source_grids, vect)
names(uss) = basename(dirname(dirname(dirname(all_unit_source_grids))))

# Coastline
zero_contour = vect(
    '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/zero_contour/zero_contour.shp'
    )

# Countries
tm_wb = vect('/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/DOC/TESTING/model_data_comparison/FIG/DART_locations_test_events/TM_WORLD_BORDERS_PACIFIC/')


