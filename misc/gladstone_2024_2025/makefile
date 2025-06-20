# Purpose: Makefile for generating the input files for the SWALS model
# 
# Run from the current directory as `make` or `make all`
#
# Remarks:
# This script does not guarantee that the model will run.
# It does not generate sources for the initial conditions (which may require compute).
# Many of the inputs are generated in GIS and only prepared for SWALS here. Inspect of the GIS files.

MODULES = modules_R_431.sh

# include the elevation makefile which has the elevation dependencies
include elevation/makefile

# The SWALS model.f90 depends on the following files
DEPS = gauges/locations.csv \
	   breakwalls/swals_breakwall_files.txt \
	   initial_stage/override_initial_stages.csv \
	   friction/swals_manning_n_files_in_preference_order.txt \
	   friction/swals_manning_n_files_in_preference_order_no_mangroves.txt \
	   inverts/swals_invert_files.txt \
	   tides/tpxo9_adjusted_for_gauges.tif \
	   multidomain_design/domains_1_3_15_15/*_edited.csv \
	   elevation

# all makes all the input files in the DEPS list and always checks for new elevation files
all: $(DEPS)

.DEFAULT_GOAL := all
.PHONY: all test elevation

# Gauge locations
gauges/locations.csv: gauges/*/* gauges/make_gauge_locations.R
	$(MODULES) && cd gauges && Rscript make_gauge_locations.R

# Breakwall list, fixed height and buffer height
breakwalls/swals_breakwall_files.txt: breakwalls/fixed_height/shapes/*.csv breakwalls/buffer_height/shapes/*.csv breakwalls/make_breakwalls_list.R
	$(MODULES) && cd breakwalls && Rscript make_breakwalls_list.R
	$(MODULES) && cd breakwalls && Rscript test_breakwalls.R

breakwalls/fixed_height/shapes/%.csv: breakwalls/fixed_height/shapes/%.shp breakwalls/fixed_height/shapes/%.prj breakwalls/fixed_height/shapes/%.dbf breakwalls/fixed_height/shapes/%.shx breakwalls/fixed_height/make_breakwalls.R
	$(MODULES) && cd breakwalls/fixed_height && Rscript make_breakwalls.R

breakwalls/buffer_height/shapes/%.csv: breakwalls/buffer_height/shapes/%.shp breakwalls/buffer_height/shapes/%.prj breakwalls/buffer_height/shapes/%.dbf breakwalls/buffer_height/shapes/%.shx breakwalls/buffer_height/make_breakwalls.R
	$(MODULES) && cd breakwalls/buffer_height && Rscript make_breakwalls.R

# Initial stages
initial_stage/override_initial_stages.csv: initial_stage/shapes/* initial_stage/shape_to_csv.R
	$(MODULES) && cd initial_stage && Rscript shape_to_csv.R

# Friction
friction/swals_manning_n_files_in_preference_order.txt: friction/data/*tif friction/data/DEA_Mangroves_friction.csv friction/data/fuels_clipped_friction.csv friction/make_friction_rasters.R
	$(MODULES) && cd friction && Rscript make_friction_rasters.R

friction/swals_manning_n_files_in_preference_order_no_mangroves.txt: friction/data/*.tif friction/data/fuels_clipped_friction_no_mangroves.csv friction/make_friction_rasters_no_mangroves.R
	$(MODULES) && cd friction && Rscript make_friction_rasters_no_mangroves.R

# Inverts
inverts/swals_invert_files.txt: inverts/*/* inverts/*.R
	$(MODULES) && cd inverts && Rscript make_inverts.R
	$(MODULES) && cd inverts && Rscript test_inverts.R

# Tidal data
tides/tpxo9_adjusted_for_gauges.tif: tides/x_y_z.csv tides/tpxo9_adjusted_for_gauges.R
	$(MODULES) && cd tides && Rscript tpxo9_adjusted_for_gauges.R

# Multidomain design + use GIS to edit the multi-domain design
multidomain_design/domains_1_3_15_15/nesting_level_%_edited.csv: multidomain_design/*.R multidomain_design/input/nesting_level_%.* multidomain_design/*level_domains_to_coarsen/*
	$(MODULES) && cd multidomain_design && Rscript create_boxes.R
	$(MODULES) && cd multidomain_design && Rscript edit_boxes.R

# Tests of inputs
test: all tides/test_adjustment.R tides/tpxo9_adjusted_for_gauges.tif
	$(MODULES) && cd tides && Rscript test_adjustment.R
