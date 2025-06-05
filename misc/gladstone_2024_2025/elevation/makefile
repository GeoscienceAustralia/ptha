# This makefile is used to generate the elevation files for the SWALS model.
#
# It is called by the main makefile in the root directory of the repository.

# The MSQ bathymetry file names
dir_msg_bathy = /g/data/w85/tsunami/DATA/ELEVATION/QLD/MSG_Gladtstone_Port/Gladstone_Data_full_density
msg_bath_name = GL020120_Gladstone_Historic_data_Interpolated_Merged_2p5m_Cube_FINAL

# Flat patches
dir_gladstone_patches = /g/data/w85/tsunami/MODELS/inundation/DATA/ELEV/QLD/gladstone_patches
SHAPEFILES = $(wildcard $(dir_gladstone_patches)/shapes/*)
FLAT_TIF_FILES = $(patsubst $(dir_gladstone_patches)/shapes/%,$(dir_gladstone_patches)/flat_rasters/%.tif,$(basename $(SHAPEFILES)))
TOP_FLAT_TIF_FILES = $(wildcard $(dir_gladstone_patches)/top_flat_rasters/*.tif)

# List all the elevation files, incomplete...
ELEV_DEPS = elevation/open_estuary_entrances/estuary_entrance_polygons_values.csv \
	elevation/swals_elevation_files_in_preference_order.txt \
    $(dir_msg_bathy)/$(msg_bath_name)_AHD_wgs.tif \
	$(FLAT_TIF_FILES) \
	$(TOP_FLAT_TIF_FILES)


elevation: $(ELEV_DEPS)

# Trick so that anything that depends on force target will always be remade
force: ;
.PHONY: force

# Convert MSQ bathymetry to a tif
$(dir_msg_bathy)/$(msg_bath_name)_AHD_wgs.tif: $(dir_msg_bathy)/convert_bag_to_tif.R $(dir_msg_bathy)/full_density_dataset/$(msg_bath_name).bag
	$(MODULES) && cd $(dir_msg_bathy) && Rscript convert_bag_to_tif.R

# Flat patches manually made for the model
$(dir_gladstone_patches)/flat_rasters/%.tif: $(dir_gladstone_patches)/shapes/%.*
	$(MODULES) && cd $(dir_gladstone_patches) && Rscript make_patch.R

# Copy the top flat rasters to the correct directory
$(dir_gladstone_patches)/top_flat_rasters/%.tif: $(dir_gladstone_patches)/flat_rasters/%.tif $(dir_gladstone_patches)/move_top_rasters.sh
	cd $(dir_gladstone_patches) && . ./move_top_rasters.sh

# force the elevation script to always run to ensure the preference list is up to date
elevation/swals_elevation_files_in_preference_order.txt: $(FLAT_TIF_FILES) $(TOP_FLAT_TIF_FILES) $(dir_msg_bathy)/$(msg_bath_name)_AHD_wgs.tif elevation/make_swals_elevation_files_preference_list.R force
	module load R && cd elevation && Rscript make_swals_elevation_files_preference_list.R

elevation/open_estuary_entrances/estuary_entrance_polygons_values.csv:
