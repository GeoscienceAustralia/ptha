# site_name is used in various 'naming' situations for the model runs and
# directory structure
site_name = basename(dirname(getwd()))

# This should expand to a vector of all the initial stage 'tif' filenames you
# want to run
initial_condition_files = normalizePath(
    Sys.glob(paste0('../EQ_SOURCE/Unit_source_data/', site_name, '/', site_name, '*.tif')))

# The runs will occur inside this folder (which will be created)
all_runs_dir = 'unit_source_tsunami'

# The output will go inside here (with sub folders corresponding to
# all_runs_dir/run_initial_condition/)
all_runs_output_base_dir = paste0('./OUTPUTS/', site_name)
# all_runs_output_base_dir = paste0('/g/data/w85/tsunami/AustPTHA/version1/unit_sources/', site_name)
# all_runs_output_base_dir = paste0('/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/', site_name, '/TSUNAMI_UNIT_SOURCES/unit_source_tsunami/')

# Elevation grid for tsunami model runs
elevation_data_file = normalizePath('../../../DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif')

# Hazard points file for tsunami model runs
hazard_points_file = normalizePath('../../../DATA/HAZARD_POINTS/merged_hazard_points_blocked.csv')

# Basic sanity checks
stopifnot(file.exists(elevation_data_file))
stopifnot(file.exists(hazard_points_file))
stopifnot(length(initial_condition_files) > 0)
stopifnot(all(file.exists(initial_condition_files)))
