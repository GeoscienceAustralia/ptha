library(sf)

# A reference scenario used to get elevation for all models -- created with `copy_elevation_rasters_locally.R`
scenario_elevation_raster_file = '../../analysis_scenarios_ID710.5/jatwc_to_inundation/elevation_in_model/all_elevation_combined.vrt'
all_elevation_tifs_dir = dirname(scenario_elevation_raster_file) # Must contain tif files with elevation for every domain 

# Tar files containing max-stage data for all scenarios -- we have to extract
# the max-stage tifs from each, and build a vrt
all_scenario_raster_tars = normalizePath(Sys.glob(c(
    '../../swals/OUTPUTS/ptha18-NSW2023b-ID710.5-sealevel110cm/random_*/ptha*/raster_output_files.tar',
    '../../swals/OUTPUTS/ptha18-NSW2023-ID1315.5-sealevel110cm/random_*/ptha*/raster_output_files.tar',
    '../../swals/OUTPUTS/ptha18-NSW2023-ID4186.3-sealevel110cm/random_*/ptha*/raster_output_files.tar'
    )))

# When computing the inundation zones, the zone is limited areas that are
# inundated at a given exceedance-rate.
#
# The exceedance-rate defining the upper limit of the inundation zone.
PTHA_EXRATE_TOL = 1/2500  # At least 1/2500
# The rasters used to determine the area inundated within this exceedance-rate
ptha_exceedance_rate_rasts = normalizePath(Sys.glob('../probabilistic_inundation/ptha18-NSW2023-MIS-sealevel110cm/highres_depth_epistemic_uncertainty/84pc/random_sum_of_source_zones-depth_exrate_0.001_0.84_sum_of_source_zones/*.tif'))

# A name for the PTHA rasters above. Used in output folder name to give
# information on PTHA raster used. Avoid spaces
ptha_exceedance_rate_rast_tag = '84pc' # 'Logic_tree_mean'

# Scenario sea-level when there is no tsunami -- needed to convert stage to
# JATWC_H values
SCENARIO_AMBIENT_SEA_LEVEL = 1.1

# For some models, we use a lower initial stage inside polygons defining 
# waterbodies that are not connected to the coast. At these sites our calculations
# need to use an alternative ambient sea level when converting raster zones to polygons
# with PTHA limits
LOWER_AMBIENT_SEA_LEVEL_IN_POLYGONS = FALSE
if(LOWER_AMBIENT_SEA_LEVEL_IN_POLYGONS){
    #
    # Example from the midwest model in WA
    #
    LOWER_AMBIENT_SEA_LEVEL_ZONES = list()
    # Add the first polygon/sea-level pair
    LOWER_AMBIENT_SEA_LEVEL_ZONES[[1]] = list(
        polygon = vect('../../initial_stage/shapes/leeman.shp'),
        sealevel = -10.0)
    # Add additional polygon/sea-level pairs
    LOWER_AMBIENT_SEA_LEVEL_ZONES[[2]] = list(
        polygon = vect('../../initial_stage/shapes/sandy_cape.shp'),
        sealevel = -1.0)
    LOWER_AMBIENT_SEA_LEVEL_ZONES[[3]] = list(
        polygon = vect('../../initial_stage/shapes/geraldton_shed.shp'),
        sealevel = -5.0)

    # Check that the provided sea levels are lower than the ambient sea level
    # (assumed by the code logic)
    for(i in 1:length(LOWER_AMBIENT_SEA_LEVEL_ZONES)){
        stopifnot(LOWER_AMBIENT_SEA_LEVEL_ZONES[[i]]$sealevel < SCENARIO_AMBIENT_SEA_LEVEL)
    }
}


# For computing the JATWC_H statistic, we only consider points reached by the
# tsunami, defined as sites with 
#     "max-stage" > ( max(SCENARIO_AMBIENT_SEA_LEVEL, ELEVATION) + WETTOL ).
# This is a reasonable critieria for all offshore sites in our modelling. 
# Care is needed for any special sites where we we set the initial sea-level to
# be less than SCENARIO_AMBIENT_SEA_LEVEL (e.g. WA, Vasse estuary, due to the
# flood gates that regularly maintain a sea-level of 40 cm AHD, see Martin et
# al's storm surge report). Such sites might not be picked up by the above
# criteria, even if they were reached by the tsunami. In typical cases that's exactly
# what we want (continuing the previous example, I didn't want sites in the
# Vasse estuary to be included in the JATWC_H sites). If that's not OK, then
# code modifications will be needed.
WETTOL = 0.001

# For the purpose of mapping the inundation zones, ignore grids in our model if
# the cell size is greater than this (arc-minutes). This prevents us creating
# results that are based on coarse resolution domains.
skip_if_cellsize_above_threshold = 0.9/(60*9*3) 
# NOTE: For the NSW model, the above choice will only catch "high-res" tiles.
# That's good because we didn't yet compute 84th percentile exrate results for other tiles.
# If this threshold were larger, and we tried to compute results where the 84th percentile tile
# didn't exist, then I'm not sure what would happen.

# JATWC MOST model (T2) excludes depth shallower than this (in meters)
jatwc_depth_limit = 20 
jatwc_dlon_dlat = 1/15 # 4 arcmin cell-size in JATWC model

# ATWS Zones
atws_zones_file = normalizePath('ATWS_ZONES/ATWS_Zones_V2_2_4/ATWS_Zones_V2_2_4.shp')
atws_zones = read_sf(dsn=atws_zones_file, layer=gsub('.shp', '', basename(atws_zones_file), fixed=TRUE))
# Was missing a projection, fix that
st_crs(atws_zones) = st_crs('EPSG:4326')

# To determine whether a domain is in a given ATWS zone, we buffer the ATWS polygon by
# "warning_zone_buffer_degrees" and then select all touching domains. The idea is to
# prevent onshore domains from being skipped if they do not touch the ATWS Zone polygon.
# A potentially negative side-effect is that some "along-coast" domains will be included
# in multiple ATWS Zone output products
warning_zone_buffer_degrees = 0.15 # NSW model needs this to be more generous than in WA

# Greenslade et al (2020) discuss relation between warning categories and the model wave height 
# statistic (termed 'JATWC_H' in this script):
#     Greenslade, D. J. M.; Uslu, B.; Allen, S. C. R.; Kain, C. L.; Wilson, K. M. &
#     Power, H. E. Evaluation of Australian tsunami warning thresholds using
#     inundation modelling Pure Appl. Geophys., Springer Science and Business Media
#     LLC, 2020, 177, 1425-1436
# 
JATWC_H_ranges = list(land_warning = c(0.55, 999999),
                      marine_warning = c(0.2, 0.55),
                      no_threat = c(-1, 0.2), 
                      major_land_warning = c(1.5, 999999),
                      minor_land_warning = c(0.55, 1.5))

# How many cores for parallel computation?
DEFAULT_MC_CORES = 104

# Fix for imperfectly recorded missing data values in rasters -- treat any
# value below this as missing data. Applied in "map_threat_levels_in_zone"
raster_na_below = -3e+38
