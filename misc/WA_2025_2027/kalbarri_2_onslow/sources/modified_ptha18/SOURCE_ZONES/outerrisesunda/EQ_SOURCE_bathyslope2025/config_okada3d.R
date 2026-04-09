#
# Configuration parameters for unit-sources on this source-zone
#

# Name of the source-zone of interest, using the source-zone names from the
# 2018 Australian Probabilistic Tsunami Hazard Assessment
site_name = 'outerrisesunda' 

# Store Okada easting/northing/z displacements on a raster with lon/lat extent
# derived from a buffered version of the depth contours extent. The buffer is
# specified below (and then in practice is increased so the extents are integers).
# The computational time will scale with the number of cells in this region
output_raster_ext_buffer = 2 # degrees
output_raster_cellsize = 2/60 # degrees

# csv file with parameters for source-zone. We will use this to get the rake
# for the source-zone unit-sources. Same as used for 2018 PTHA, see here for a version controlled variant:
#    https://github.com/GeoscienceAustralia/ptha/blob/master/R/examples/austptha_template/DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv
sourcezone_parameter_file = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv'
stopifnot(file.exists(sourcezone_parameter_file))

# Location for outputs
output_base_dir = './'
dir.create(output_base_dir, recursive=TRUE, showWarnings=FALSE)

# Directory with sourcezone contours used in Australian PTHA. Accessed from here (zip-file): 
#    https://thredds.nci.org.au/thredds/catalog/fj6/PTHA/AustPTHA_1/DATA/catalog.html
SOURCEZONE_CONTOURS_BASEDIR = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/SOURCEZONE_CONTOURS/'
# Directory with sourcezone downdip lines used in Australian PTHA. Accessed
# from the same location as SOURCEZONE_CONTOURS_BASEDIR (above)
SOURCEZONE_DOWNDIP_LINES_BASEDIR = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/SOURCEZONE_DOWNDIP_LINES/'

# Elevation raster used to compute bathymetric gradients for vertical motion
# due to horizontal components, and for Kajiura filtering.
elevation_raster_file = '/g/data/w85/tsunami/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
# We need (elevation_rater * elevation_raster_scale) to give the elevation in
# units of meters above sea level (so bathymetry is negative).
elevation_raster_scale = 1.0 
# Is the elevation raster in geographic spherical coordinates?
isLonLat_elevation_raster = isLonLat(raster(elevation_raster_file))
stopifnot(isLonLat_elevation_raster) # We need geographic spherical

# A vector with shapefile names for all contours that we want to convert to
# unit sources
all_sourcezone_shapefiles = paste0(SOURCEZONE_CONTOURS_BASEDIR, '/', site_name, '.shp')
all_sourcezone_downdip_shapefiles =  paste0(SOURCEZONE_DOWNDIP_LINES_BASEDIR, '/', site_name, '_downdip.shp')
stopifnot(length(all_sourcezone_shapefiles) == 1)
stopifnot(length(all_sourcezone_downdip_shapefiles) == 1)

# Get the rake and desired length/width from a file
all_sourcezone_par = read.csv(sourcezone_parameter_file, header=TRUE, stringsAsFactors=FALSE)
source_rows = which(all_sourcezone_par$sourcename == site_name)[1]

# Desired unit source geometric parameters
desired_subfault_length = as.numeric(all_sourcezone_par$approx_unit_source_length[source_rows]) # km
desired_subfault_width = as.numeric(all_sourcezone_par$approx_unit_source_width[source_rows]) # km

# A vector with the desired rake angle (one entry per sourcezone)
sourcezone_rake = rep(as.numeric(all_sourcezone_par$rake[source_rows]),
    len=length(all_sourcezone_shapefiles)) # degrees

# Desired spacing of sub-unit-source points
# Lower values (e.g. 1000) may be required for accuracy in unit sources
# near the trench, because shallow deformation tends to be quite localised.
# For deeper unit sources, a much coarser point spacing can be used without
# sacrificing accuracy. 
# Hence we use different values for the 'shallow' sub-unit-source points (i.e.
# < 50km down dip) and the deeper ones.
# The computational effort approximately scales with the inverse square of
# the point density. 
shallow_subunitsource_point_spacing = 600 # m
deep_subunitsource_point_spacing = 4000 #m

# Taper edges of unit_source slip with circular filter having this radius (m)
# This can be useful to avoid features of the Okada solution associated with
# slip discontinuities at the rupture edges. 
# E.G. For ruptures with shallow (but non-zero) top depth, the Okada solution
# suggests a high 'ridge' of deformation just above the top-edge, which is
# entirely due to the discontinuity in the slip. Slip tapering will smooth out
# such features.
slip_edge_taper_width = 10000

# For computational efficiency, only compute the okada deformation at
# distances <= okada_distance_factor x (depth of sub-unit-source point) 
# This can save computational effort for shallow unit sources.
# But be careful if using a wide subunitsource_point_spacing.
okada_distance_factor = 20 # Inf 

# Spatial scale for sub-cell point integration
# During the Okada computation, points with "abs(deformation) > 10% of max(abs(deformation))"
# will have deformations re-computed as the average of the 16 Okada point values
# around point p. These 16 points have coordinates:
#     points = expand.grid(p[1] + cell_integration_scale[1]*c(-1,-1/3,1/3,1), 
#                          p[2] + cell_integration_scale[2]*c(-1,-1/3,1/3,1))
# If 'cell_integration_scale' is close to half the grid size, then this is an approximation
# of the within-pixel average Okada deformation. We do this because near the trench,
# the Okada deformation might not be smooth [e.g. when rupture depth --> 0], and this
# reduces the chance of artificial 'spikes' in the Okada deformation.
# In the code below, this is only applied along the 'top' row of unit-sources
# where the trench depth might --> 0.
cell_integration_scale = c(1500, 1500)

# Number of cores for parallel parts. Values > 1 will only work on shared
# memory linux machines.
MC_CORES = 104 #16 #detectCores()
