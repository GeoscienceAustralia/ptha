#
# Configuration parameters for unit-sources on this source-zone
#

site_name = basename(dirname(getwd()))

# csv file with parameters for source-zone. We will use this to get the rake
# for the source-zone unit-sources
sourcezone_parameter_file = 
    '../../../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv'
stopifnot(file.exists(sourcezone_parameter_file))

output_base_dir = paste0('/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/',
    site_name, '/EQ_SOURCE/')
dir.create(output_base_dir, recursive=TRUE, showWarnings=FALSE)

# A vector with shapefile names for all contours that we want to convert to
# unit sources
all_sourcezone_shapefiles = paste0(
    '../../../DATA/SOURCEZONE_CONTOURS/', site_name, '.shp')
all_sourcezone_downdip_shapefiles =  paste0(
    '../../../DATA/SOURCEZONE_DOWNDIP_LINES/', site_name, '_downdip.shp')
stopifnot(length(all_sourcezone_shapefiles) == 1)
stopifnot(length(all_sourcezone_downdip_shapefiles) == 1)

# Desired unit source geometric parameters
desired_subfault_length = 50 # km
desired_subfault_width = 50 # km

# Get the rake (from a file)
all_sourcezone_par = read.csv(config$sourcezone_parameter_file, header=TRUE, 
    stringsAsFactors=FALSE)
source_rows = which(all_sourcezone_par$sourcename == site_name)[1]
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

# elevation raster (required for Kajiura filtering). Should give elevation in m, 
# with the ocean having elevation < 0. Should have a lon/lat spatial projection. 
# Set to NULL to not use Kajiura filtering.
#elevation_raster = NULL 
## A realistic example would look like:
elevation_raster = raster('../../../DATA/ELEV/GEBCO_2014_1m/GEBCO_2014_1minx1min_W-39.9958333-E320.0041667.tif')
## Note that for Kajiura filtering, a minimum depth of 10m will be assumed 
## (to avoid passing negative depths to the Kajiura smoothing routine)

# For computational efficiency, only apply Kajiura filtering in a box
# containing all points where the unit source deformation exceeds
# kajiura_use_threshold. Set to zero to apply Kajiura filter everywhere.
#
# Use of a small positive number can be faster.
# Since the unit sources have 1m slip, use of e.g. 1e-03 suggests an
# error of < 1cm to the free surface, even if the slip were 10m. 
# In practice there might be greater difference because our Kajiura routine
# involves interpolation to/from cartesian coordinates. Interpolation creates
# slight diffusion, and changes to the Kajiura box will affect the
# interpolation and so also affect this, though not in a systematic way.
kajiura_use_threshold = 1.0e-03

# When applying the kajiura filter, the data is regridded onto a grid with
# spacing=kajiura_gridspacing. The latter should be small compared to the
# horizontal distance over which the free surface deformation changes
# significantly (and small compared with the distance of
# tsunami_source_cellsize). If this is not small enough, artefacts
# can be observed especially when summing tsunami sources.
# A numerically easier alternative is to apply kajiura AFTER summing
# the sources [see script in 'combine_tsunami_sources' folder]
kajiura_grid_spacing = 500 # m

# Cell size for output rasters
# The computation time will scale inversely with tsunami_source_cellsize^2
# Here we use a relatively coarse discretization, for demonstration purposes
tsunami_source_cellsize = 2/60 # degrees. 

# Spatial scale for sub-cell point integration
# During the Okada computation, points with "abs(deformation) > 10% of max(abs(deformation)"
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
MC_CORES = 16

# Option to illustrate 3d interactive plot creation
#
# Only make the 3d interactive plot if you can use interactive graphics and
# have rgl (i.e. use FALSE on NCI). 
make_3d_interactive_plot = FALSE

# Make a multi-page pdf plot of the sources
make_pdf_plot = FALSE

# Option to reduce the size of RDS output
# TRUE should be fine for typical usage
minimise_tsunami_unit_source_output = TRUE


