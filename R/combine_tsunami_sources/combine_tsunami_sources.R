# Combine unit source tsunami initial conditions to make tsunami initial
# conditions for earthquake events.
library(rptha)
fse = new.env()
source('fit_simulate_earthquake.R', local=fse)

## Input parameters ###########

# Folder containing one directory for each sourcename. Inside the latter
# directories are tif files for each unit source (and no other tif files)
unit_source_dirname = '../source_contours_2_unit_sources/Unit_source_data'

# sourcename. This should be the name of a directory inside
# unit_source_dirname, and also the name of a discretised_source (among those
# contained in all_discretized_source_RDS)
sourcename = 'alaska'

# RDS filename containing all discretized source information. The object
# therein should include a list entry corresponding to sourcename
all_discretized_source_RDS = 
    '../source_contours_2_unit_sources/all_discretized_sources.RDS'

# Earthquake parameters
desired_Mw = 9.0
mu = 3e+10
#
# Corner wavenumbers come from S_{NCF} regression relations in Davies et al.
# (2015)
#
physical_corner_wavenumbers = 10**c(-0.54*desired_Mw + 2.03, -0.41*desired_Mw + 1.18) 

## end input ##################


discretized_source = readRDS(all_discretized_source_RDS)[[sourcename]]
# To make a stochastic slip field, we need to get dimensions for the unit sources
# This can be estimated from the statistics
discretized_source_statistics = discretized_source_approximate_summary_statistics(
    discretized_source)

nx = discretized_source$discretized_source_dim['strike']
ny = discretized_source$discretized_source_dim['dip']

# Read the raster corresponding to each row in the discretized_source_statistics
unit_source_raster_files = paste(unit_source_dirname, '/', sourcename, '/', 
    sourcename, '_', discretized_source_statistics$downdip_number, '_',
    discretized_source_statistics$alongstrike_number, '.tif', sep="")

if(!(all(file.exists(unit_source_raster_files)))){
    stop('Could not find some unit source raster files')
}

unit_source_rasters = lapply(as.list(unit_source_raster_files), f<-function(x) raster(x))

# Get average dx/dy for unit sources, where dx is along-strike and dy is
# down-dip
mean_dx = mean(discretized_source_statistics$length)
mean_dy = sum(discretized_source_statistics$width * 
    discretized_source_statistics$length) / 
    sum(discretized_source_statistics$length)

# Record full dx/dy in matrices
dx = matrix(NA, ncol=nx, nrow=ny)
dy = matrix(NA, ncol=nx, nrow=ny)
for(i in 1:nrow(discretized_source_statistics)){
    nr = discretized_source_statistics$downdip_number[i] 
    nc = discretized_source_statistics$alongstrike_number[i]
    dx[nr, nc] = discretized_source_statistics$length[i]
    dy[nr, nc] = discretized_source_statistics$width[i]
}

#
template_slip_matrix = dx * 0
slip_matrix = dx * 0

peak_slip_row = 4
peak_slip_col = floor(nx/2)

numerical_corner_wavenumbers = physical_corner_wavenumbers * 
    c(dx[peak_slip_row, peak_slip_col], dy[peak_slip_row, peak_slip_col])

desired_M0 = M0_2_Mw(desired_Mw, inverse=TRUE)

source_info = list()
for(j in 1:3){
    # Ensure peak slip occurs in desired location
    template_slip_matrix = template_slip_matrix * 0
    template_slip_matrix[peak_slip_row, peak_slip_col] = 1
    slip_matrix = fse$simulate_sffm(numerical_corner_wavenumbers, 
        template_slip_matrix)

    # Ensure M0 is correct
    # We need slip * dx * dy * mu = M0
    initial_moment = sum(slip_matrix * dx * dy * 1e+06 * mu)
    slip_matrix = slip_matrix/initial_moment * desired_M0
    stopifnot(abs(sum(slip_matrix * dx * dy * 1e+06 * mu) - desired_M0) < (1.0e-06 * desired_M0))

    # Make a raster for nice output plots
    slip_raster = raster(slip_matrix, xmn=0, xmx=nx*mean_dx, ymn=0, ymx=ny*mean_dy)

    # Compute the tsunami source
    source_raster = unit_source_rasters[[1]] * 0
    for(i in 1:nrow(discretized_source_statistics)){
        nr = discretized_source_statistics$downdip_number[i] 
        nc = discretized_source_statistics$alongstrike_number[i]
        source_raster = source_raster + slip_matrix[nr, nc] * unit_source_rasters[[i]]
    }

    # Save the information
    source_info[[j]] = list(slip_matrix = slip_matrix, slip_raster = slip_raster, 
        source_raster = source_raster)
}

# Apply kajiura to the raster. This can be less artefact prone than applying
# it first to the unit sources, and then summing. Although mathematically the
# two are equivalent, numerically the unit source approach is much more prone to
# artefacts

# FIXME: Integrate this code into the main ptha package
source('kajiura_smooth_raster.R')
for(i in 1:length(source_info)){
    source_info[[i]]$source_raster_smooth = 
        kajiura_smooth_raster(
            source_info[[i]]$source_raster,
            new_origin=c(209, 58),
            elevation_raster_file = '../../../../DATA/ELEV/GEBCO_08/gebco_08.nc',
            kj_filter_grid_dxdy = 2000,
            kj_filter_def_threshold=1.0e-02,
            kj_cartesian_buffer = 10000,
            minimum_kj_depth = 10,
            elevation_extraction_x_offset=-360,
            spherical_input = TRUE)
}

png('Deformation_plot.png', width=15, height=9, units='in', res=300)
par(mfrow=c(2,3))
par(mar=c(3,2,1,1))
for(i in 1:3) plot(source_info[[i]]$slip_raster, xlab='Along-strike distance (km)', 
    ylab='Up-dip distance (km)', xlim=c(100, 500))
for(i in 1:3) {
    persp(source_info[[i]]$source_raster_smooth, col='skyblue', shade=TRUE, 
        border=NA, phi=30, theta=-90, scale=FALSE, expand=0.2)
}
dev.off()
