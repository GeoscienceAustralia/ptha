## ---- initialise ----

# Main 'driver' script to create the unit sources
#
# Gareth Davies, Geoscience Australia 2015/16
#
library(rptha)
library(raster)

###############################################################################
#
# Main input parameters 
#
###############################################################################

# A vector with shapefile names for all contours that we want to convert to
# unit sources
site_name = basename(dirname(getwd()))

all_sourcezone_shapefiles = paste0('../../../DATA/SOURCEZONE_CONTOURS/', site_name, '.shp')
all_sourcezone_downdip_shapefiles =  paste0('../../../DATA/SOURCEZONE_DOWNDIP_LINES/', site_name, '_downdip.shp')

# Desired unit source geometric parameters
desired_subfault_length = 50 # km
desired_subfault_width = 50 # km

# A vector with the desired rake angle (one entry per sourcezone)
sourcezone_rake = rep(90, len=length(all_sourcezone_shapefiles)) # degrees

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
# E.G. For ruptures with shallow top depth, the Okada solution suggests a high
# 'ridge' of deformation just above the top-edge, which is entirely due to the
# discontinuity in the slip. Slip tapering will smooth out such features.
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
# If 'cell_integration_scale' is close to the grid size, then this is an approximation
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


## ---- takeCommandLineParameter ----

if(interactive() == FALSE){

    #
    # Optionally take an input argument when run from Rscript. This should be
    # an integer giving the index of the shapefile we want to run
    #
    # This can be useful to allow the code to be run in batch on NCI
    # with 1 job per shapefile. 
    #
    input_arguments = commandArgs(trailingOnly=TRUE)
    if(length(input_arguments) != 1){
        print('Problem with input arguments')
        print(input_arguments)
        stop()
    }else{
        source_index = as.numeric(input_arguments)
    }

    # Get a vector with all contours that we want to convert to unit sources
    all_sourcezone_shapefiles = all_sourcezone_shapefiles[source_index]
    all_sourcezone_downdip_shapefiles = all_sourcezone_downdip_shapefiles[source_index]
    sourcezone_rake = sourcezone_rake[source_index]

}

## ---- makeDiscretizedSources ----

# Capture plots that occur as source is made in pdf
pdf('UnitSources.pdf', width=10, height=10)

# Loop over all source contour shapefiles, and make the discretized source zone
discretized_sources = list()
discretized_sources_statistics = list()

print('Making discretized sources ...')
for(source_shapefile_index in 1:length(all_sourcezone_shapefiles)){

    source_shapefile = all_sourcezone_shapefiles[source_shapefile_index]
    source_downdip_lines = all_sourcezone_downdip_shapefiles[source_shapefile_index]
    
    # Extract a name for the source
    sourcename = gsub('.shp', '', basename(source_shapefile))
    
    # Create unit sources for source_shapefile
    discretized_sources[[sourcename]] = 
        discretized_source_from_source_contours(source_shapefile, 
            desired_subfault_length, desired_subfault_width, make_plot=TRUE,
            downdip_lines = source_downdip_lines)

    # Get unit source summary stats
    #discretized_sources_statistics[[sourcename]] = 
    #    #discretized_source_approximate_summary_statistics(
    #    discretized_source_summary_statistics(
    #        discretized_sources[[sourcename]],
    #        default_rake = sourcezone_rake[source_shapefile_index],
    #        make_plot=TRUE)
    
    usg = unit_source_grid_to_SpatialPolygonsDataFrame(discretized_sources[[sourcename]]$unit_source_grid)
    proj4string(usg) = '+init=epsg:4326'
    writeOGR(usg, dsn='unit_source_grid', layer=sourcename, driver='ESRI Shapefile', overwrite=TRUE)
}

saveRDS(discretized_sources, 'all_discretized_sources.RDS')


dev.off() # Save pdf plot

## ---- makeTsunamiSources ----

###############################################################################
#
# Step 2: Make tsunami unit sources
#
###############################################################################

dir.create('Unit_source_data', showWarnings=FALSE)

for(sourcename_index in 1:length(names(discretized_sources))){

    sourcename = names(discretized_sources)[sourcename_index]

    # Get the discretized source
    ds1 = discretized_sources[[sourcename]]

    ## Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)

    # Ensure tsunami extent exactly aligns with a degree
    # (in practice this will help us align pixels with our propagation model)
    tsunami_extent = rbind(floor(source_lonlat_extent[c(1,3)] - c(2,2)), 
                           ceiling(source_lonlat_extent[c(2,4)] + c(2,2)))

    tsunami_surface_points_lonlat = expand.grid(
        seq(tsunami_extent[1,1], tsunami_extent[2,1], by = tsunami_source_cellsize),
        seq(tsunami_extent[1,2], tsunami_extent[2,2], by = tsunami_source_cellsize))

    # If elevation data is provided, lookup the depths at the tsunami surface points.
    if(!is.null(elevation_raster)){

        use_kajiura_filter = TRUE

        # Need to ensure that we look up points at longitudes which are within
        # the raster longitude range
        raster_longitude_midpoint = 0.5 * 
            (extent(elevation_raster)@xmin + extent(elevation_raster)@xmax)
    
        ltspl = length(tsunami_surface_points_lonlat[,1])
        tmp_tsp = adjust_longitude_by_360_deg(tsunami_surface_points_lonlat, 
            matrix(raster_longitude_midpoint, ncol=2, nrow=ltspl))

        # Process in chunks to reduce memory usage
        chunk_inds = floor(seq(1, ltspl + 1, len=10))
        surface_point_ocean_depths = tmp_tsp[,1]*NA
        for(i in 1:(length(chunk_inds)-1)){
            inds = chunk_inds[i]:(chunk_inds[i+1]-1)
            surface_point_ocean_depths[inds] = extract(elevation_raster, tmp_tsp[inds,1:2])
            gc()
        }

        # Convert negative elevation to depth, and ensure a minimum depth of 10m
        # for Kajiura filter 
        surface_point_ocean_depths = pmax(-surface_point_ocean_depths, 10)
        rm(tmp_tsp); gc()

    }else{
        # In this case depths are not provided, and Kajiura filtering is not used
        use_kajiura_filter = FALSE
        surface_point_ocean_depths = NULL

    }

    # Make indices for unit sources in parallel computation.
    # If j varies fastest then the shallow unit sources
    # will be submitted early, which will be efficient if
    # they have more interior points (if the spacing is based on the depth)
    ij = expand.grid(j = 1:ds1$discretized_source_dim[2], 
                     i = 1:ds1$discretized_source_dim[1])

    print('Making tsunami sources in parallel...')

    myrake = sourcezone_rake[sourcename_index]

    gc()
    
    source_output_dir = paste0('Unit_source_data/', sourcename, '/')
    dir.create(source_output_dir, showWarnings=FALSE, recursive=TRUE)

    library(parallel)

    # Function to facilitate running in parallel with mcmapply
    parallel_fun<-function(ind){

        # Make a single tsunami unit source 
        down_dip_index = ij$i[ind]
        along_strike_index = ij$j[ind]

        # Set the sub-unit-source point spacing based on the minimum sourcezone depth
        di = down_dip_index:(down_dip_index+1)
        sj = along_strike_index:(along_strike_index+1)
        depth_range = range(ds1$unit_source_grid[di,3,sj])*1000
        approx_dx = min(
            max(shallow_subunitsource_point_spacing, min(depth_range)), 
            deep_subunitsource_point_spacing)
        approx_dy = approx_dx

        # Use within-pixel integration for Okada along the top-row of unit-sources
        local_cell_integration_scale = cell_integration_scale * (down_dip_index == 1)
      
        tsunami_ = make_tsunami_unit_source(
            down_dip_index, 
            along_strike_index, 
            discrete_source=ds1, 
            rake=myrake,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            approx_dx = approx_dx, 
            approx_dy = approx_dy, 
            depths_in_km=TRUE, 
            kajiura_smooth=use_kajiura_filter, 
            surface_point_ocean_depths=surface_point_ocean_depths,
            kajiura_grid_spacing=kajiura_grid_spacing, 
            kajiura_where_deformation_exceeds_threshold=kajiura_use_threshold,
            minimal_output=minimise_tsunami_unit_source_output, 
            verbose=FALSE,
            dstmx=okada_distance_factor,
            edge_taper_width=slip_edge_taper_width,
            cell_integration_scale=local_cell_integration_scale)

        # Save as RDS 
        output_RDS_file =  paste0(source_output_dir, sourcename, '_', 
            down_dip_index, '_', along_strike_index, '.RDS')
        saveRDS(tsunami_, file = output_RDS_file)

        tsunami_source_raster_filename = paste0(source_output_dir, sourcename, '_', 
            down_dip_index, '_', along_strike_index, '.tif')

        # Make a raster
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(tsunami_source_cellsize, tsunami_source_cellsize))

        gc()

        return(output_RDS_file)
    }


    if(MC_CORES > 1){
        all_tsunami_files = mcmapply(parallel_fun, ind=as.list(1:length(ij[,1])),
            mc.cores=MC_CORES, mc.preschedule=TRUE, SIMPLIFY=FALSE)
    }else{
        all_tsunami_files = mapply(parallel_fun, ind=as.list(1:length(ij[,1])), 
            SIMPLIFY=FALSE)
    }


    if(make_pdf_plot){
        # Finally -- read in all the results and make some plots
        all_tsunami = lapply(as.list(all_tsunami_files), f<-function(x) readRDS(x))
       
        all_rasters = paste0(source_output_dir, '/',
            gsub('.RDS', '', basename(unlist(all_tsunami_files))), '.tif')

        all_tsunami_rast = lapply(as.list(all_rasters), f<-function(x) raster(x))

        # Plotting -- make a pdf for checking the sources
        if(make_pdf_plot) plot_all_tsunami_unit_sources(sourcename, all_tsunami, all_tsunami_rast, ds1)

        rm(all_tsunami, all_tsunami_rast); gc()
    }
}

###############################################################################
#
# Optional plotting (interactive)
#
###############################################################################

scatter3d<-function(x, y, z, add=FALSE, ...){
    library(rgl)
    colfun = colorRamp(rainbow(255))
    
    col_01 = (z - min(z))/(max(z) - min(z)+1.0e-20)
    colz = colfun(col_01)
    colz = rgb(colz[,1],colz[,2], colz[,3], maxColorValue=255)
    plot3d(x, y, z, col = colz, add=add, ...)
}

if(make_3d_interactive_plot){
    # NOTE: The next line will need to be changed interactively
    sourcename = 'alaska'
    all_tsunami = lapply(
        Sys.glob(paste0('Unit_source_data/', sourcename , '/', sourcename, '*.RDS')), 
        readRDS)

    print('Computing unit sources for plotting in parallel...')

    ds1 = discretized_sources[[sourcename]]
    origin = ds1$unit_source_grid[1,1:2,1]
    source_lonlat_extent = extent(ds1$depth_contours)

    ## Get surface points for tsunami source
    tsunami_extent = rbind(floor(source_lonlat_extent[c(1,3)] - c(2,2)), 
                           ceiling(source_lonlat_extent[c(2,4)] + c(2,2)))

    tsunami_surface_points_lonlat = expand.grid(
        seq(tsunami_extent[1,1], tsunami_extent[2,1], by = tsunami_source_cellsize),
        seq(tsunami_extent[1,2], tsunami_extent[2,2], by = tsunami_source_cellsize))


    ## Compute interior points for all unit sources for plotting purposes
    unit_source_indices = expand.grid(1:ds1$discretized_source_dim[1], 
        1:ds1$discretized_source_dim[2])
    unit_source_index_list = list()
    for(i in 1:length(unit_source_indices[,1])){
        unit_source_index_list[[i]] = c(unit_source_indices[i,1], 
            unit_source_indices[i,2])
    }

    library(parallel)


    # Make extra unit source points 'for show'
    us = mcmapply(unit_source_interior_points_cartesian,
        discretized_source=list(ds1), 
        unit_source_index = unit_source_index_list,
        origin=list(origin),
        approx_dx = list(5000), approx_dy = list(5000), 
        mc.preschedule=TRUE, mc.cores=MC_CORES, SIMPLIFY=FALSE)

    ## Make a 3D plot of the points inside the unit source
    #for(i in 1:length(us)){
    #    plot3d_unit_source_interior_points_cartesian(us[[i]], add=(i>1))
    #}

    ## Make points of tsunami source FOR PLOTTING.
    ## Origin is the same as unit sources above
    tsunami_source_points_4plot = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat = origin)

    ## Combine all unit sources
    zstore = all_tsunami[[1]]$smooth_tsunami_displacement*0
    for(i in 1:length(all_tsunami)){
        zstore = zstore + all_tsunami[[i]]$smooth_tsunami_displacement
    }

    # Make a 3D plot of the points inside the unit source

    for(i in 1:length(us)){
        plot3d_unit_source_interior_points_cartesian(us[[i]], add=(i>1), 
            add_zero_plane=FALSE)
    }

    #ti = 1

    scatter3d(tsunami_source_points_4plot[,1], tsunami_source_points_4plot[,2],
        zstore*1.0e+05, add=TRUE, size=7)
}
