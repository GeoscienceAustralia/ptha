## ---- initialise ----

# Main 'driver' script to create the unit sources
#
# Gareth Davies, Geoscience Australia 2015
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
all_sourcezone_shapefiles = Sys.glob('./CONTOURS/*.shp') # Matches all shapefiles in CONTOURS

# Desired unit source geometric parameters
desired_subfault_length = 100 # km
desired_subfault_width = 50 # km

# A vector with the desired rake angle (one entry per sourcezone)
sourcezone_rake = rep(90, len=length(all_sourcezone_shapefiles)) # degrees

# Desired spacing of sub-unit-source points
# Lower values (e.g. 1000) may be required for accuracy in unit sources
# along the trench, because shallow deformation tends to be quite localised.
# For deeper unit sources, a much coarser point spacing can be used without
# sacrificing accuracy. 
# Hence we use different values for the 'shallow' sub-unit-source points (i.e.
# along the trench) and the deeper ones.
# The computational effort approximately scales with the inverse square of
# the point density. 
shallow_subunitsource_point_spacing = 1000 # m
deep_subunitsource_point_spacing = 6000 #m

# For computational efficiency, only compute the okada deformation at
# distances <= okada_distance_factor x (depth of sub-unit-source point) 
# This can save computational effort for shallow unit sources.
# But be careful if using a wide subunitsource_point_spacing.
okada_distance_factor = 50 # Inf 

# elevation raster (required for Kajiura filtering). Should give elevation in m, 
# with the ocean having elevation < 0. Should have a lon/lat spatial projection. 
# Set to NULL to not use Kajiura filtering
elevation_raster = NULL 
# elevation_raster = raster('../RAW/GEBCO/gebco_08.nc')

# For computational efficiency, only apply Kajiura filtering in a box
# containing all points where the unit source deformation exceeds
# kajiura_use_threshold. Set to zero to apply Kajiura filter everywhere.
# Use of a small positive number can be faster.
kajiura_use_threshold = 1.0e-04

# When applying the kajiura filter, the data is regridded onto a grid with
# spacing=kajiura_gridspacing. The latter should be small compared to the
# horizontal distance over which the deformation changes significantly
kajiura_grid_spacing = 1000 # m

# Cell size for output rasters
# The computation time will scale inversely with this squared
# Here we use a relatively coarse discretization, for demonstration purposes
tsunami_source_cellsize = 4/60 # degrees. 

# Number of cores for parallel parts. Values > 1 will only work on shared
# memory linux machines.
MC_CORES = 12 

# Option to illustrate 3d interactive plot creation
#
# Only make the 3d interactive plot if you can use interactive graphics and
# have rgl (i.e. use FALSE on NCI). 
make_3d_interactive_plot = FALSE 

# Option to reduce the size of RDS output
# Set to FALSE if make_3d_interactive_plot = TRUE, or if you have other
# needs to use the detailed outputs.
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
    
    # Extract a name for the source
    sourcename = gsub('.shp', '', basename(source_shapefile))
    
    # Create unit sources for source_shapefile
    discretized_sources[[sourcename]] = 
        discretized_source_from_source_contours(source_shapefile, 
            desired_subfault_length, desired_subfault_width, make_plot=TRUE)

    # Get unit source summary stats
    discretized_sources_statistics[[sourcename]] = 
        discretized_source_approximate_summary_statistics(
            discretized_sources[[sourcename]],
            default_rake = sourcezone_rake[source_shapefile_index],
            make_plot=TRUE)
}


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

        # FIXME:  This assumes the elevation raster ranges from [-180 , 180]
        #
        # Make some temporary tsunami surface points with longitudes in
        # [-180,180] so we can lookup the elevation raster
        ltspl = length(tsunami_surface_points_lonlat[,1])
        tmp_tsp = adjust_longitude_by_360_deg(tsunami_surface_points_lonlat, 
            matrix(0, ncol=2, nrow=ltspl))

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

    approx_dx = (ij$i > 1)*deep_subunitsource_point_spacing + 
        (ij$i == 1)*shallow_subunitsource_point_spacing
    approx_dy = approx_dx 

    myrake = sourcezone_rake[sourcename_index]

    gc()
    
    raster_dir = paste0('Unit_source_data/', sourcename, '/')
    dir.create(raster_dir, showWarnings=FALSE, recursive=TRUE)

    library(parallel)

    # Function to facilitate running in parallel with mcmapply
    parallel_fun<-function(ind){
        # Make a single tsunami unit source 
        down_dip_index = ij$i[ind]
        along_strike_index = ij$j[ind]
      
        tsunami_ = make_tsunami_unit_source(
            down_dip_index, 
            along_strike_index, 
            discrete_source=ds1, 
            rake=myrake,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            approx_dx = approx_dx[ind], 
            approx_dy = approx_dy[ind], 
            depths_in_km=TRUE, 
            kajiura_smooth=use_kajiura_filter, 
            surface_point_ocean_depths=surface_point_ocean_depths,
            kajiura_grid_spacing=kajiura_grid_spacing, 
            kajiura_where_deformation_exceeds_threshold=kajiura_use_threshold,
            minimal_output=minimise_tsunami_unit_source_output, 
            dstmx=okada_distance_factor)

        # Save as RDS 
        output_RDS_file =  paste0('Unit_source_data/', sourcename, '_', 
            down_dip_index, '_', along_strike_index, '.RDS')
        saveRDS(tsunami_, file = output_RDS_file)

        tsunami_source_raster_filename = paste0(raster_dir, sourcename, '_', 
            down_dip_index, '_', along_strike_index, '.tif')

        # Make a raster
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(tsunami_source_cellsize, tsunami_source_cellsize))

        return(output_RDS_file)
    }

    if(MC_CORES > 1){
        all_tsunami_files = mcmapply(parallel_fun, ind=as.list(1:length(ij[,1])),
            mc.cores=MC_CORES, mc.preschedule=FALSE, SIMPLIFY=FALSE)
    }else{
        all_tsunami_files = mapply(parallel_fun, ind=as.list(1:length(ij[,1])), 
            SIMPLIFY=FALSE)
    }


    # Finally -- read in all the results and make some plots
    all_tsunami = lapply(as.list(all_tsunami_files), f<-function(x) readRDS(x))
   
    all_rasters = paste0(raster_dir, '/',
        gsub('.RDS', '', basename(unlist(all_tsunami_files))), '.tif')

    all_tsunami_rast = lapply(as.list(all_rasters), f<-function(x) raster(x))

    # Plotting -- make a pdf for checking the sources
    plot_all_tsunami_unit_sources(sourcename, all_tsunami, all_tsunami_rast, ds1)

    rm(all_tsunami, all_tsunami_rast); gc()
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
        Sys.glob(paste0('Unit_source_data/', sourcename, '*.RDS')), 
        readRDS)

    print('Computing unit sources for plotting in parallel...')

    ds1 = discretized_sources[[sourcename]]
    origin = ds1$unit_source_grid[1,1:2,1]

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
        mc.preschedule=FALSE, mc.cores=MC_CORES, SIMPLIFY=FALSE)

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
