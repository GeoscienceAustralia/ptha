# Main 'driver' script to create the unit sources
#
# Gareth Davies, Geoscience Australia 2015
#
library(rptha)

# Main input parameters 
desired_subfault_length = 100 # km
desired_subfault_width = 50 # km
MC_CORES = 12 # Number of cores for parallel parts
tsunami_source_cellsize = 1/60 # degrees

# Get a vector with all contours that we want to convert to unit sources
all_sourcezone_shapefiles = Sys.glob('./CONTOURS/*.shp')

# Only make the 3d interactive plot if you can use interactive graphics and
# have rgl (i.e. use FALSE on NCI). 
make_3d_interactive_plot = FALSE 

###############################################################################
#
# Step 1: Make discretized source zones for all ruptures
#
###############################################################################

# Capture plots that occur as source is made in pdf
pdf('UnitSources.pdf', width=10, height=10)

# Loop over all source contour shapefiles, and make the discretized source zone
discretized_sources = list()
discretized_sources_statistics = list()

print('Making discretized sources ...')
for(interface_shapefile in all_sourcezone_shapefiles){

    # Extract a name for the source
    sourcename = gsub('.shp', '', basename(interface_shapefile))
    
    # Create discretized sources for interface_shapefile
    discretized_sources[[sourcename]] = 
        discretized_source_from_source_contours(interface_shapefile, 
            desired_subfault_length, desired_subfault_width, make_plot=TRUE)

    # Get discretized source approximate summary stats
    # More accurate results can be obtained from a similar function without
    # _approximate_ in the name, but this is a bit more intensive
    discretized_sources_statistics[[sourcename]] = 
        discretized_source_approximate_summary_statistics(
            discretized_sources[[sourcename]],
            make_plot=TRUE)
}


dev.off() # Save pdf plot

# Save the R image (can be useful for testing/debugging)
dir.create('Rimages', showWarnings=FALSE)
save.image('Rimages/post_sources.Rdata')

###############################################################################
#
# Step 2: Make tsunami unit sources
#
###############################################################################

dir.create('Unit_source_data', showWarnings=FALSE)

for(sourcename in names(discretized_sources)){

    # Get the discretized source
    ds1 = discretized_sources[[sourcename]]

    ## Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)

    # Ensure tsunami extent exactly aligns with a degree
    # (in practice this will help us align pixels with our propagation model)
    tsunami_extent = rbind(floor(source_lonlat_extent[c(1,3)] - c(2,2)), 
                           ceiling(source_lonlat_extent[c(2,4)] + c(2,2)))

    # Create a grid of points where the surface deformation is computed
    tsunami_surface_points_lonlat = expand.grid(
        seq(tsunami_extent[1,1], tsunami_extent[2,1], by = tsunami_source_cellsize),
        seq(tsunami_extent[1,2], tsunami_extent[2,2], by = tsunami_source_cellsize))

    # Make unit source indices for parallel computation.
    #
    # Each unit source is indexed by its position along-strike, and down-dip
    # Typically shallow unit sources require more integration points, so
    # are computationally expensive. Thus, in a parallel computation, we want
    # to start those ones first.
    #
    # If j varies fastest then the shallow unit sources will be submitted early
    # So make j the first argument in expand.grid
    ij = expand.grid(j = 1:ds1$discretized_source_dim[2], 
                     i = 1:ds1$discretized_source_dim[1])

    # Approximate spacing for sub-unit-source integration points
    approx_dx = (ij$i > 1)*5000 + 1000
    approx_dy = approx_dx

    print('Making tsunami sources in parallel...')

    # i and j vary, other arguments are recycled
    library(parallel)
    all_tsunami = mcmapply(
        make_tsunami_unit_source, 
        i = as.list(ij$i), 
        j = as.list(ij$j), 
        discrete_source=list(ds1), 
        rake=list(90), # Pure thrust
        tsunami_surface_points_lonlat = list(tsunami_surface_points_lonlat),
        approx_dx = as.list(approx_dx), 
        approx_dy = as.list(approx_dy), 
        depths_in_km=list(TRUE),
        mc.cores=MC_CORES, # Parallel arguments here and below
        mc.preschedule=FALSE, 
        SIMPLIFY=FALSE)

    print('Saving outputs...')
    # Save results to a file
    saveRDS(all_tsunami,
        file = paste0('Unit_source_data/', sourcename, '.RDS'))

    # Make output rasters
    #  Construct filename like 'OutputDir/source_dipIndex_strikeIndex.tif'
    raster_dir = paste0('Unit_source_data/', sourcename, '/') 
    down_dip_index = unlist(lapply(all_tsunami, f<-function(x) x$i))
    along_strike_index = unlist(lapply(all_tsunami, f<-function(x) x$j))
    dir.create(raster_dir, showWarnings=FALSE)
    tsunami_source_raster_filenames = paste0(raster_dir, sourcename, '_', 
        down_dip_index, '_', along_strike_index, '.tif')

    # Make the rasters in parallel
    all_tsunami_rast = mcmapply(tsunami_unit_source_2_raster, 
        all_tsunami, filename = as.list(tsunami_source_raster_filenames),
        saveonly = list(FALSE),
        mc.cores=MC_CORES, mc.preschedule=FALSE, SIMPLIFY=FALSE)

    ## Plotting -- make a pdf for checking the sources
    plot_all_tsunami_unit_sources(sourcename, all_tsunami, all_tsunami_rast, ds1)
}

###############################################################################
#
# Optional plotting (interactive)
#
###############################################################################

# Convenience plotting
scatter3d<-function(x, y, z, add=FALSE, ...){
    library(rgl)
    colfun = colorRamp(rainbow(255))
    
    col_01 = (z - min(z))/(max(z) - min(z)+1.0e-20)
    colz = colfun(col_01)
    colz = rgb(colz[,1],colz[,2], colz[,3], maxColorValue=255)
    plot3d(x, y, z, col = colz, add=add, ...)
}

if(make_3d_interactive_plot){

    sourcename = 'alaska'
    all_tsunami = readRDS(paste0('Unit_source_data/', sourcename, '.RDS'))

    print('Computing unit sources for plotting in parallel...')

    ds1 = discretized_sources[[sourcename]]
    origin = ds1$unit_source_grid[1,1:2,1]

    ## Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)
    tsunami_surface_points_lonlat = all_tsunami[[1]]$tsunami_surface_points_lonlat


    ## Compute interior points for all unit sources for plotting purposes

    unit_source_indices = expand.grid(1:ds1$discretized_source_dim[1], 
        1:ds1$discretized_source_dim[2])
    unit_source_index_list = list()
    for(i in 1:length(unit_source_indices[,1])){
        unit_source_index_list[[i]] = c(unit_source_indices[i,1], 
            unit_source_indices[i,2])
    }

    library(parallel)

    us = mcmapply(unit_source_interior_points_cartesian,
        discretized_source=list(ds1), 
        unit_source_index = unit_source_index_list,
        origin=list(origin),
        approx_dx = list(NULL), approx_dy = list(NULL), 
        mc.preschedule=FALSE, mc.cores=MC_CORES, SIMPLIFY=FALSE)

    ## Make points of tsunami source FOR PLOTTING.
    ## Origin is the same as unit sources above
    tsunami_source_points_4plot = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat = origin)
    ## Combine all unit sources
    zstore = all_tsunami[[1]]$tsunami_source$zdsp*0
    for(i in 1:length(all_tsunami)){
        zstore = zstore + all_tsunami[[i]]$tsunami_source$zdsp
    }

    # Make a 3D plot of the points inside the unit source

    for(i in 1:length(us)){
        plot3d_unit_source_interior_points_cartesian(us[[i]], add=(i>1), 
            add_zero_plane=FALSE)
    }

    scatter3d(tsunami_source_points_4plot[,1], tsunami_source_points_4plot[,2],
        #all_tsunami[[ti]]$ts$zdsp*1.0e+05, add=TRUE, size=7)
        zstore*1.0e+05, add=TRUE, size=7)
}
