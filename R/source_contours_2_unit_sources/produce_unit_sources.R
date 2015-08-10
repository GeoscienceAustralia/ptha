# Main 'driver' script to create the unit sources
#
# Gareth Davies, Geoscience Australia 2015
#

# Get main functions
#source('unit_sources.R')
#source('tsunami_sources.R')
library(rptha)

###############################################################################
#
# Step 1: Make discretized source zones for all ruptures
#
###############################################################################

desired_subfault_length = 100 # km
desired_subfault_width = 50 # km
MC_CORES = 12 # Number of cores for parallel parts
tsunami_source_pixels = 900

# Get a vector with all contours that we want to convert to unit sources
all_sourcezone_shapefiles = 
    Sys.glob('../../MAKE_SOURCE_CONTOURS/OUTPUT_DATA/CONTOURS_FULL/*.shp')[1]

# Capture plots that occur as source is made in pdf
pdf('UnitSources.pdf', width=10, height=10)

# Loop over all source contour shapefiles, and make the discretized source zone
discretized_sources = list()
discretized_sources_statistics = list()

print('Making discretized sources ...')
for(interface_shapefile in all_sourcezone_shapefiles){

    # Extract a name for the source
    sourcename = gsub('.shp', '', basename(interface_shapefile))
    
    # Create unit sources for interface_shapefile
    discretized_sources[[sourcename]] = 
        discretized_source_from_source_contours(interface_shapefile, 
            desired_subfault_length, desired_subfault_width, make_plot=TRUE)

    # Get unit source summary stats
    discretized_sources_statistics[[sourcename]] = 
        discretized_source_summary_statistics(discretized_sources[[sourcename]],
            make_plot=TRUE)
}


dev.off() # Save pdf plot

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
    tsunami_surface_points_lonlat = expand.grid(
        seq(source_lonlat_extent[1]-2, source_lonlat_extent[2]+2, len=tsunami_source_pixels),
        seq(source_lonlat_extent[3]-2, source_lonlat_extent[4]+2, len=tsunami_source_pixels))

    # Make indices for unit sources in parallel computation.
    # If j varies fastest then the shallow unit sources
    # will be submitted early, which will be efficient if
    # they have more interior points (if the spacing is based on the depth)
    ij = expand.grid(j = 1:ds1$discretized_source_dim[2], 
                     i = 1:ds1$discretized_source_dim[1])

    print('Making tsunami sources in parallel...')

    # i and j vary, other arguments are recycled
    library(parallel)
    all_tsunami = mcmapply(
        make_tsunami_unit_source, 
        i = as.list(ij$i), j = as.list(ij$j), discrete_source=list(ds1), 
        tsunami_surface_points_lonlat = list(tsunami_surface_points_lonlat),
        approx_dx = list(NULL), approx_dy = list(NULL), scale_dxdy=list(1.0), 
        depths_in_km=list(TRUE),
        mc.cores=MC_CORES, mc.preschedule=FALSE, SIMPLIFY=FALSE)

    # Save results to a file
    saveRDS(all_tsunami,
        file = paste0('Unit_source_data/', sourcename, '.RDS'))

    ## Make output rasters
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

#' Convenience plotting
#'
#' @export
scatter3d<-function(x, y, z, colramp = 'cpt-city/ds9/rainbow.cpt', add=FALSE, 
    ...){
    library(rgl)
    library(colorRampPC)

    colfun = colorRampPC(colramp)

    colz = colfun( (z - min(z))/diff(range(z)))
    plot3d(x, y, z, col = colz, add=add, ...)
}

make_plot = FALSE
if(make_plot){

    sourcename = 'aleutians'
    all_tsunami = readRDS(paste0('Unit_source_data/', sourcename, '.RDS'))

    print('Computing unit sources for plotting in parallel...')

    ds1 = discretized_sources[[sourcename]]
    origin = ds1$unit_source_grid[1,1:2,1]

    ## Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)
    tsunami_surface_points_lonlat = expand.grid(
        seq(source_lonlat_extent[1]-2, source_lonlat_extent[2]+2, len=900),
        seq(source_lonlat_extent[3]-2, source_lonlat_extent[4]+2, len=900))


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

    ## Make a 3D plot of the points inside the unit source
    #for(i in 1:length(us)){
    #    plot3d_unit_source_interior_points_cartesian(us[[i]], add=(i>1))
    #}

    ## Make points of tsunami source FOR PLOTTING.
    ## Origin is the same as unit sources above
    tsunami_source_points_4plot = spherical_to_cartesian2d_coordinates(
        tsunami_surface_points_lonlat, origin_lonlat = origin)
    ## Combine all unit sources
    zstore = all_tsunami[[1]]$ts$zdsp*0
    for(i in 1:length(all_tsunami)){
        zstore = zstore + all_tsunami[[i]]$ts$zdsp
    }

    # Make a 3D plot of the points inside the unit source

    for(i in 1:length(us)){
        plot3d_unit_source_interior_points_cartesian(us[[i]], add=(i>1), 
            add_zero_plane=FALSE)
    }

    #ti = 1

    scatter3d(tsunami_source_points_4plot[,1], tsunami_source_points_4plot[,2],
        #all_tsunami[[ti]]$ts$zdsp*1.0e+05, add=TRUE, size=7)
        zstore*1.0e+05, add=TRUE, size=7)
}
