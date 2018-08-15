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

source('config.R')


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
    
    usg = unit_source_grid_to_SpatialPolygonsDataFrame(
        discretized_sources[[sourcename]]$unit_source_grid)
    proj4string(usg) = '+init=epsg:4326'
    writeOGR(usg, dsn=paste0(output_base_dir, 'unit_source_grid'),
        layer=sourcename, driver='ESRI Shapefile', overwrite=TRUE)
}

saveRDS(discretized_sources, paste0(output_base_dir, 'all_discretized_sources.RDS'))


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
        # for Kajiura filter. NOTE: When sources are on-land it may be better to increase
        # this 10m limit to avoid running out of memory (because it affects the spacing of points
        # in the kajiura filter). The only time I saw this was the 'makran2' source in PTHA18
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
    
    source_output_dir = paste0(output_base_dir, 'Unit_source_data/', sourcename, '/')
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
    sourcename = site_name #### 'alaska'
    all_tsunami = lapply(
        Sys.glob(paste0(output_base_dir, sourcename, '/Unit_source_data/', 
            sourcename , '/', sourcename, '*.RDS')),
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
