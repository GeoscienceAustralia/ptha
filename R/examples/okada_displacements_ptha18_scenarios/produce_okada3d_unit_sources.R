## ---- initialise ----

# Main 'driver' script to create the unit sources
#
# Gareth Davies, Geoscience Australia 2019.
#
# The code was modified from the PTHA18 scripts which make tsunami initial
# conditions - because the latter only consider the vertical component, and
# apply a Kajiura filter - whereas here we save the unfiltered
# east/north/vertical displacements.
#
library(rptha)
library(raster)

###############################################################################
#
# Main input parameters -- e.g. giving the sourcezone name
#
###############################################################################

source('config.R')


## ---- makeDiscretizedSources ----

# Capture plots that occur as source is made in pdf
#pdf('UnitSources.pdf', width=10, height=10)

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
            desired_subfault_length, desired_subfault_width, make_plot=FALSE,
            downdip_lines = source_downdip_lines)

    # Get unit source summary stats
    usg = unit_source_grid_to_SpatialPolygonsDataFrame(
        discretized_sources[[sourcename]]$unit_source_grid)
    proj4string(usg) = '+init=epsg:4326'
    writeOGR(usg, dsn=paste0(output_base_dir, 'unit_source_grid'),
        layer=sourcename, driver='ESRI Shapefile', overwrite=TRUE)
}

saveRDS(discretized_sources, paste0(output_base_dir, 'all_discretized_sources.RDS'))


#dev.off() # Save pdf plot

## ---- makeTsunamiSources ----

###############################################################################
#
# Step 2: Make okada unit sources
#
###############################################################################

for(sourcename_index in 1:length(names(discretized_sources))){

    sourcename = names(discretized_sources)[sourcename_index]

    # Get the discretized source
    ds1 = discretized_sources[[sourcename]]

    ## Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)
    tsunami_extent = matrix(output_raster_extent, ncol=2)
    tsunami_surface_points_lonlat = expand.grid(
        seq(tsunami_extent[1,1], tsunami_extent[2,1], by = output_raster_cellsize),
        seq(tsunami_extent[1,2], tsunami_extent[2,2], by = output_raster_cellsize))

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
            minimal_output=FALSE, # If TRUE then Okada vector is not stored! 
            verbose=FALSE,
            dstmx=okada_distance_factor,
            edge_taper_width=slip_edge_taper_width,
            cell_integration_scale=local_cell_integration_scale)

        # Save as RDS 
        output_RDS_file =  paste0(source_output_dir, sourcename, '_', 
            down_dip_index, '_', along_strike_index, '.RDS')
        saveRDS(tsunami_, file = output_RDS_file)

        #
        # Raster outputs
        #

        # Store vertical displacement as raster
        tsunami_source_raster_filename = paste0(source_output_dir, sourcename, '_vertical_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='vertical_displacement')

        # Store easting displacement as raster
        tsunami_source_raster_filename = paste0(source_output_dir, sourcename, '_easting_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='easting_displacement')

        # Store northing displacement as raster
        tsunami_source_raster_filename = paste0(source_output_dir, sourcename, '_northing_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='northing_displacement')

        gc()

        return(output_RDS_file)
    }


    if(MC_CORES > 1){
        all_tsunami_files = mcmapply(parallel_fun, ind=as.list(1:length(ij[,1])),
            mc.cores=MC_CORES, mc.preschedule=FALSE, SIMPLIFY=FALSE)
    }else{
        all_tsunami_files = mapply(parallel_fun, ind=as.list(1:length(ij[,1])), 
            SIMPLIFY=FALSE)
    }

}

