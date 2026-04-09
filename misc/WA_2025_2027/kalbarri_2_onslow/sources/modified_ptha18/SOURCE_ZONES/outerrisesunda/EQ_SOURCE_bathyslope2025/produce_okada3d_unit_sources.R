## ---- initialise ----

# Main 'driver' script to create the unit sources
#
# The code was modified from the PTHA18 scripts which make tsunami initial
# conditions - because the latter only consider the vertical component
#
library(rptha)
library(parallel)

###############################################################################
#
# Main input parameters -- e.g. giving the sourcezone name
#
###############################################################################

source('config_okada3d.R')
source('get_elevation_and_gradients_at_output_grid.R')


## ---- makeDiscretizedSources ----

# Loop over all source contour shapefiles, and make the discretized source zone
discretized_sources = list()

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
    proj4string(usg) = 'epsg:4326'
    writeOGR(usg, dsn=paste0(output_base_dir, 'unit_source_grid'),
        layer=sourcename, driver='ESRI Shapefile', overwrite=TRUE)
}

saveRDS(discretized_sources, paste0(output_base_dir, 'all_discretized_sources.RDS'))


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

    # Get surface points for tsunami source
    source_lonlat_extent = extent(ds1$depth_contours)

    # Ensure tsunami extent exactly aligns with a degree
    # (in practice this will help us align pixels with our propagation model)
    output_raster_extent = rbind(
        floor(  source_lonlat_extent[c(1,3)] - c(1,1)*output_raster_ext_buffer), 
        ceiling(source_lonlat_extent[c(2,4)] + c(1,1)*output_raster_ext_buffer))
    tsunami_extent = matrix(output_raster_extent, ncol=2)
    tsunami_surface_points_lonlat = expand.grid(
        seq(tsunami_extent[1,1], tsunami_extent[2,1], by = output_raster_cellsize),
        seq(tsunami_extent[1,2], tsunami_extent[2,2], by = output_raster_cellsize))

    # Get bathymetric gradients
    elev_slope = get_elevation_and_gradients_at_output_grid(tsunami_surface_points_lonlat, 
        elevation_raster_file, elevation_raster_scale, isLonLat_elevation_raster); gc()
    # Save the elevation and its gradients as a tif for QC
    for(nm in names(elev_slope)){
        writeRaster(elev_slope[[nm]], 
            file=paste0(output_base_dir, '/elevation_', nm, '.tif'), 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }

    # Make indices for unit sources in parallel computation.
    # If j varies fastest then the shallow unit sources
    # will be submitted early, which will be efficient if
    # they have more interior points (if the spacing is based on the depth)
    ij = expand.grid(j = 1:ds1$discretized_source_dim[2], 
                     i = 1:ds1$discretized_source_dim[1])

    myrake = sourcezone_rake[sourcename_index]

    gc()
    
    source_output_dir = paste0(output_base_dir, 'Unit_source_data/', sourcename, '/')
    dir.create(source_output_dir, showWarnings=FALSE, recursive=TRUE)

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
        
        # Coordinate near-ish to middle of unit source for spherical-to-cartesian
        # transformation inside kajiura.
        lonlat_origin_for_kajiura = c(mean(ds1$unit_source_grid[di,1,sj]), mean(ds1$unit_source_grid[di,2,sj]))

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
            kajiura_smooth=FALSE, 
            surface_point_ocean_depths=NULL,
            kajiura_grid_spacing=NA, 
            kajiura_where_deformation_exceeds_threshold=NA,
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
        tsunami_source_raster_filename_z = paste0(source_output_dir, sourcename, '_vertical_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename_z, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='vertical_displacement')

        # Store easting displacement as raster
        tsunami_source_raster_filename_e = paste0(source_output_dir, sourcename, '_easting_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename_e, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='easting_displacement')

        # Store northing displacement as raster
        tsunami_source_raster_filename_n = paste0(source_output_dir, sourcename, '_northing_displacement_', 
            down_dip_index, '_', along_strike_index, '.tif')
        tsunami_unit_source_2_raster(
            tsunami_, tsunami_source_raster_filename_n, saveonly=TRUE,
            tsunami_surface_points_lonlat = tsunami_surface_points_lonlat,
            res=c(output_raster_cellsize, output_raster_cellsize),
            field_to_output='northing_displacement')

        # Compute vertical displacement due to both uplift, and horizontal components along depth gradients
        vertical_displacement_unit_source = raster(tsunami_source_raster_filename_z) +
            (-1.0)*elev_slope$slope_E*raster(tsunami_source_raster_filename_e) + 
            (-1.0)*elev_slope$slope_N*raster(tsunami_source_raster_filename_n) 
        # Save it
        tsunami_source_raster_filename_vert = paste0(source_output_dir, sourcename, '_vertical_inc_bathyslope_',
            down_dip_index, '_', along_strike_index, '.tif')
        writeRaster(vertical_displacement_unit_source, 
            file=tsunami_source_raster_filename_vert, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

        # Apply a kajiura filter to the combined vertical displacement
        smoothed = kajiura_smooth_raster(vertical_displacement_unit_source, 
            new_origin = lonlat_origin_for_kajiura, 
            elevation_raster = elev_slope$elev_m_above_sea_level, 
            spherical_input=TRUE,
            # The kajiura_filter throws an error if the volume change caused by KJ filter is too great.
            # Allow for errors slightly larger than default.
            volume_change_error_threshold = 0.03)
        # Save it
        tsunami_source_raster_filename_vert_ks = paste0(source_output_dir, sourcename, '_vertical_inc_bathyslope_kajiura_',
            down_dip_index, '_', along_strike_index, '.tif')
        writeRaster(smoothed, 
            file=tsunami_source_raster_filename_vert_ks, 
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)


        return(output_RDS_file)
    }
          
    # Run the code in parallel if there are more than 1 cores available
    if(MC_CORES == 1){
        all_tsunami_files = mapply(parallel_fun, ind=as.list(1:length(ij[,1])), 
            SIMPLIFY=FALSE)
    }else{
        print('Making tsunami sources in parallel...')

        # Beware mcmapply(...) is extremely slow on NCI although it works fine on my home machine.
        # Seems efficient to setup a proper cluster
        cl = makeCluster(MC_CORES)
        nothing = clusterCall(cl, function(){ library(rptha); library(parallel) })  
        nothing = clusterExport(cl, varlist=ls())
        all_tsunami_files = clusterApplyLB(cl, as.list(1:length(ij[,1])), parallel_fun)
        stopCluster(cl)
    }
}

