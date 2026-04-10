#' Get elevation gradients and elevation at the output_grid points.
#'
#' Apply conversions so that the elevation is in m, negative below MSL.
#' The output gradients are always 'change in elevation in m'/'distance in m'
#'
#' @param elevation_raster_file Elevation raster used to compute bathymetric gradients for vertical motion
#' due to horizontal components.
#' @param elevation_raster_scale We need (elevation_rater * elevation_raster_scale) to give the elevation in
#' units of meters above sea level (so bathymetry is negative). 
#' @param isLonLat_elevation_raster If isLonLat_elevation_raster then assume the elevation_raster_file has horizontal
#' units in degrees. Otherwise assume they are in meters.
#'
get_elevation_and_gradients_at_output_grid<-function(output_grid, elevation_raster_file, elevation_raster_scale, isLonLat_elevation_raster,
    output_proj4_string='epsg:4326'){

    elev_rast = raster(elevation_raster_file)
    dx = abs(res(elev_rast)) # abs( [dx, dy] ) or abs([dlon, dlat])

    #
    # North component of slope by central differences
    #
    elevation_output_grid_plus_N  = extract(elev_rast, cbind(output_grid[,1], output_grid[,2] + dx[2]), method='bilinear')
    elevation_output_grid_minus_N = extract(elev_rast, cbind(output_grid[,1], output_grid[,2] - dx[2]), method='bilinear')

    # Distance between 'plus_N' and 'minus_N' locations
    if(isLonLat_elevation_raster){
        print('Assuming elevation raster horizontal units are in degrees')
        d_distance = distHaversine(cbind(output_grid[,1], output_grid[,2] + dx[2]), 
                                   cbind(output_grid[,1], output_grid[,2] - dx[2]))
    }else{
        print('Assuming elevation raster horizontal units are in meters')
        d_distance = 2*dx[2]
    }
    output_bathy_slope_N = (elevation_output_grid_plus_N - elevation_output_grid_minus_N)*elevation_raster_scale/d_distance
    rm(elevation_output_grid_plus_N, elevation_output_grid_minus_N, d_distance); gc()

    #
    # East component of slope by central differences
    #
    elevation_output_grid_plus_E = extract(elev_rast,  cbind(output_grid[,1]+dx[1], output_grid[,2]), method='bilinear')
    elevation_output_grid_minus_E = extract(elev_rast, cbind(output_grid[,1]-dx[1], output_grid[,2]), method='bilinear')
    if(isLonLat_elevation_raster){
        d_distance = distHaversine(cbind(output_grid[,1]+dx[1], output_grid[,2]), 
                                   cbind(output_grid[,1]-dx[1], output_grid[,2]))
    }else{
        d_distance = 2*dx[1]
    }
    output_bathy_slope_E = (elevation_output_grid_plus_E - elevation_output_grid_minus_E)*elevation_raster_scale/d_distance
    rm(elevation_output_grid_plus_E, elevation_output_grid_minus_E, d_distance); gc()

    slope_N = rasterFromXYZ(
        xyz=cbind(output_grid[,1], output_grid[,2], output_bathy_slope_N),
        crs=CRS(output_proj4_string))
    slope_E = rasterFromXYZ(
        xyz=cbind(output_grid[,1], output_grid[,2], output_bathy_slope_E),
        crs=CRS(output_proj4_string))

    # Also extract the elevation for QC
    elevation = extract(elev_rast, cbind(output_grid[,1], output_grid[,2]), method='bilinear') * elevation_raster_scale
    elev_m_above_sea_level = rasterFromXYZ(
        xyz=cbind(output_grid[,1], output_grid[,2], elevation),
        crs=CRS(output_proj4_string))

    return(list(slope_N = slope_N, slope_E = slope_E, elev_m_above_sea_level = elev_m_above_sea_level))
}
