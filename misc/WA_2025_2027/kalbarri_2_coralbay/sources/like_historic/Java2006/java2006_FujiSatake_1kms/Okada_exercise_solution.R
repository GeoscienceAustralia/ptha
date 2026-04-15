#
# Code to compute rupture from multiple earthquake sources, 
# using the 1km/s source model from
# Fujii, Y. & Satake, K. Source of the July 2006 West Java tsunami estimated from
# tide gauge records Geophysical Research Letters, 2006, 33
# 

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

#
# Input parameters -- rupture geometry -- 1km/s source
#
lat = c( -9.515,  -9.096,  -9.661,  -9.242,  -9.807,  -9.388,  -9.953, -9.534, -10.1, -9.681)  # top-east edge
lon = c(107.274, 107.420, 107.705, 107.851, 108.136, 108.282, 108.568, 108.714, 109.000, 109.146) # top-east edge
width = 50.0 + 0*lon # rupture width in km
length= 50.0 + 0*lon # rupture length in km
strike= 289.0 + 0*lon # Strike in degrees clockwise from N
dip = 10.0 + 0*lon 
depth = rep(c(3, 11.7), 5) # Depth at top
rake  = 95 + 0*lon # Degrees -- rupture rake
slip = c(0.50, 0.19, 0.67, 0.30, 1.81, 0.95, 2.42, 2.73, 0.26, 3.69)
rupture_name   = as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'FujiiSatake06_Java2006' # Name at the start of output files

# Adjust lon/lat to be middle of top edge
new_lonlat = destPoint(cbind(lon, lat), b=strike, d=length/2*1000)
lon = new_lonlat[,1] 
lat = new_lonlat[,2]

#
# Flags to denote meaning of lon, lat and depth
#
# Lon/Lat location: 
# Use -1 for 'middle of bottom edge of rupture', 0 for 'centroid', and 1 for
# 'middle of top edge of rupture'
lon_lat_location = 1 
# Depth location:
# Use -1 for 'depth at bottom edge of rupture', 0 for 'depth at centroid', and
# 1 for 'depth at top edge of rupture'.
depth_location = 1 

# Elevation raster used to compute bathymetric gradients for vertical motion
# due to horizontal components.
elevation_raster = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/MODELS/AustPTHA_c/DATA/ELEV/merged_dem/merged_gebco_ga250_dem_patched.tif'
# We need (elevation_rater * elevation_raster_scale) to give the elevation in
# units of meters above sea level (so bathymetry is negative). 
elevation_raster_scale = 1.0 
# If isLonLat_elevation_raster then assume the elevation_raster has horizontal
# units in degrees. Otherwise assume they are in meters.
isLonLat_elevation_raster = isLonLat(raster(elevation_raster))
# Beware isLonLat(..) may not work if the elevation raster doesn't have
# proper projection information. In that case, better to force the answer, e.g.
#   isLonLat_elevation_raster = TRUE


#
# Input parameters -- output grid
#
output_res = 0.5/60 # 30 second resolution
output_lons = seq(103, 113, by=output_res)
output_lats = seq(-14, -3, by=output_res)
output_proj4_string = "epsg:4326" # code for WGS84 


#
# END INPUT PARAMETERS
#


# Make a 'grid' of points for output
output_grid = expand.grid(output_lons, output_lats)

# Offsets (in m east/north) to get from output_grid to centroid
deg2rad = pi/180
dlon = width/2 * cos(dip*deg2rad) * cos(strike*deg2rad) * 1000 * lon_lat_location 
dlat = -width/2 * cos(dip*deg2rad) * sin(strike*deg2rad) * 1000 * lon_lat_location
# Centroid depth (in km)
centroid_depth = depth + width/2 * sin(dip*deg2rad) * depth_location

# Output data in a convenient way
rup2csv<-function(){

    # Estimate the centroid longitude/latitude
    newlon = lon + destPoint(cbind(0, lat), b=90, d=dlon, f=0)[,1]
    newlat = destPoint(cbind(lon, lat), b=0,  d=dlat, f=0)[,2]

    output = data.frame(lonc=newlon, latc=newlat, depthc=centroid_depth,
        length=length, width=width, strike=strike, rake=rake, dip=dip, slip=slip)

    write.csv(output, 'rupture_summary.csv', row.names=FALSE)
}
rup2csv()

get_elevation_and_gradients_at_output_grid<-function(output_grid, elevation_raster, elevation_raster_scale, isLonLat_elevation_raster){
    #
    # Get elevation gradients at the output_grid points.
    # Apply conversions so that the elevation is in m, negative below MSL.
    # The output gradients are always 'change in elevation in m'/'distance in m'
    #
    elev_rast = raster(elevation_raster)
    dx = abs(res(elev_rast)) # abs( [dx, dy] ) or abs([dlon, dlat])

    #
    # North component of slope by central differences
    #
    elevation_output_grid_plus_N  = extract(elev_rast, cbind(output_grid[,1], output_grid[,2] + dx[2]), method='bilinear')
    elevation_output_grid_minus_N = extract(elev_rast, cbind(output_grid[,1], output_grid[,2] - dx[2]), method='bilinear')

    # Distance between 'plus_N' and 'minus_N' locations
    if(isLonLat_elevation_raster){
        print('Assuming elevation raster horizontal units are in degrees')
        d_distance = distHaversine(cbind(output_grid[,1], output_grid[,2] + dx[2]), cbind(output_grid[,1], output_grid[,2] - dx[2]))
    }else{
        print('Assuming elevation raster horizontal units are in meters')
        d_distance = 2*dx[2]
    }
    output_bathy_slope_N = (elevation_output_grid_plus_N - elevation_output_grid_minus_N)*elevation_raster_scale/d_distance
    rm(elevation_output_grid_plus_N, elevation_output_grid_minus_N, d_distance); gc()

    #
    # East component of slope by central differences
    #
    elevation_output_grid_plus_E = extract(elev_rast, cbind(output_grid[,1]+dx[1], output_grid[,2]), method='bilinear')
    elevation_output_grid_minus_E = extract(elev_rast, cbind(output_grid[,1]-dx[1], output_grid[,2]), method='bilinear')
    if(isLonLat_elevation_raster){
        d_distance = distHaversine(cbind(output_grid[,1]+dx[1], output_grid[,2]), cbind(output_grid[,1]-dx[1], output_grid[,2]))
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

elev_slope = get_elevation_and_gradients_at_output_grid(output_grid, elevation_raster, elevation_raster_scale, isLonLat_elevation_raster); gc()
# Save the elevation and its gradients as a tif for QC
for(nm in names(elev_slope)){
    writeRaster(elev_slope[[nm]], file=paste0('elevation_', nm, '.tif'), options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}

# Loop over all elementary sources and compute the Okada unit source
for(i in 1:length(lon)){

    # Make sure centroid depth is valid
    if(centroid_depth[i] - width[i]/2*sin(dip[i]*deg2rad) < 0){
        stop('Error: Rupture protruding from the earth')
    }

    # Convert output grid to a 'locally cartesian' coordinate
    # system with the rupture lon,lat at the origin. We must
    # do this because Okada's solution is in cartesian coordinates
    output_grid_cartesian = spherical_to_cartesian2d_coordinates(
        coords_lonlat=output_grid,
        origin_lonlat=c(lon[i],lat[i]))

    # Apply okada to the 'cartesian-converted' coordinates
    okada_def = okada_tsunami(
        elon=dlon[i], # Earthquake centroid cartesian x coordinate
        elat=dlat[i], # Earthquake centroid cartesian y coordinate
        edep=centroid_depth[i], # Centroid depth in km!
        strk = strike[i], # Strike
        dip = dip[i], # Dip
        lnth = length[i], # Rupture length (km)
        wdt = width[i], # Rupture width (km)
        disl1 = 1.0 * cos(rake[i]*deg2rad), # Along-strike slip component, unit-slip
        disl2 = 1.0 * sin(rake[i]*deg2rad), # Along-thrust slip component, unit-slip
        rlon = output_grid_cartesian[,1], # Output x coordinate
        rlat = output_grid_cartesian[,2] # Output y coordinate
        )

    #
    # Save all unit source components individually
    # 
    for(comp in c('zdsp', 'ndsp', 'edsp')){
        # Convert the vertical displacement to a raster
        okada_rast = rasterFromXYZ(
            xyz=cbind(output_grid[,1], output_grid[,2], okada_def[[comp]]),
            crs=CRS(output_proj4_string))

        # Save as tif to the current directory
        writeRaster(okada_rast, 
            file=paste0(rupture_basename, '_', rupture_name[i], '_unit_source_', comp, '.tif'),
            format='GTiff',
            options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
    }

    # Compute a unit source that gives the vertical uplift from BOTH zdsp
    # and the horizontal components over a slope
    raster_file_zdsp = paste0(rupture_basename, '_', rupture_name[i], '_unit_source_zdsp.tif')
    raster_file_edsp = paste0(rupture_basename, '_', rupture_name[i], '_unit_source_edsp.tif')
    raster_file_ndsp = paste0(rupture_basename, '_', rupture_name[i], '_unit_source_ndsp.tif')
    unit_source_including_horiz = raster(raster_file_zdsp) +
        # Minus sign because we use gradients of elevation (negative below sea
        # level) rather than depth (as in Tanioka and Satake)
        (-1.0)*(raster(raster_file_edsp)*elev_slope$slope_E + raster(raster_file_ndsp)*elev_slope$slope_N)

    writeRaster(unit_source_including_horiz, 
        file=paste0(rupture_basename, '_', rupture_name[i], '_unit_source_verticalplushorizontal.tif'),
        format='GTiff',
        options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
}


#
# Make a raster for the whole sum
#
all_raster_files_zdsp = paste0(rupture_basename, '_', rupture_name, '_unit_source_zdsp.tif')
all_raster_files_vplush = paste0(rupture_basename, '_', rupture_name, '_unit_source_verticalplushorizontal.tif')
okada_rast_sum = raster(all_raster_files_zdsp[1])*0
okada_rast_sum_zdsp_only = okada_rast_sum * 0
for(i in 1:length(all_raster_files_zdsp)){
    # Z component only
    okada_rast_sum_zdsp_only = okada_rast_sum_zdsp_only + slip[i] * raster(all_raster_files_zdsp[i])

    # Including effect of horizontal components. Minus sign because we use
    # gradients of elevation (negative below sea level) rather than depth (as
    # in Tanioka and Satake)
    okada_rast_sum = okada_rast_sum + slip[i] * raster(all_raster_files_vplush[i])
}

# Save outputs
writeRaster(okada_rast_sum_zdsp_only,
    file=paste0(rupture_basename, '_SUM.tif'), 
    format='GTiff',
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)
writeRaster(okada_rast_sum,
    file=paste0(rupture_basename, '_SUM_including_horizontal_components.tif'), 
    format='GTiff',
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)


