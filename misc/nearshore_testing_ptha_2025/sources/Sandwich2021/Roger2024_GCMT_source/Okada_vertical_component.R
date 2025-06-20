#
# Code to compute rupture from multiple earthquake sources, 
# using the source model from
# Roger, J. H. M.; Jamelot, A.; Hébert, H.; Power, W.; Gusman, A. & Thomas, B.
# E. O. The South Sandwich Tsunami of 12 August 2021: An Underestimated
# Widespread Tsunami Hazard Around the World Journal of Geophysical Research:
# Oceans, American Geophysical Union (AGU), 2024, 129

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

#
# Input parameters -- rupture geometry
#
lat = c(-58.03, -59.52, -60.42)
lon = c(-23.81, -23.79, -25.08) # rupture location in DEGREES E, WGS84, centroid
width = c(57.2,   58.1,   68.1) # rupture width in km
length=c(114.4,  116.2,  136.2) # rupture length in km
strike=c(168.0,  177.0,  216.0) # Strike in degrees clockwise from N
dip = c(20, 24, 17)
depth = c(20, 20, 20) # Depth at centroid
rake  = c(93, 78, 104) # Degrees -- rupture rake
slip = c(2.49, 2.53, 2.96)
rupture_name   = as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'RogerGCMT_2024_Sandwich2021' # Name at the start of output files

#
# Flags to denote meaning of lon, lat and depth
#
# Lon/Lat location: 
# Use -1 for 'middle of bottom edge of rupture', 0 for 'centroid', and 1 for
# 'middle of top edge of rupture'
lon_lat_location = 0 # This is a guess, not clear from the paper
# Depth location:
# Use -1 for 'depth at bottom edge of rupture', 0 for 'depth at centroid', and
# 1 for 'depth at top edge of rupture'.
depth_location = 0 # A guess, not clear from paper

#
# Input parameters -- output grid
#
output_res = 2/60 # 
output_lons = seq(-30, -20, by=output_res)
output_lats = seq(-62, -54, by=output_res)
output_proj4_string = "+init=epsg:4326" # code for WGS84 

#
# END INPUT PARAMETERS
#

# Make a 'grid' of points for output
output_grid = expand.grid(output_lons, output_lats)

deg2rad = pi/180

# Offsets to make get lon/lat/depth at centroid
dlon = width/2 * cos(dip*deg2rad) * cos((strike)*deg2rad) * 1000 * lon_lat_location
dlat = -width/2 * cos(dip*deg2rad) * sin((strike)*deg2rad) * 1000 * lon_lat_location
# Define centroid depths
centroid_depth = depth + width/2 * sin(dip*deg2rad) * depth_location

# Output data in a convenient way
rup2csv<-function(){

    # Estimate the centroid longitude/latitude
    newlon = lon + destPoint(cbind(0, lat), b=90, d=dlon, f=0)[,1]
    newlat = destPoint(cbind(lon, lat), b=0,  d=dlat, f=0)[,2]

    output = data.frame(lonc=newlon, latc=newlat, depthc=centroid_depth,
        length=length, width=width, strike=strike, rake=rake, slip=slip)

    write.csv(output, 'rupture_summary.csv', row.names=FALSE)
}
rup2csv()


# Name for each unit-source
unit_source_raster_filenames = paste0(rupture_basename, '_', 
    rupture_name, '_unit_source.tif')

# Loop over all ruptures
for(i in 1:length(lon)){


    # Make sure centroid depth is valid
    if(centroid_depth[i] - width[i]/2*sin(dip[i]*deg2rad) < 0){
        stop('Error: Rupture protruding from the earth')
    }

    # Convert output grid to a 'locally cartesian' coordinate
    # system with the rupture centroid at the origin. We must
    # do this because Okada's solution is in cartesian coordinates
    output_grid_cartesian = spherical_to_cartesian2d_coordinates(
        coords_lonlat=output_grid,
        origin_lonlat=c(lon[i],lat[i]))

    # Apply okada to the 'cartesian-converted' coordinates. Assume
    # unit slip, so we compute unit-sources. Scale the results later. 
    unit_slip = 1.0
    okada_def = okada_tsunami(
        elon=dlon[i], # Earthquake centroid cartesian x coordinate
        elat=dlat[i], # Earthquake centroid cartesian y coordinate
        edep=centroid_depth[i], # Centroid depth in km!
        strk = strike[i], # Strike
        dip = dip[i], # Dip
        lnth = length[i], # Rupture length (km)
        wdt = width[i], # Rupture width (km)
        disl1 = unit_slip * cos(rake[i]/180 * pi), # Along-strike slip component, unit-slip
        disl2 = unit_slip * sin(rake[i]/180 * pi), # Along-thrust slip component, unit-slip
        rlon = output_grid_cartesian[,1], # Output x coordinate
        rlat = output_grid_cartesian[,2] # Output y coordinate
        )

    #
    # Save all rasters individually
    # 

    # Convert the vertical displacement to a raster
    okada_rast = rasterFromXYZ(
        xyz=cbind(output_grid[,1], output_grid[,2], okada_def$zdsp),
        crs=CRS(output_proj4_string))

    # Save as tif to the current directory
    writeRaster(okada_rast, 
        file=unit_source_raster_filenames[i],
        format='GTiff',
        options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

}

# Make a raster for the whole source

okada_rast_sum = raster(unit_source_raster_filenames[1])*0
for(i in 1:length(unit_source_raster_filenames)){
    okada_rast_sum = okada_rast_sum + slip[i] * raster(unit_source_raster_filenames[i])
}

writeRaster(okada_rast_sum,
    file=paste0(rupture_basename, '_SUM.tif'), 
    format='GTiff',
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

