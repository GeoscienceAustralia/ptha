#
# Code to compute rupture from multiple earthquake sources, 
# using the S2 source model from
# Wei, Y.; Fritz, H. M.; Titov, V. V.; Uslu, B.; Chamberlin, C. & Kalligeris, N. Source Models and Near-Field Impact of the 1 April 2007 Solomon Islands Tsunami Pure and Applied Geophysics, Springer Science and Business Media LLC, 2015, 172, 657–682

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

mydat = read.table('Wei_S2_data.txt', header=TRUE)
#
# Input parameters -- rupture geometry
#
lat = mydat$lat # rupture location in DEGRESS N, WGS84, bottom-edge of fault
lon = mydat$lon # rupture location in DEGREES E, WGS84, bottom-edge of fault
width = mydat$width # rupture width in km
length = mydat$length # rupture length in km
strike = mydat$strike  # Strike in degrees clockwise from N
dip = mydat$dip # Dip
depth = mydat$depth # Depth at top-edge of fault
rake  = mydat$rake # Degrees -- rupture rake
slip = mydat$slip
rupture_name   = as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'Wei_S2_Solomon2007_source' # Name at the start of output files

#
# Flags to denote meaning of lon, lat and depth
#
# Lon/Lat location: 
# Use -1 for 'middle of bottom edge of rupture', 0 for 'centroid', and 1 for
# 'middle of top edge of rupture'
lon_lat_location = -1
# Depth location:
# Use -1 for 'depth at bottom edge of rupture', 0 for 'depth at centroid', and
# 1 for 'depth at top edge of rupture'.
depth_location = 1

#
# Input parameters -- output grid
#
output_res = 2/60 # 
output_lons = seq(153, 160, by=output_res)
output_lats = seq(-12, -4, by=output_res)
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

