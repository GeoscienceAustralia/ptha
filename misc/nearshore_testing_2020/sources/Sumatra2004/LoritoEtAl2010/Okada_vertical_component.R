
# Code to compute rupture from multiple earthquake sources, using
# example from 
# Lorito, S.; Piatanesi, A.; Cannelli, V.; Romano, F. & Melini, D. Kinematics and source zone properties of the 2004~Sumatra-Andaman earthquake and tsunami: Nonlinear joint inversion of tide gauge, satellite altimetry, and GPS data Journal of Geophysical Research, Wiley-Blackwell, 2010, 115

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

#
# Input parameters -- rupture geometry
#
#
lat = c(2.671, 1.915, 3.312, 2.583, 4.206, 3.667, 5.473, 5.133, 7.051, 6.750, 8.416, 8.1667, 10.085, 10.075, 11.637, 11.750, 13.067, 13.250) # rupture location in DEGRESS N, WGS84

lon = c(95.845, 95.210, 94.964, 94.300, 94.183, 93.400, 93.582, 92.850, 93.084, 92.350, 92.524, 91.733, 92.089, 91.312, 92.259, 91.562, 92.602, 91.988) # rupture location in DEGREES E, WGS84 
width =   c(113.728, 110.273, 113.624, 110.137, 111.504, 106.331, 103.985, 90.568, 101.266, 89.367, 101.569, 93.561, 98.863, 88.3, 95.999, 80.698, 92.323, 73.763) # rupture width in km
length = c( 137.154, 137.211, 109.662, 123.126, 125.25, 158.965, 182.4, 183.497, 163.41, 172.928, 176.108, 181.21, 167.488, 204.321, 156.887, 166.563, 167.403, 172.847) # rupture length in km
strike = c( 301.8, 301.8, 309.14, 309.5, 328.16, 329.9, 336.16, 341.05, 342.91, 341.39, 334.56, 334.66, 357.33, 0.77, 10.71, 12.08, 15.94, 17.23) # Strike in degrees clockwise from N
dip = c(13.22, 5.2, 13.23, 5.21, 13.48, 5.4, 14.48, 6.34, 14.88, 6.42, 14.83, 6.14, 15.25, 6.5, 15.71, 7.12, 16.36, 7.79) # Dip
depth = c(10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1, 10.1, 0.1) # Depth at upper edge
# NOTE: here I force rupture to the trench, although it's not what was done in the paper. In practice the finite-depth will
# produce a 'spike' near the trench (albeit that might not be resolved in this grid).
depth[depth<0.2] = 0
rake           =  c( rep(95, 12), rep(135, 6) ) # Degrees -- rupture rake
## Slip from PURE TSUNAMI inversion
## slip = c(1.71, 0, 1.44, 1.68, 4.62, 0, 9.25, 3.88, 0.3, 0, 0, 5.75, 5.46, 0, 1.98, 10.21, 12.1, 0.43, 90.01, 17.36, 17.48, 8.87, 40.08, 0, 53.22, 4.14, 0)
## Slip from TSUNAMI AND GEODETIC inversion
slip = c(1, 21, 34, 0, 22, 33, 8, 16, 3, 0, 12, 18, 3, 1, 8, 0, 6, 0)
rupture_name   =    as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'Lorito_sumatra2004_sources' # Name at the start of output files
make_plot = FALSE

if(make_plot){
    # Quick check that the slip matches the picture in the paper
    plot(lon, lat, 
         col=colorRampPalette(c('white', 'yellow', 'orange', 'red', 'black'))(341)[floor(slip*10)+1], 
         asp=1, pch=19, cex=2)
    points(lon, lat, cex=2)
}

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

#
# Input parameters -- output grid
#
output_res = 2/60 # 5 arc-minute resolution, consistent with the study
output_lons = seq(89, 97, by=output_res) + 2
output_lats = seq(0,  15, by=output_res)
output_proj4_string = "+init=epsg:4326" # code for WGS84 

#
# Input parameters -- multi-source png figure
#
png_width = 12
png_height= 8
png_xlim = range(output_lons)
png_ylim = range(output_lats) 



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
        length=length, width=width, strike=strike, dip=dip, rake=rake, slip=slip)

    write.csv(output, 'rupture_summary.csv', row.names=FALSE)
}
rup2csv()


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

    # Apply okada to the 'cartesian-converted' coordinates
    okada_def = okada_tsunami(
        elon=dlon[i], # Earthquake centroid cartesian x coordinate
        elat=dlat[i], # Earthquake centroid cartesian y coordinate
        edep=centroid_depth[i], # Centroid depth in km!
        strk = strike[i], # Strike
        dip = dip[i], # Dip
        lnth = length[i], # Rupture length (km)
        wdt = width[i], # Rupture width (km)
        disl1 = slip[i] * cos(rake[i]/180 * pi), # Along-strike slip component, unit-slip
        disl2 = slip[i] * sin(rake[i]/180 * pi), # Along-thrust slip component, unit-slip
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
        file=paste0(rupture_basename, '_', rupture_name[i], '.tif'),
        format='GTiff',
        options=c('COMPRESS=DEFLATE'), overwrite=TRUE)

}

# Make a raster for the whole sum
all_raster_files = paste0(rupture_basename, '_', rupture_name, '.tif')

okada_rast_sum = raster(all_raster_files[1])*0
for(i in 1:length(all_raster_files)){
    okada_rast_sum = okada_rast_sum + raster(all_raster_files[i])
}

writeRaster(okada_rast_sum,
    file=paste0(rupture_basename, '_SUM.tif'), 
    format='GTiff',
    options=c('COMPRESS=DEFLATE'), overwrite=TRUE)


if(make_plot){
    #
    # Plot all individual raster files
    #
    rasts = lapply(as.list(all_raster_files), raster) # Read all rasters into 1 list
    names(rasts) = all_raster_files 
    rasts_range = range(unlist(lapply(rasts, f<-function(x) range(as.matrix(x)))))
    n = length(all_raster_files)
    nr = floor(sqrt(n)) # = 2 if n >1, otherwise 1
    nc = ceiling(n/nr) # Number of columns in plot

    png(paste0(rupture_basename, '_quick_plot.png'), 
        width=png_width, height=png_height, 
        units='in', res=300)
    par(mfrow=c(nr, nc)) # Divide into rows, columns
    par(mar=c(3,4,4,5)) # Adjust individual plot margins
    for(i in 1:length(rasts)){
        plot(rasts[[i]], main=all_raster_files[i],
            col=colorRampPalette(
                c('blue', 'green', 'yellow', 'red'))(100),
            zlim=rasts_range, 
            xlim=png_xlim, ylim=png_ylim)
    }
    dev.off()

    png('sum_raster_quick_plot.png', width=8, height=6, units='in', res=200)
    plot(okada_rast_sum, col=colorRampPalette(
        c('blue', 'green', 'yellow', 'red'))(100))
    #text(lon, lat, rupture_name, col='white')
    dev.off()

}
