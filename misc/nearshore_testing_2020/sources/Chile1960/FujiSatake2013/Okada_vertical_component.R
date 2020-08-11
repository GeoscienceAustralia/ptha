#
# Code to compute rupture from multiple earthquake sources, using
# example from 
# Fujii, Y. & Satake, K. Slip Distribution and Seismic Moment of the 2010 and 1960 Chilean Earthquakes Inferred from Tsunami Waveforms and Coastal Geodetic Data Pure and Applied Geophysics, 2013, 170, 1493-1509 

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

#
# Input parameters -- rupture geometry
#
lat = -c(45.2, 44.30611, 43.41222, 42.51832, 41.62443, 40.73054, 39.83665, 38.94276, 38.04887, 45.25488, 44.36099, 43.46709, 42.5732, 41.67931, 40.78542, 39.89153, 38.99764, 38.10374, 45.30976, 44.41586, 43.52197, 42.62808, 41.73419, 40.8403, 39.94641, 39.05251, 38.15862) # rupture location in DEGRESS N, WGS84
lon = 360 - c(75.9, 75.75513, 75.61026, 75.46538, 75.32051, 75.17564, 75.03077, 74.8859, 74.74102, 75.31006, 75.16518, 75.02031, 74.87544, 74.73057, 74.58569, 74.44082, 74.29595, 74.15108, 74.72011, 74.57524, 74.43037, 74.28549, 74.14062, 73.99575, 73.85088, 73.70601, 73.56113) # rupture location in DEGREES E, WGS84 
width          =   rep(50, length(lon)) # km -- rupture width
length         =   rep(100, length(lon)) # km -- rupture width
dip            =   rep(20, length(lon)) # Degrees -- rupture dip
rake           =  rep(105, length(lon))  # Degrees -- rupture rake
depth = c( rep(0, 9), rep(17.1, 9), rep(34.2, 9) )# km -- depth of rupture (top)
strike         =    rep(7, length(lon)) # Degrees 
## Slip from PURE TSUNAMI inversion
## slip = c(1.71, 0, 1.44, 1.68, 4.62, 0, 9.25, 3.88, 0.3, 0, 0, 5.75, 5.46, 0, 1.98, 10.21, 12.1, 0.43, 90.01, 17.36, 17.48, 8.87, 40.08, 0, 53.22, 4.14, 0)
## Slip from TSUNAMI AND GEODETIC inversion
slip = c(3.05, 12.76, 21.42, 6.48, 4.11, 2.34, 5.54, 4.16, 0.76, 17.09, 18.12, 14.81, 27.38, 24.69, 25.93, 30.07, 10.76, 1.14, 2.26, 4.55, 5.14, 5.51, 9.59, 9.41, 17.53, 1, 1.83)
rupture_name   =    as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'Fuji_chile1960_sources' # Name at the start of output files
make_plot = FALSE

# Adjust lon/lat to be on the central top edge
new_lonlat = destPoint(cbind(lon, lat), b=strike, d=length/2*1000)
lon = new_lonlat[,1] + 360
lat = new_lonlat[,2]

if(make_plot){
    # Quick check that the slip matches the picture in the paper
    plot(lon, lat, col=colorRampPalette(c('white', 'yellow', 'orange', 'red', 'black'))(300)[floor(slip*10)], asp=1, pch=19, cex=2)
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
output_lons = seq(282, 290, by=output_res)
output_lats = seq(-47,  -32, by=output_res)
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


#
# Plot all individual raster files
#
if(make_plot){
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
