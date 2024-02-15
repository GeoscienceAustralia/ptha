#
# Code to compute rupture from multiple earthquake sources, using
# example from 
# Lorito, S., F. Romano, S. Atzori, X. Tong, A. Avallone, J. McCloskey, M. Cocco, E. Boschi, and A. Piatanesi (2011). Limited overlap between the seismic gap and coseismic slip of the great 2010 Chile earthquake, Nature Geoscience, 4, 173–177, doi:10.1038/ngeo1073 

#
# Approach: We compute the rupture for each source separately.
# We then sum the deformation to get the final result
# 

library(rptha)

chile2010 = read.csv('LoritoEtAl2011_source_model.csv')
#
# Input parameters -- rupture geometry
#
lat = chile2010$lat # rupture location in DEGRESS N, WGS84
lon = chile2010$lon # rupture location in DEGREES E, WGS84 
depth = chile2010$depth + 0.0
width = chile2010$width + 0.0
length = chile2010$length + 0.0
slip = chile2010$slip + 0.0
strike = chile2010$strike + 0.0
dip = chile2010$dip + 0.0
rake = chile2010$rake + 0.0
rupture_name   =    as.character(1:length(lon)) # A name, which will appear in the filename
rupture_basename = 'Lorito_chile2010_sources' # Name at the start of output files
make_plot = FALSE

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
lon_lat_location = 0
# Depth location:
# Use -1 for 'depth at bottom edge of rupture', 0 for 'depth at centroid', and
# 1 for 'depth at top edge of rupture'.
depth_location = 0

#
# Input parameters -- output grid
#
output_res = 2/60 # 5 arc-minute resolution, consistent with the study
output_lons = seq(282, 290, by=output_res)
output_lats = seq(-44,  -28, by=output_res)
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
