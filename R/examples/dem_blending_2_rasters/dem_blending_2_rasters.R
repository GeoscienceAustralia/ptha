#
# Experimental DEM merging.
#
# We have 2 DEMs -- a spatially extensive, lower-quality one, and another that
# is higher-quality but less spatially extensive. We'd like to make a new DEM
# (say covering the high-quality DEM, plus a more extensive region around it)
# which transitions smoothly between the high and low quality DEMs.
#
# The idea in this script is to interpolate the low-res-DEM error along contour lines. 
#
# From experience with a few examples, I suspect the topography should often be
# corrected along contour lines (e.g. because the shelf occurs too soon, or too
# late, in the low-res-DEM). It seems reasonable to make a DEM merging approach
# along these lines.
#
# Note -- this represents an "ad-hoc guess" as to how the elevation data should
# be corrected. There is no guarantee it will produce accurate results -- we need
# measurements to do that. Case-by-case modification may be required (or just
# use a different approach). 
# 

library(rptha)
rasterOptions(maxmemory = 1e+09)

#
# INPUTS
#

# The "good quality" bathymetry -- this is the standard 10m Victorian coastal bathymetry, reprojected to WGS84
high_res_file = 'Victorian-coast_Bathy_10m_EPSG4326.tif'
# The "low quality but spatially extensive" bathymetry. This is GA250m, converted to Geotif.
low_res_file = 'GA250m/ER_Mapper_ers/ausbath_09_v4_ex_ex.ers.tif'
# Extent of the smooth output raster
smooth_DEM_extent = extent(140.4, 150.4, -39.75, -37.11)
# Filename for smooth output raster
output_tif_file = 'Victoria_smooth_between_GA250_and_Vic10m.tif'

# Define contour_levels.
# We want more contours in shallow water, coarsening further out.
# Here is an ad-hoc approach to get more contours in shallow water.
cl = unique(ceiling(cumsum(2:360)/3) + 0.01) # intermediate variable only
cl = cl[cl < 10000]
# The GA250 data is integer, so no point making contours with sub-integer
# spacing (actually it can lead to artefacts)
contour_levels = c(# Bathymetry
                   -rev(cl),  
                   # Coastline
                   0, 
                   # Topography
                   cl
                   )

# Spacing of points on contours, which are then used for interpolation of the correction
contour_point_spacing = 100 # m

# Define the distance (m) after which we transition linearly to zero correction.
# Our contours may cover long distances in which there is a gap in the high-quality-DEM.
# Over short distances we are happy to interpolate the correction from the other contour points,
# but over long-enough distances the correction should smoothly transition to zero.
# This defines that distance (m).
# For example, if it is 50m, then any point within 50m of a non-gap will use regular interpolation,
# and from 50m-100m we will transition linearly to zero correction (and > 100m will not be corrected).
transition_error_to_zero_distance = 10000

# Part of the code runs in parallel -- use a Fork cluster with this many cores.
MC_CORES = 1 # 12

#
# END INPUTS
#

r_high = raster(high_res_file)
r_GA250 = raster(low_res_file)
r_low = crop(r_GA250, smooth_DEM_extent)

#
# Resample the low-res raster at high resolution
# We use GDAL for speed, so have to mess around with files
#
temp_low_res_raster_file = paste0(tempfile(), '.tif')
writeRaster(r_low, temp_low_res_raster_file, options=c('COMPRESS=DEFLATE'))

r_high_res = res(r_high)
r_high_ext = extent(r_high)
new_temp_high_res_raster_file = paste0(tempfile(), '.tif')
gdal_resample_command = paste0('gdalwarp ', 
    ' -te ', r_high_ext@xmin, ' ', r_high_ext@ymin, ' ', r_high_ext@xmax, ' ', r_high_ext@ymax, 
    ' -tr ', r_high_res[1], ' ', r_high_res[2], 
    ' -r bilinear ', 
    temp_low_res_raster_file, ' ', new_temp_high_res_raster_file)
system(gdal_resample_command)
r_lowhigh = raster(new_temp_high_res_raster_file)
# "Error raster"
r_diff = r_high - r_lowhigh # FIXME: This is slow due to use of "raster" methods
# Cleanup
rm(r_lowhigh)
unlink(new_temp_high_res_raster_file)
unlink(temp_low_res_raster_file)

options("max.contour.segments" = 1e+08)
r_low_contours = rasterToContour(r_low, maxpixels=Inf, levels = contour_levels)
# Points have prescribed spacing in km, and include @data which gives the index
# of the corresponding line in r_low_contours
r_low_contours_p = approxSpatialLines(r_low_contours, spacing=contour_point_spacing/1000, longlat=TRUE, 
    distinguish_disjoint_line_segments=TRUE)

writeOGR(r_low_contours, dsn='tmp_contours', layer='tmp_contours', 
         driver='ESRI Shapefile', overwrite=TRUE)

# Get the error and the original elevation at the points
contour_err = extract(r_diff, r_low_contours_p)
contour_orig = extract(r_low, r_low_contours_p)

# This will be used to interpolate along the contours
fill_contour_NAs <-function(x, y, z, zero_distance=transition_error_to_zero_distance){

    k = which(!is.na(z))
    # Get distance along line, s
    point_mat = cbind(x, y)
    n = dim(point_mat)[1]
    s = distHaversine(point_mat[-1,], point_mat[-n,])
    s = c(0, cumsum(s))

    # Interpolate the correction
    elev_approx = approx(s[k], z[k], xout=s, rule=2)$y

    # If the distance of a point from the s[k], z[k] is too great, then
    # we'd prefer to transition the interpolation to zero
    
    if(length(k) != length(x)){
        # For each point that is not on a high-res-DEM cell, get the distance along the contour
        # to a high-res-DEM-cell point
        not_k = which(is.na(z)) 
        nearest_s_lower = approx(s[k], s[k], method='constant', xout=s[not_k], f=0, rule=2)$y
        nearest_s_upper = approx(s[k], s[k], method='constant', xout=s[not_k], f=1, rule=2)$y
        not_k_distance = pmin(abs(s[not_k] - nearest_s_lower), abs(s[not_k] - nearest_s_upper))
        # Make a correction. The following number = 
        #     0 for not_k_distance < transition_error_to_zero_distance
        #     varies linearly from 0 to 1 with an additional transition_error_to_zero_distance
        not_k_correction = pmin(1, pmax(not_k_distance - zero_distance, 0)/zero_distance)
        # Update the correction, so we smoothly transition to zero far from the high-res DEM
        elev_approx[not_k] = elev_approx[not_k] * (1 - not_k_correction)
    }

    return(elev_approx)
}

# Do the interpolation along each connected contour
contour_err_filled = contour_err # Hold interpolated values
r_low_contours_p_coords = coordinates(r_low_contours_p) # Useful
unique_IDS = unique(r_low_contours_p@data[,1]) # Loop over each of these

# Convenient to run the contour interpolation in parallel
parallel_fun<-function(i){
    # Find points on this contour
    k = which(r_low_contours_p@data[,1] == unique_IDS[i])
    out = contour_err_filled[k]
    # Do the interpolation, so long as there is enough points, and some NA values.
    if(sum(!is.na(contour_err[k])) < 2){
        # Need at least 2 non-NA points to do the interpolation
        out = 0
    }else if(sum(is.na(contour_err[k])) == 0){
        # next
    }else{
        out = fill_contour_NAs(r_low_contours_p_coords[k,1], 
            r_low_contours_p_coords[k,2], contour_err[k])
    }
    # Return values and indices
    return(cbind(k, out))
}

# Setup cluster
library(parallel)
cl = makeForkCluster(nnodes=MC_CORES)
# Load libraries on cluster
clusterEvalQ(cl=cl, {
    library(rptha)
    rasterOptions(maxmemory = 1e+09)
            })

# Export data to cluster
clusterExport(cl, 
              varlist=c('unique_IDS', 'r_low_contours_p', 'contour_err_filled', 
                        'contour_err', 'fill_contour_NAs', 'transition_error_to_zero_distance'))
# Run in parallel
# Beware fairly heavy memory
contour_filled_inds_values = parLapply(cl=cl, 1:length(unique_IDS), parallel_fun)
stopCluster(cl)

# Pack into contour_err_filled
for(i in 1:length(contour_filled_inds_values)){
    k = contour_filled_inds_values[[i]][,1]
    contour_err_filled[k] = contour_filled_inds_values[[i]][,2]
}

# Make a raster with the dimensions of r_low, with cells = approximation error
# Use delaunay triangulation to do the interpolation
r_low_xyz = rasterToPoints(r_low)

# To do delaunay triangulation, the coordinates should be close to cartesian.
# Over small areas of the earth, we can do this by centering the longitudes and
# then adjusting by cos(median_latitude/180*pi)
median_longitude = median(r_low_xyz[,1])
median_latitude = median(r_low_xyz[,2])
long_stretch = cos(median_latitude/180*pi)
r_low_contours_p_coords_stretched = r_low_contours_p_coords[,1:2]
r_low_contours_p_coords_stretched[,1] = (r_low_contours_p_coords_stretched[,1] - median_longitude)*long_stretch
r_low_xyz_coords_stretched = r_low_xyz[,1:2]
r_low_xyz_coords_stretched[,1] = (r_low_xyz_coords_stretched[,1] - median_longitude) * long_stretch

# Finally do the delauny triangulation, using the stretched coordinates
approx_err = triangular_interpolation(xy=r_low_contours_p_coords_stretched[,1:2], vals=contour_err_filled, 
    newPts=r_low_xyz_coords_stretched[,1:2], useNearestNeighbour=FALSE)

# Apply the approximation error
new_xyz = r_low_xyz
new_xyz[,3] = new_xyz[,3] + approx_err

# Finally, where we have high-res data, we should use it
on_high_res = extract(r_high, new_xyz[,1:2])
k = which(!is.na(on_high_res))
new_xyz[k,3] = on_high_res[k]

# Write to file
new_rast_low = rasterFromXYZ(new_xyz, res=res(r_low), crs=CRS(proj4string(r_low)))
writeRaster(new_rast_low, file=output_tif_file, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)


