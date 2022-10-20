#
# For this problem the bathymetry needs to be constructed with a Matlab script.
#
# I don't have Matlab. But the problem is straightforward enough.
#     0. Do some datum adjustments of xyz data
#     1. Grid the xyz data
#     2. Add in the shoal
#
library(rptha)

# Read xyz file provided with the benchmark description
xyz = read.csv('TWB 071708_5CM_XYZ_METERS.txt', header=FALSE)

# There is an offset applied, based on the mean of points where x < 10.2
vertical_offset = mean(xyz[xyz[,1] < 10.2, 3])

zero_survey_5cm = xyz
zero_survey_5cm[,3] = zero_survey_5cm[,3] - vertical_offset

# They grid at a given number of points. 
# npts = 1000
#xseq = seq(min(zero_survey_5cm[,1]), max(zero_survey_5cm[,1]), length=npts+1)    
#yseq = seq(min(zero_survey_5cm[,2]), max(zero_survey_5cm[,2]), length=npts+1)
# To save filespace I will do 5cm gridding
#xseq = seq(min(zero_survey_5cm[,1]), max(zero_survey_5cm[,1]), by=0.05)    
#yseq = seq(min(zero_survey_5cm[,2]), max(zero_survey_5cm[,2]), by=0.05)    

# Useful to extend the computational domain so that boundary conditions can be imposed well
# See Macais et al. 2020
xseq = seq(-9, 44.6, by=0.05)    
yseq = seq(-13, 13, by=0.05)    

target_xy = expand.grid(xseq, yseq)

# Grid the xyz data -- missing the cone
zvals = triangular_interpolation(zero_survey_5cm[,1:2], zero_survey_5cm[,3], target_xy, useNearestNeighbour=TRUE)
bg_rast = rasterFromXYZ(cbind(target_xy, zvals))

# Add in the cone
xc = 17
yc = 0
hs = 0.45
rs = 3.0
xy_rad = sqrt((target_xy[,1] - xc)**2 + (target_xy[,2]-yc)**2)
cone_vals = hs * pmax(0, 1 - xy_rad/rs)
cone_rast = rasterFromXYZ(cbind(target_xy, cone_vals))

# Combine the two datasets
combined_rast = bg_rast + cone_rast

writeRaster(combined_rast, 'bathy_with_cone.tif', options=c('COMPRESS=ZSTD', 'PREDICTOR=3'), overwrite=TRUE)
#writeRaster(combined_rast, 'bathy_with_cone.tif', options=c('COMPRESS=DEFLATE', 'PREDICTOR=3'), overwrite=TRUE)
