# The file 'merged_gebco_GA250_dem.tif' uses GA250 wherever it exists, and GEBCO2014 elsewhere.
# However, there are a few locations where GA250 is clearly wrong, and GEBCO2014 is better.
# Thus we 'patch' the merged DEM in these locations with GEBCO2014

input_args = c(
    'merged_gebco_ga250_dem.tif', 
    '../GEBCO_2014_1m/GEBCO_2014_1minx1min_W-39.9958333-E320.0041667.tif', 
    'patch_polygon')

library(rgdal)
library(raster)

# Tricks to avoid blowing out the JOBFS space on NCI
dir.create('tmp', showWarnings=FALSE)
rasterOptions(tmpdir='tmp')
Sys.setenv(TMPDIR = 'tmp')

r1 = raster(input_args[1])
r2 = raster(input_args[2])
p1 = readOGR(dsn=input_args[3], layer=input_args[3])

# Find indices of cells in r1 which should be updated
cells_to_patch = extract(r1, p1, cellnumbers=TRUE)

# Loop over all patch polygons, and apply the patch
for(i in 1:length(p1)){
    cells_to_patch_xy = xyFromCell(r1, cells_to_patch[[i]][,1])
    new_cell_values = extract(r2, cells_to_patch_xy)
    r1[cells_to_patch[[i]][,1]] = new_cell_values
}

writeRaster(r1, 'merged_gebco_ga250_dem_patched.tif', format='GTiff', options=c('COMPRESS=DEFLATE'))

