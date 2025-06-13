#' Test the tidal adjustment tif at specified points
library(terra)


raster <- rast("tpxo9_adjusted_for_gauges.tif")
test_points <- read.table("x_y_z.csv", skip=1, header=TRUE, sep=",")

# extract the raster values at the (points$x, points$y) locations
locations <- data.frame(x=test_points$x, y=test_points$y)
print(locations)
test_points$z_tif <- extract(raster, locations, ID=FALSE)

# compare the extracted values with the z values
test_points$diff <- -test_points$z - test_points$z_tif$z

# assert that all differences are below a tolerance
tolerance <- 0.001  # m
stopifnot(all(abs(test_points$diff) < tolerance))
print("All differences are below tolerance. Maximum difference:")
print(max(abs(test_points$diff)))