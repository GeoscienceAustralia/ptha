#' Check that the input coordinates and log them
#' 
library(terra)


input_shapes = c(Sys.glob("../../breakwalls/*/*/shp"), Sys.glob("../../inverts/*/*/shp"))
for (i in seq_along(input_shapes)) {
    x <- vect(input_shapes[i])
    writeLines(c(x, crs(x)), "input_coords.txt")
}

input_rasters = c(readLines("../../elevation/swals_elevation_rasters.txt"), "../../friction/*.tif")
for (i in seq_along(input_rasters)) {
    x <- rast(input_rasters[i])
    writeLines(c(x, crs(x)), "input_coords.txt")
}
