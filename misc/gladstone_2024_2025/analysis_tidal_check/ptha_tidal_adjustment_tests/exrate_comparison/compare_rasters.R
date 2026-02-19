# Compute the relative and absolute difference between each corresponding raster in two direcories
library(terra)


compare_rasters <- function(dir1, dir2, output_dir) {
    # Get the list of rasters in each directory
    rasters1 <- list.files(dir1, pattern = "\\.tif$", full.names = TRUE)
    rasters2 <- list.files(dir2, pattern = "\\.tif$", full.names = TRUE)

    # Check that the number of rasters in each directory is the same
    if (length(rasters1) != length(rasters2)) {
        stop("The number of rasters in each directory is not the same")
    }

    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # Loop over the rasters and compute the differences
    for (i in seq_along(rasters1)) {
        raster1 <- rast(rasters1[i])
        raster2 <- rast(rasters2[i])

        # Compute the relative difference
        # relative_diff <- abs(raster1 - raster2) / abs(raster1)

        # # Compute the absolute difference
        # absolute_diff <- abs(raster1 - raster2)

        # Compute the difference
        diff <- raster1 - raster2

        # Write the differences to disk
        filename <- basename(rasters1[i])
        # writeRaster(relative_diff, file.path(output_dir, paste0("relative_diff_", filename)))
        # writeRaster(absolute_diff, file.path(output_dir, paste0("absolute_diff_", filename)))
        writeRaster(diff, file.path(output_dir, paste0("diff_", filename)))
    }
}

# Set the directories for comparison
dir1 <- "sea_level_vary/highres_depth_with_variance/sea_level_vary-depth-LogicTreeMean-sum_of_source_zones"
dir2 <- "sea_level_1.459/highres_depth_with_variance/sea_level_1.459-depth-LogicTreeMean-sum_of_source_zones"
output_directory <- "sea_level_1.459/highres_depth_with_variance/sea_level_1.459-depth-LogicTreeMean-sum_of_source_zones-comparison"

# Compare the rasters
compare_rasters(dir1, dir2, output_directory)
