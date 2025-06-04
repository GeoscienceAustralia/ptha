#' Read the tpxo9 raster, project, lower and negate the elevation, fill empty cells westwards, set easterly positive cells to NA and smoothly deform the raster to fit the x_y_z data.
library(terra)


get_background_sea_level <- function(raster){
    # Get mean cell value on the eastern boundary (last column)
    eastern_mean <- mean(values(raster, col=ncol(raster), ncol=1), na.rm=TRUE)
    return(eastern_mean)
}

# Create a Gaussian to superimpose on the tpxo9 raster to adjust one gauge.
get_gaussian_bump <- function(x, y, dz, raster, sigma_km=20) {
    dist <- rast(raster)
    point <- vect(cbind(x, y), crs=crs(raster))
    dist <- distance(dist, point, unit="km")
    gaussian <- dz * exp(-dist^2 / (2 * sigma_km^2))
    return(gaussian)
}

adjust_to_gauges <- function(raster, xyz_file){
    xyz <- read.table(xyz_file, skip=4, header=TRUE, sep=",")

    xyz$offset <- as.numeric(xyz$offset)
    xyz$offset[is.na(xyz$offset)] <- 0

    mse <- Inf
    a_tol <- 0.0000001
    while (mse > a_tol) {
        for(i in 1:nrow(xyz)){
            # Note raster is negated (-- = +)
            dz <- xyz[i, "z"] + xyz[i, "offset"] + unlist(extract(raster, cbind(xyz[i, "x"], xyz[i, "y"])))
            gaussian <- get_gaussian_bump(xyz[i, "x"], xyz[i, "y"], dz, raster, sigma_km=xyz[i, "sigma_km"])
            raster <- raster - gaussian
        }

        mse <- mean((xyz$z + xyz$offset + unlist(extract(raster, cbind(xyz$x, xyz$y))))^2)
        print(paste0("Mean Squared Error to Gauges: ", mse))
    }

    return(raster)
}

fill_west <- function(raster){
    for(j in 1:nrow(raster)){
        row <- values(raster, row=j, nrows=1)
        # first !is.na value from the west
        first_numeric <- which(!is.na(row))[1]
        is_some_leading_na <- (!is.na(first_numeric) & first_numeric > 1)
        if (is_some_leading_na){
            first_val <- row[first_numeric]
            row[is.na(row)] <- first_val
            raster[j,] <- row
        }
    }
    return(raster)
}

lin_interpolate <- function(raster, inner_ext, background_sea_level){
    # create a vect from ext
    buffer_vec <- vect(inner_ext, crs=crs(raster))
    inner_raster <- mask(raster, buffer_vec)
    dist_to_inner <- distance(inner_raster)
    
    # max distance is in corner cases. Use a sqrt(2) factor to scale the distance to the background sea level.
    max_dist <- max(as.matrix(dist_to_inner), na.rm=TRUE)
    ratio_of_background <- min((dist_to_inner / max_dist * sqrt(2)), 1.0)

    outer_raster <- raster + (background_sea_level - raster) * ratio_of_background
    outer_raster <- mask(outer_raster, buffer_vec, inverse=TRUE)
    inner_raster <- subst(inner_raster, NA, 0)
    outer_raster <- subst(outer_raster, NA, 0)
    raster <- inner_raster + outer_raster
    return(raster)
}

# Read the tpxo9 raster and negate
filename <- "/g/data/w85/tsunami/MODELS/tidal_prediction/TPXO9/get_tidal_extrema_on_grid/outputs/max_tidal_stage.tif"
tpxo9 <- -rast(filename)

print("Fill westward")
tpxo9 <- fill_west(tpxo9)

print("Get the background sea level")
background_sea_level <- get_background_sea_level(tpxo9)

# print("Adjust the raster to fit the x_y_z data at the gauges and AHD at points.")
# xyz_file <- "x_y_z.csv"
# tpxo9 <- adjust_to_gauges(tpxo9, xyz_file)

print("Linearly decay the raster in an internal buffer zone of 1 degree to the background sea level.")
extent <- ext(tpxo9)
buffer <- 1
buffer_ext <- ext(c(extent$xmin+buffer, extent$xmax-buffer, extent$ymin+buffer, extent$ymax-buffer))
tpxo9 <- lin_interpolate(tpxo9, buffer_ext, background_sea_level)

print("Extend raster to match model extent with a 1 pixel buffer. See ../swals/model_multidomain_design_mod.f90")
x_buffer <- xres(tpxo9)
y_buffer <- yres(tpxo9)
full_extent <- ext(100-x_buffer, 320+x_buffer, -79-y_buffer, 68+y_buffer)
tpxo9 <- extend(tpxo9, full_extent)

print("replace all NA with background sea level")
tpxo9[is.na(tpxo9)] <- background_sea_level

print(tpxo9)
writeRaster(tpxo9, "tpxo9_adjusted_for_gauges.tif", overwrite=TRUE)
