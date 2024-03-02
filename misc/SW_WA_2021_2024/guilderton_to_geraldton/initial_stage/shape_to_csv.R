library(terra)

# get all shapes
shape_files <- list.files(path = "shapes", pattern = "\\.shp$", full.names = TRUE)

# create a csv of x, y coords from a given shapefile
create_csv <- function(shape_file) {
    spat_vec <- vect(shape_file)
    coords <- crds(spat_vec)
 
    basename <- tools::file_path_sans_ext(shape_file)
    file_write <- paste0(basename, ".csv")
    write.csv(coords, file = file_write, row.names = FALSE)
}

lapply(shape_files, create_csv)
