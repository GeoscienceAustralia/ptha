library(OpenStreetMap)
library(raster)

# Get the OSM data, with caching to speed it up
get_osm_backdrop_reproj <- function(lower_left, upper_right, dbox, backdrop, zoom) {
    osm_cachefile <- paste0(
        "osm_tile_cache/",
        upper_right[1], "_",
        upper_right[2], "_",
        lower_left[1], "_",
        lower_left[2], "_",
        "zoom_", zoom,
        ".RDS"
    )
    print(osm_cachefile)
    print(file.exists(osm_cachefile))

    # if the exact cache file exists
    if (file.exists(osm_cachefile)) {
        osm_backdrop_reproj <- readRDS(osm_cachefile)
        return(osm_backdrop_reproj)
    }

    # # first check if we can get it from a larger cached
    # candidate_osm_files <-
    # for (file in candidate_osm_files) {
    #     osm_backdrop_reproj <- readRDS(file)
    #     requested_extent <- ext(lower_left[1], upper_right[1], lower_left[2], upper_right[2])

    # }

    # # Check if cached extent covers the requested extent
    # if (ext(osm_backdrop_reproj) >= requested_extent) {
    #     # Crop to requested extent
    #     cropped_backdrop <- crop(osm_backdrop_reproj, requested_extent)
    #     return(cropped_backdrop)
    # }


    if (!file.exists(osm_cachefile)) {
        # Get a piece that is bigger than we need, to avoid gaps in plots when
        # the plot limits are extended by R (which depends on the figure size
        # too)
        print("Downloading from Open Street Map.")

        # add buffer around. Swapping e.g. upper_right for upperLeft
        if (length(dbox) == 4) {
            # assumes dbox as c(xmin, xmax, ymin, ymax)
            upperLeft <- c(upper_right[2] + dbox[4], lower_left[1] - dbox[1])
            lowerRight <- c(lower_left[2] - dbox[3], upper_right[1] + dbox[2])
        } else {
            dmax <- max(dbox)
            upperLeft <- c(upper_right[2] + dmax, lower_left[1] - dmax)
            lowerRight <- c(lower_left[2] - dmax, upper_right[1] + dmax)
        }

        osm_backdrop <- try(openmap(
            upperLeft = upperLeft,
            lowerRight = lowerRight,
            type = backdrop,
            zoom = zoom
        ))
        if (is(osm_backdrop, "try-error")) {
            stop("Download failed")
        }
        # Reproject it
        osm_backdrop_reproj <- openproj(osm_backdrop, proj4string(CRS("+init=EPSG:4326")))

        dir.create(dirname(osm_cachefile), recursive = TRUE)
        saveRDS(osm_backdrop_reproj, file = osm_cachefile)
        print("Saved Open Street Map tiles.")
    } else {
        # Read from cache
        osm_backdrop_reproj <- readRDS(osm_cachefile)
    }
    return(osm_backdrop_reproj)
}
