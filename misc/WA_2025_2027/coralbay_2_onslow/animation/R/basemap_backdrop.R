library(basemaps)
library(terra)

# Get the OSM data, with caching to speed it up
get_osm_backdrop_reproj <- function(lower_left, upper_right, dbox, backdrop, map_service = "esri") {
    osm_cachefile <- paste0(
        "osm_tile_cache/",
        upper_right[1], "_",
        upper_right[2], "_",
        lower_left[1], "_",
        lower_left[2], "_",
        map_service, "-",
        backdrop,
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
        dmax <- max(dbox)
        extent <- ext(c(
            lower_left[1] - dmax,
            upper_right[1] + dmax,
            lower_left[2] - dmax,
            upper_right[2] + dmax
        ))
        bbox_vect <- vect(extent, crs = "+proj=longlat +datum=WGS84 +lon_wrap=180")

        osm_backdrop <- try(basemap_terra(
            ext = bbox_vect,
            map_service = map_service,
            map_type = backdrop
        ))
        if (is(osm_backdrop, "try-error")) {
            stop("Download failed")
        }
        # Reproject it
        osm_backdrop_reproj <- project(osm_backdrop, "+proj=longlat +datum=WGS84 +lon_wrap=180")

        dir.create(dirname(osm_cachefile), recursive = TRUE)
        saveRDS(osm_backdrop_reproj, file = osm_cachefile)
        print("Saved Imagery tiles.")
    } else {
        # Read from cache
        osm_backdrop_reproj <- readRDS(osm_cachefile)
    }
    return(osm_backdrop_reproj)
}
