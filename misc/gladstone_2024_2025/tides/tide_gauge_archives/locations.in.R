
#' One named list for each location containing the data and the threshold
#'
#' @param threshold_tidal_plane_MSL: Enter the value of HAT used in the model at each location from the official semidiurnal tidal plane model 2024: https://www.msq.qld.gov.au/tmr_site_msq/_/media/tmronline/msqinternet/msqfiles/home/tides/2024-semidiurnal-tidal-planes.pdf
#' This is the value of HAT minus the value of MSL (both reported in LAT).
#' @param msl_in_LAT: The mean sea level (m above Queensland Port Datum) comes from https://www.msq.qld.gov.au/tmr_site_msq/_/media/tmronline/msqinternet/msqfiles/home/tides/msl-predictions-2024.pdf
#' Note that Queensland Port Datum (called LAT) is used in both the MSL and the station datum. This is also reported in the semidiurnal tidal plane model 2024, but with lower precision.
#' @param threshold_TPXO9_MSL: The max value from the TPXO9 model at the location. Inspect the raster of the TPXO9 model at the location to find the max value.
locations <- list(
    auckland_point = list(
        name = "Auckland Point",
        file_write_name_stem = "auckland_point",
        files = Sys.glob("data/auckland_point/*.csv"),
        metadata_file = Sys.glob("data/auckland_point/*.ashx"),
        msl_in_LAT = 2.406,
        threshold_tidal_plane_MSL = NA,
        threshold_TPXO9_MSL =  2.22
    ),
    fishermans_landing = list(
        name = "Fishermans Landing",
        file_write_name_stem = "fishermans_landing",
        files = Sys.glob("data/fishermans_landing/*.csv"),
        metadata_file = Sys.glob("data/fishermans_landing/*txt"),
        msl_in_LAT = NA,
        threshold_tidal_plane_MSL = 5.16-2.47,
        threshold_TPXO9_MSL = 2.22,
        bounds_AHD = c(-3.0, Inf)  # On inspection, one point it -5.0m AHD, which is clearly an error.
    ),
    port_alma = list(
        name = "Port Alma",
        file_write_name_stem = "port_alma",
        files = Sys.glob("data/port_alma/*.csv"),
        metadata_file = Sys.glob("data/port_alma/*.ashx"),
        msl_in_LAT = 2.966,
        threshold_tidal_plane_MSL = 5.96-2.95,
        threshold_TPXO9_MSL = 2.674
    ),
    rosslyn_bay = list(
        name = "Rosslyn Bay",
        file_write_name_stem = "rosslyn_bay",
        files = Sys.glob("data/rosslyn_bay/*.csv"),
        metadata_file = Sys.glob("data/rosslyn_bay/*.ashx"),
        msl_in_LAT = 2.494,
        threshold_tidal_plane_MSL = 5.21-2.48,
        threshold_TPXO9_MSL = 2.556
    ),
    south_trees = list(
        name = "South Trees",
        file_write_name_stem = "south_trees",
        files = Sys.glob("data/south_trees/*.csv"),
        metadata_file = Sys.glob("data/south_trees/*txt"),
        msl_in_LAT = NA,  # No data published
        threshold_tidal_plane_MSL = 4.67-2.26,
        threshold_TPXO9_MSL =  2.15
    ),
    mackay = list(
        name = "Mackay",
        file_write_name_stem = "mackay",
        files = Sys.glob("data/mackay/*.csv"),
        metadata_file = Sys.glob("data/mackay/*.ashx"),
        msl_in_LAT = 3.082,
        threshold_tidal_plane_MSL = 6.62 - 3.07,
        threshold_TPXO9_MSL = 3.11
    ),
    bundaberg = list(
        name = "Bundaberg",
        file_write_name_stem = "bundaberg",
        files = Sys.glob("data/bundaberg/*.csv"),
        metadata_file = Sys.glob("data/bundaberg/*.ashx"),
        msl_in_LAT = 1.794,
        threshold_tidal_plane_MSL = 3.68-1.78,
        threshold_TPXO9_MSL = 1.65
    )
)
