#' Kamchatka movie inputs
#' 
#' Two tibbles expected: scenario and shots as described below
#' These are used by make_intro, make_frames.R and merge_frames.R
library(ncdf4)
library(tibble)
library(dplyr)

#' Details about the scenario
scenario <- list(
  multidomain_dir = "../../swals/OUTPUTS/kamchatka2025_42699_movie-full-ambient_sea_level_0/RUN_20250811_095631327",
  intro_title = "Kamchatka 2025 Tsunami Simulation",
  title = "Kamchatka 2025 Tsunami Simulation",
  dir = "kamchatka",
  info_text = "On 29/07/2025 23:24:52 (UTC) a Mw 8.8 earthquake occured 25 km deep off the east coast of Kamchatka. It caused a tsunami throughout the Pacific Ocean. This simulation recreates the tsunami using a similar scenario from the PTHA18 database as it progresses towards Queensland.",
  info_text_2 = "The impacts in the region were negligable and no onshore inundation was recorded. However, tide gauges observed deviations in the sea level which are well matched by the simulations.",
  info_text_3 = "",
  zlim = c(-0.15, 0.15),
  fps = 12
)

#' Shots tibble (df) - scene order and duration
#'   
#' Play selected scenes from the scenes.in.R for these durations
#' dir: the scene's directory
#' start: the index of first frame to play (indexed by the netCDF file)
#' end: the final frame to play (indexed by the netCDF file)
#' file_pattern: the sprintf format to read files inside the dir. e.g. made by
#' make_intro.R and R/make_image_2d.R
shots <- rbind(
    tibble(dir = "intro", start = 1, end = 152, file_pattern = "intro_%03d.png"),
    tibble(dir = "pacific", start = 1, end = 433, file_pattern = "frame_%03d.png"),
    tibble(dir = "coral_sea", start = 120, end = 240, file_pattern = "frame_%03d.png"),
    tibble(dir = "capricornia", start = 120, end = 433, file_pattern = "frame_%03d.png"),
    tibble(dir = "../outro", start = 1, end = 96, file_pattern = "outro_%03d.png")
) %>%
    # see how long they'll be with the framerate
    mutate(sec = (end - start) / scenario$fps)


nc_file <- file.path(scenario$multidomain_dir, "RUN_ID00000000010000000001_00001_20250811_095631.418", "Grid_output_ID00000000010000000001.nc")

# read times saved in the netcdf output.
# Beware, it may have some floating point issues
get_output_times <- function(nc_file) {
  fid <- nc_open(nc_file)
  OUTPUT_TIMES <- ncvar_get(fid, "time")
  # round to the nearest microsecond to avoid floating point issues
  OUTPUT_TIMES <- round(OUTPUT_TIMES, 3)
  nc_close(fid)
  return(OUTPUT_TIMES)
}

file_path = file.path(scenario$dir, "times.csv")
dir.create(scenario$dir)
if (file.exists(file_path)) {
  OUTPUT_TIMES <- read.csv(file_path)
} else {
  OUTPUT_TIMES <- get_output_times(nc_file)
  write.table(OUTPUT_TIMES, file_path, row.names = FALSE, col.names = "seconds")
}

scenario$time <- list(seconds = OUTPUT_TIMES$seconds)
scenario$time_index <- seq_along(scenario$time$seconds)
