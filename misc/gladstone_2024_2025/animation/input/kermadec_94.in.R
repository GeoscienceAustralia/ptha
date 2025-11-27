#' Kermadec-Tonga Mw 9.4 movie inputs
#' 
#' Two tibbles expected: scenario and shots as described below
#' These are used by make_intro, make_frames.R and merge_frames.R

library(ncdf4)
library(tibble)
library(dplyr)

#' Details about the scenario
scenario <- list(
  multidomain_dir = "../../swals/OUTPUTS/movie_kermadectonga2_0043427_Mw_94-full-ambient_sea_level_0/RUN_20251103_111445644",  # directory containing the simulation
  intro_title = "Gladstone Tsunami Hazards",
  title = "Kermadec Mw 9.4 Tsunami Simulation",
  dir = "kermadec_94",  # directory to make and store output frames and movie
  info_text = "Geoscience Australia modelled hundreds of earthquake tsunami scenarios to assess the hazard in and around Gladstone.",
  info_text_2 = "This simulation shows one example - a hypothetical tsunami from a magnitude 9.4 earthquake between New Zealand and Tonga. It is a large but plausible event.",
  info_text_3 = "There is widespread evidence of a large 15th century tsunami in the south-west Pacific, recorded in coastal sediments. A Mw 9.4 earthquake has been suggested as the source.",
  zlim = c(-2., 2.),  # z limit for colourmap of wave heights
  fps = 24  # frames per second
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
  tibble(dir = "intro", start = 1, end = 744, file_pattern = "intro_%03d.png"),
  tibble(dir = "sw_pacific_full", start = 1, end = 720, file_pattern = "frame_%03d.png"),
  tibble(dir = "capricornia", start = 541, end = 720, file_pattern = "frame_%03d.png"),
  tibble(dir = "agnes_waters", start = 721, end = 948, file_pattern = "frame_%03d.png"),
  tibble(dir = "gladstone_harbour", start = 780, end = 1212, file_pattern = "frame_%03d.png"),
  tibble(dir = "yeppoon_offshore", start = 840, end = 1080, file_pattern = "frame_%03d.png"),
  tibble(dir = "boyne_island", start = 1081, end = 1321, file_pattern = "frame_%03d.png"),
  # tibble(dir = "lady_elliot", start = 1322, end = 1562, file_pattern = "frame_%03d.png"),
  tibble(dir = "capricornia", start = 1322, end = 1801, file_pattern = "frame_%03d.png"),
  tibble(dir = "sw_pacific", start = 1802, end = 2162, file_pattern = "frame_%03d.png"),
  tibble(dir = "../outro", start = 1, end = 192, file_pattern = "outro_%03d.png")
) %>%
  # see how long they'll be with the framerate
  mutate(sec = (end - start) / scenario$fps)

# Use any netCDF file from the simulation to extract which times are saved
nc_file <- file.path(scenario$multidomain_dir, "RUN_ID00000000010000000001_00001_20251103_111445.716", "Grid_output_ID00000000010000000001.nc")

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

file_path <- file.path(scenario$dir, "times.csv")
dir.create(scenario$dir)
if (file.exists(file_path)) {
  OUTPUT_TIMES <- read.csv(file_path)
} else {
  OUTPUT_TIMES <- get_output_times(nc_file)
  write.table(OUTPUT_TIMES, file_path, row.names = FALSE, col.names = "seconds")
}

scenario$time <- list(seconds = OUTPUT_TIMES$seconds)
scenario$time_index <- seq_along(scenario$time$seconds)
