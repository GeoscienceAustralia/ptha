#' Solomon Mw 9.0 movie inputs
#' 
#' Two tibbles expected: scenario and shots as described below
#' These are used by make_intro, make_frames.R and merge_frames.R
library(ncdf4)
library(tibble)
library(dplyr)

scenario <- list(
  multidomain_dir = "../../swals/OUTPUTS/movie_solomon2_0015016_Mw_90-full-ambient_sea_level_0/RUN_20251111_103624872",
  intro_title = "Gladstone Tsunami Hazards",
  title = "Solomon Mw 9.0 Tsunami Simulation",
  dir = "solomon_90",
  info_text = "Geoscience Australia modelled hundreds of earthquake tsunami scenarios to assess the hazard in and around Gladstone.",
  info_text_2 = "This simulation shows one example - a hypothetical tsunami from a magnitude 9.0 earthquake near the Solomon Islands. It is a large but plausible event.",
  info_text_3 = "",
  zlim = c(-2.0, 2.0),
  fps = 24
)

shots <- rbind(
  tibble(dir = "intro", start = 1, end = 456, file_pattern = "intro_%03d.png"),
  tibble(dir = "coral_sea", start = 1, end = 576, file_pattern = "frame_%03d.png"),
  tibble(dir = "capricornia", start = 360, end = 720, file_pattern = "frame_%03d.png"),
  tibble(dir = "agnes_waters", start = 500, end = 740, file_pattern = "frame_%03d.png"),
  tibble(dir = "gladstone_harbour", start = 560, end = 920, file_pattern = "frame_%03d.png"),
  tibble(dir = "boyne_island", start = 600, end = 840, file_pattern = "frame_%03d.png"),
  tibble(dir = "yeppoon_offshore", start = 600, end = 960, file_pattern = "frame_%03d.png"),
  tibble(dir = "capricornia", start = 961, end = 1393, file_pattern = "frame_%03d.png"),
  tibble(dir = "coral_sea", start = 1394, end = 1802, file_pattern = "frame_%03d.png"),
  tibble(dir = "../outro", start = 1, end = 192, file_pattern = "outro_%03d.png")
) %>%
  # see how long they'll be with the framerate
  mutate(sec = (end - start) / scenario$fps)


nc_file <- file.path(scenario$multidomain_dir, "RUN_ID00000000010000000001_00001_20251111_103624.945", "Grid_output_ID00000000010000000001.nc")

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
dir_save <- "solomon_90"
file_path <- file.path(dir_save, "times.csv")
dir.create(dir_save)
if (file.exists(file_path)) {
  OUTPUT_TIMES <- read.csv(file_path)
} else {
  OUTPUT_TIMES <- get_output_times(nc_file)
  write.table(OUTPUT_TIMES, file_path, row.names = FALSE, col.names = "seconds")
}

scenario$time <- list(seconds = OUTPUT_TIMES$seconds)
scenario$time_index <- seq_along(scenario$time$seconds)
