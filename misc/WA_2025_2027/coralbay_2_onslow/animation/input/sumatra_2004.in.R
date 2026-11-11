
# Define Sumatra 2004 scenario animation shots inputs as data frame
library(ncdf4)
library(tibble)
library(dplyr)

scenario <- list(
  multidomain_dir = "/g/data/w85/tsunami/MODELS/inundation/WA_tsunami_inundation_DFES/coralbay_2_onslow/swals/OUTPUTS/historic/sumatra2004_FujiSatake2007_timevarying_movie-full-ambient_sea_level_0/RUN_20260818_105632177",
  intro_title = "2004 Indian Ocean Tsunami",
  title = "2004 Indian Ocean Tsunami",
  dir = "sumatra_2004",
  info_text = "The Indian Ocean tsunami was generated on 26/12/2004 by a Magnitude 9.2 earthquake near Sumatra, Indonesia.\n\n
  Geoscience Australia modelled the tsunami from the earthquake source to the coast of Western Australia.",
  info_text_2 = "",
  info_text_3 = "",
  fps = 24,
  zlim = c(-0.5, 0.5)
)
stopifnot(dir.exists(scenario$multidomain_dir))
scenario$sample_nc_file <- file.path(scenario$multidomain_dir, "RUN_ID00000000010000000001_00018_20260818_105632.306", "Grid_output_ID00000000010000000001.nc")

shots <- rbind(
    tibble(dir = "intro", start = 1, end = 160, file_pattern = "intro_%03d.png"),
    tibble(dir = "indian_ocean", start = 1, end = 120*4, file_pattern = "frame_%03d.png"),
    tibble(dir = "nw_australia", start = 120*4+1, end = 120*7, file_pattern = "frame_%03d.png"),
    tibble(dir = "cape_range", start = 120*4.5+1, end = 120*7, file_pattern = "frame_%03d.png"),
    tibble(dir = "onslow", start = 120*5+1, end = 120*8, file_pattern = "frame_%03d.png"),
    tibble(dir = "coral_bay", start = 120*4.5+1, end = 120*8, file_pattern = "frame_%03d.png"),
    tibble(dir = "midwest_gascoyne", start = 120*5+1, end = 120*8, file_pattern = "frame_%03d.png"),
    tibble(dir = "indian_ocean", start = 120*8+1, end = 120*10+1, file_pattern = "frame_%03d.png")
    #,
    # tibble(dir = "../outro", start = 1, end = 192, file_pattern = "outro_%03d.png")
) %>%
    # see how long they'll be with their framerate
    mutate(sec = (end - start) / scenario$fps)


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
  OUTPUT_TIMES <- get_output_times(scenario$sample_nc_file)
  write.table(OUTPUT_TIMES, file_path, row.names = FALSE, col.names = "seconds")
}

scenario$time <- list(seconds = OUTPUT_TIMES$seconds)
scenario$time_index <- seq_along(scenario$time$seconds)
