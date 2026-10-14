
# Define Sumatra 2004 scenario animation shots inputs as data frame
library(ncdf4)
library(tibble)
library(dplyr)

scenario <- list(
  multidomain_dir = "/g/data/w85/tsunami/MODELS/misc/2026_Puysegur_scenario/OUTPUTS/puysegur2_row_003046_Mw_84_HS",
  title = "Puysegur trench scenario",
  dir = "puysegur_84",
  fps = 24,
  zlim = c(-0.5, 0.5)
)
stopifnot(dir.exists(scenario$multidomain_dir))
scenario$sample_nc_file <- file.path(scenario$multidomain_dir, "RUN_ID00000000000000000001_00001_20261013_112752.767", "Grid_output_ID00000000000000000001.nc")

shots <- rbind(
    tibble(dir = "se_australia", start = 1, end = 961, file_pattern = "intro_%03d.png")
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
