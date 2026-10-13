#' Copy frames into frames directory with fresh counting then
#' combine into a movie with ffmpeg
library(tibble)
system("module load ffmpeg")

# load camera shots data frame
source("input/puysegur_84.in.R")

frame_dir <- file.path(scenario$dir, "se_australia")

# merge shots into movie
cmd <- paste0(
  "ffmpeg -framerate ", scenario$fps,
  " -i ", frame_dir, "/frame_%03d.png",
  " -c:v libx264 -pix_fmt yuv420p ",
  file.path(scenario$dir, paste0(basename(scenario$dir), ".mp4"))
)

cat(cmd, "\n")

system(cmd)

