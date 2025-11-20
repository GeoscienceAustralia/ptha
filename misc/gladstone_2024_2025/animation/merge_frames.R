#' Copy frames into frames directory with fresh counting then
#' combine into a movie with ffmpeg
library(tibble)


frame_rate <- 24

# load camera shots data frame
# source("input/kermadec_94.in.R")
source("input/solomon_90.in.R")


frame_dir <- file.path(scenario$dir, "frames")
dir.create(frame_dir)

# consistent frame count for output
i <- 0

print(paste0("Copying all png files into ", frame_dir, " with single frame count."))
for (scene_idx in seq_len(nrow(shots))) {
    scene <- shots[scene_idx, ]

    for (t in seq(scene$start, scene$end)) {
        img <- file.path(scenario$dir, scene$dir, sprintf(scene$file_pattern, t))
        if (!file.exists(img)) stop(img)
        new_name <- file.path(frame_dir, sprintf("frame_%05d.png", i))

        file.copy(img, new_name)
        i <- i + 1
    }
}

# merge shots into movie
system(
    paste0(
        "ffmpeg -framerate ", frame_rate, " -i ", frame_dir, "/frame_%05d.png -c:v libx264 -pix_fmt yuv420p ", scenario$dir, "/", scenario$dir, ".mp4"
    )
)
