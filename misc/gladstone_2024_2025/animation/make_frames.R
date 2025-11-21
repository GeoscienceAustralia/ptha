#' Make a tsunami animation
#' total frame width=1920, height=1080
library(parallel)

source("R/make_image_2d.R")
source("R/osm_backdrop.R")
source("R/combine_frame.R")
source("R/make_ic.R")
source("input/custom_colours.R")

# input scenario data frame called scenario
# source("input/kamchatka.in.R")
# source("input/kermadec_94.in.R")
source("input/solomon_90.in.R")

# Data frame of scenes to zoom into
source("input/scenes.in.R")

# e.g. scene_dir <- "sw_pacific_full"
args <- commandArgs(trailingOnly = TRUE)
scene_dir <- args[1]
scene <- filter(scenes, dir == scene_dir)
print(scene$title)

scene$dir_write <- file.path(scenario$dir, scene$dir)
dir.create(scene$dir_write, recursive = TRUE, showWarnings = FALSE)

# make larger images so domain lines don't render
image_scale <- 1
if (scene$dir == "sw_pacific_full") image_scale <- 2

plot_config <- list(
  cols = ga_ocean_cols(n = 1024, alpha = 1.0),
  zlim = scenario$zlim,
  width = 1300 * image_scale,
  height = 950 * image_scale,
  add_scale_bar = FALSE,
  add_legend = TRUE,
  cex.axis = 2 * image_scale,
  legend.axis = 2,
  legend_line = 6,
  par = list(
    cex.axis = 2 * image_scale,
    cex.lab = 2 * image_scale,
    mgp = c(3, 3, 0),
    oma = c(0, 2, 0, 4)
  )
)

# only make required frames for shots with this scene
this_shots <- left_join(scene, shots, by = "dir")
if (nrow(this_shots) == 0) {
  stop(paste0("No frames required for requested dir: ", scene$dir))
}
time_index <- unique(unlist(mapply(seq, this_shots$start, this_shots$end)))
if (!all(time_index %in% scenario$time_index)) {
  stop("The shots requested a time index not in netCDF file.")
}
scenario$time_index <- time_index

# combine lower left and upper right pairs
ur <- c(scene$ur_x, scene$ur_y)
ll <- c(scene$ll_x, scene$ll_y)

# domains to plot
if (is.null(unlist(scene$domains))) {
  site_grids <- NULL
} else {
  domain_nc_pattern <- paste0(
    scenario$multidomain_dir, "/RUN*0000", unlist(scene$domains), "_*/Grid*.nc"
  )
  site_grids <- unlist(lapply(domain_nc_pattern, Sys.glob))
}

# get backdrop once (download once with internet and cache)
osm_backdrop_reproj <- get_osm_backdrop_reproj(
  ll,
  ur,
  dbox = unlist(scene$dbox),
  backdrop = "esri-imagery",
  zoom = scene$osm_zoom
)

# don't plot if crossing anti-meridian
if (scene$dir == "pacific") osm_backdrop_reproj <- NULL

if (grepl("login", system("hostname", intern = TRUE))) {
  stop("Don't proceed if on a login node.")
}

# make images
image_result <- mclapply(
  scenario$time_index,
  make_image_2d,
  dir_save = scene$dir_write,
  ll = ll,
  ur = ur,
  multidomain_dir = scenario$multidomain_dir,
  backdrop = osm_backdrop_reproj,
  plot_config = plot_config,
  cities = scene$cities[[1]],
  site_grids = site_grids,
  mc.cores = detectCores(logical = TRUE)
)

#  add header and side panel
print("Merging header, side panels and image")
pb <- txtProgressBar(min = 1, max = length(scenario$time_index), style = 2)
for (i in scenario$time_index) {
  combine_frame(i, scenario, scene)
  setTxtProgressBar(pb, i - scenario$time_index[1])
}

print(scene$title)
