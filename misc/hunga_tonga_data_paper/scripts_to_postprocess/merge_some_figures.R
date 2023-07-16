#
# For the paper, we combine a few figures into a single figure
#
library(magick) # Great tutorial here: https://cran.r-project.org/web/packages/magick/vignettes/intro.html
source('global_variables.R')

#
# Figure with locations of MSLP and tide sensors
#
mslp_fig = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/BOM_MSLP_gauge_locations.png'))
tg_fig = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_locations.png'))
regional_setting_fig = image_append(c(mslp_fig, tg_fig), stack=TRUE)
image_write(regional_setting_fig, path=paste0(OUTPUT_GRAPHICS_DIR, '/combined_gauge_locations.png'), format='png')

#
# Figure with tide-gauge maxima in Australia
#
big_fig = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_peak_size.png'))
inset_fig = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/tide_gauge_peak_size_inset.png'))
combined_fig = image_append(c(big_fig, inset_fig))
combined_fig_title = image_annotate(combined_fig, "Maximum high-pass filtered sea level (m), January 15-16", size=100, gravity='north')
image_write(combined_fig_title, path=paste0(OUTPUT_GRAPHICS_DIR, '/combined_tide_gauge_peak_size_fig.png'), format='png')
