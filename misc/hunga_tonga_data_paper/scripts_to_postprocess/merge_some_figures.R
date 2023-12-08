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

#
# Figure with the arrival time of the MSLP residual maxima
#
big_fig = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/lamb_wave_arrival_and_pressure_maxima_arrival.png'))
inset_fit = image_read(paste0(OUTPUT_GRAPHICS_DIR, '/outlier_pressure_time_series.png'))
inset_fit = image_border(inset_fit, "black", "5x5")
combined_fig = image_composite(big_fig, image_scale(inset_fit, "x480"), offset="+970+700")
image_write(combined_fig, path=paste0(OUTPUT_GRAPHICS_DIR, '/combined_lamb_wave_arrival_pressure_maxima_and_ts.png'), format='png')
