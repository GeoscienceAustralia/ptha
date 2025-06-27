library(geosphere)
#
# Sumatra 2004, Fuji and Sakate
#
source_data_file = 'rupture_summary.csv'
source_data = read.csv(source_data_file)
source_tifs = paste0('Fuji_andaman2004_unit_sources_', 1:nrow(source_data), '_KAJIURA_SMOOTHED.tif')
source_tifs = normalizePath(source_tifs)

# The inversion we use assumes a rupture propagation speed of 1km/s
# from the hypocentre and a 3min rise time, which is the best fit solution
# from the paper. 
# Here we use GCMT to define the epicentre (similar but perhaps not identical
# to what Fuji and Satake used -- they don't provide coordinates, would have to
# look up Harvard CMT), and distance-from-epicentre to compute the start time.
start_time = distHaversine(source_data[,1:2], c(95.78, 3.3))/1000
rise_time = 180 # 3 min rise-time
end_time = start_time + rise_time

# Realistic case
swals_forcing_data = data.frame(source_tifs=source_tifs, slip=source_data$slip, start_time=start_time, end_time=end_time)
write.csv(swals_forcing_data, file='FujiSatake2007_time_varying_forcing_realistic.csv', row.names=FALSE)

# Case where rupture starts immediately (for testing)
swals_forcing_data = data.frame(source_tifs=source_tifs, slip=source_data$slip, start_time=0*start_time, end_time=0*end_time + 1.0e-03)
write.csv(swals_forcing_data, file='FujiSatake2007_time_varying_forcing_instantaneous.csv', row.names=FALSE)

