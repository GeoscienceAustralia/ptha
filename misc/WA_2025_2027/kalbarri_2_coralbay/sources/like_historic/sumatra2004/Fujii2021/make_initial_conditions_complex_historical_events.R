library(geosphere)
#
# Sumatra 2004, Fujii et al 2021:
#
source_data_file = 'rupture_summary.csv'
source_data = read.csv(source_data_file)
source_tifs = normalizePath(paste0('Fuji21_andaman2004_unit_sources_', 1:nrow(source_data), '.tif'))

# The inversion we use assumes a rupture propagation speed of 1.3km/s
# from the hypocentre and a 1min rise time, which is the best fit solution
# from the paper. 
# Here we use GCMT to define the epicentre, and distance-from-epicentre to compute the start time.
start_time = distHaversine(source_data[,1:2], c(95.78, 3.3))/1300
rise_time = 60 # 1 min rise-time
end_time = start_time + rise_time

# Realistic case -- here the slip is included in the unit sources
swals_forcing_data = data.frame(source_tifs=source_tifs, slip=source_data$slip*0 + 1, start_time=start_time, end_time=end_time)
write.csv(swals_forcing_data, file='Fujii21_time_varying_forcing_realistic.csv', row.names=FALSE)

# Case where rupture starts immediately (for testing)
swals_forcing_data = data.frame(source_tifs=source_tifs, slip=source_data$slip*0 + 1, start_time=0*start_time, end_time=0*end_time + 1.0e-03)
write.csv(swals_forcing_data, file='Fujii21_time_varying_forcing_instantaneous.csv', row.names=FALSE)

