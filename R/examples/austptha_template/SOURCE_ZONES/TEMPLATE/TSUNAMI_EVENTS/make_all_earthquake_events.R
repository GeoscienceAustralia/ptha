library(rptha)
# Get code for summing unit sources (only required for matching tide gauge
# files and unit-sources)
sum_sources = new.env()
source('sum_tsunami_unit_sources.R', local=sum_sources)
# Get key information on tsunami unit sources
tsunami_unit_source_config = new.env()
source('../TSUNAMI_UNIT_SOURCE/config.R', local=tsunami_unit_source_config)

# Make events with magnitudes ranging from Mw_min to Mw_max, separated by dMw
Mw_min = 7.2
Mw_max = 9.8
dMw = 0.1


# Read the discretized sources geometry
discretized_sources_file = normalizePath('../EQ_SOURCE/all_discretized_sources.RDS')
all_discretized_sources = readRDS(discretized_sources_file)

source_zone_name=names(all_discretized_sources)[1]

# For each unit source, make some summary statistics. These are required to
# construct the earthquake events.
unit_source_statistics = discretized_source_approximate_summary_statistics(
    all_discretized_sources[[1]])

#source_tif_names = paste0(
#    normalizePath(paste0('../EQ_SOURCE/Unit_source_data/', source_zone_name)),
#    '/', source_zone_name, '_', 
#    unit_source_statistics$downdip_number, '_',
#    unit_source_statistics$alongstrike_number, 
#    '.tif')

# Append the filenames of the unit source initial conditions to the
# unit_source_statistics
source_tif_names = tsunami_unit_source_config$initial_condition_files
if(!all(file.exists(source_tif_names))){
    stop('Could not find all unit_source tifs')
}
# Make sure the order is the same as in unit_source_statistics
source_tif_match = paste0(source_zone_name, '_', unit_source_statistics$downdip_number, 
    '_', unit_source_statistics$alongstrike_number, '.tif')
mtch = match(source_tif_match, basename(source_tif_names))
if( (any(is.na(mtch)) | (length(source_tif_match) != length(source_tif_names))){
    stop('Problem re-ordering unit-source initial condition files')
}
source_tif_names = source_tif_names[source_tif_match]

unit_source_statistics$initial_condition_file = source_tif_names

# Append the filenames of the unit-source tsunami gauges to the
# unit_source_statistics
source_gauge_names = Sys.glob(
    paste0(tsunami_unit_source_config$all_runs_output_basedir, '/'
        tsunami_unit_source_config$all_runs_dir, '/*/*/Gauge*.nc'))

source_gauge_names = sum_sources$sort_tide_gauge_files_by_unit_source_table(
    unit_source_statistics, source_gauge_names)
unit_source_statistics$tide_gauge_file = source_gauge_names

# Compute all the earthquake events, using the global PTHA paper methodology
all_eq_events = get_all_earthquake_events(
    unit_source_statistics = unit_source_statistics,
    Mmin=Mw_min, 
    Mmax=Mw_max, 
    dMw=dMw, 
    source_zone_name = source_zone_name)
# Add an annual rate variable that we can change later
all_eq_events$rate_annual = rep(-1, length(all_eq_events[,1]))

# Make the plot
pdf(paste0('event_size_scaling_', source_zone_name, '.pdf'), 
    width=10, height=10)
plot_earthquake_event_properties(all_eq_events)
dev.off()

# Save key unit source statistics to csv
write.csv(unit_source_statistics,
    paste0('unit_source_statistics_', source_zone_name, '.csv'),
    row.names=FALSE, quote=FALSE)
# Also save unit source statistics to netcdf
unit_source_statistics_file = paste0(getwd(), '/unit_source_statistics_', source_zone_name, '.nc')
uss_attr = list('discretized_sources_file' = discretized_sources_file,
    'parent_script_name' = parent_script_name())
write_table_to_netcdf(unit_source_statistics, 
    basename(unit_source_statistics_file),
    global_attributes_list = uss_attr,
    add_session_info_attribute=TRUE)

# Save earthquake events to csv
write.csv(all_eq_events,
    paste0('all_eq_events_', source_zone_name, '.csv'),
    row.names=FALSE, quote=FALSE)
# Also save earthquake events to netcdf
eq_events_attr = c(
    list('unit_source_statistics_file' = unit_source_statistics_file,
        'slip_type' = 'uniform slip'), 
    uss_attr)
write_table_to_netcdf(all_eq_events,
    paste0('all_eq_events_', source_zone_name, '.nc'),
    global_attributes_list = eq_events_attr,
    add_session_info_attribute=TRUE)

