#
# Make netcdf files describing the earthquake events, for the cases of:
#
# 1) Fixed size uniform slip earthquake events
# 2) Stochastic slip, variable size earthquake events
# 3) Uniform slip, variable size earthquake events
#

library(rptha)
#
# Get code for summing unit sources (only required for matching tide gauge
# files and unit-sources)
#
sum_sources = new.env()
source('sum_tsunami_unit_sources.R', local=sum_sources)
#
# Get key information on tsunami unit sources
#
tsunami_unit_source_config = new.env()
source('../TSUNAMI_UNIT_SOURCE/config.R', local=tsunami_unit_source_config, 
    chdir=TRUE)
#
# Get local parameters code
#
config_env = new.env()
source('config.R', local=config_env)

# Make events with magnitudes ranging from Mw_min to Mw_max, separated by dMw
# These values do not have to correspond to Mw ranges that are assigned non-zero 
# probability -- but they should fully contain such possible values
Mw_min = config_env$Mw_min 
Mw_max = config_env$Mw_max 
dMw = config_env$dMw 

# Make at least N stochastic and/or variable_uniform slip events for each
# fixed-size uniform event
number_stochastic_events_for_each_uniform_event = 
    config_env$number_stochastic_events_for_each_uniform_event
# ... but ensure that there are at least M >> N stochastic slip events in each
# magnitude category
minimum_number_stochastic_events_for_each_magnitude = 
    config_env$minimum_number_stochastic_events_for_each_magnitude

#
# Read the discretized sources geometry
#
discretized_sources_file = normalizePath(
    '../EQ_SOURCE/all_discretized_sources.RDS')
all_discretized_sources = readRDS(discretized_sources_file)

source_zone_name=names(all_discretized_sources)[1]

#
# For each unit source, make some summary statistics. These are required to
# construct the earthquake events.
#
unit_source_statistics = discretized_source_approximate_summary_statistics(
    all_discretized_sources[[1]])

#
# Append the filenames of the unit source initial conditions to the
# unit_source_statistics
#
source_tif_names = tsunami_unit_source_config$initial_condition_files
if(!all(file.exists(source_tif_names))){
    stop('Could not find all unit_source tifs')
}
# Make sure the filename order is the same as in unit_source_statistics
source_tif_match = paste0(source_zone_name, '_', unit_source_statistics$downdip_number, 
    '_', unit_source_statistics$alongstrike_number, '.tif')
mtch = match(source_tif_match, basename(source_tif_names))
if( any(is.na(mtch)) | (length(source_tif_match) != length(source_tif_names))){
    stop('Problem re-ordering unit-source initial condition files')
}
source_tif_names = source_tif_names[mtch]

unit_source_statistics$initial_condition_file = source_tif_names

#
# Append the filenames of the unit-source tsunami gauges to the
# unit_source_statistics
#
source_gauge_names = Sys.glob(
    paste0(tsunami_unit_source_config$all_runs_output_base_dir, '/', 
        tsunami_unit_source_config$all_runs_dir, '/*/*/Gauge*.nc'))

source_gauge_names = sum_sources$sort_tide_gauge_files_by_unit_source_table(
    unit_source_statistics, source_gauge_names)
unit_source_statistics$tide_gauge_file = source_gauge_names

#
# Compute all the earthquake events, using the global PTHA paper methodology
#
all_eq_events = get_all_earthquake_events(
    unit_source_statistics = unit_source_statistics,
    Mmin=Mw_min, 
    Mmax=Mw_max, 
    dMw=dMw, 
    source_zone_name = source_zone_name)

# Add an annual rate variable that we can change later, with a value that is
# obviously double-precision [so the netcdf data type is correctly inferred]
all_eq_events$rate_annual = rep(config_env$null_double, 
    length(all_eq_events[,1]))
all_eq_events$rate_annual_upper_ci = rep(config_env$null_double, 
    length(all_eq_events[,1]))
all_eq_events$rate_annual_lower_ci = rep(config_env$null_double, 
    length(all_eq_events[,1]))


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
unit_source_statistics_file = paste0(getwd(), '/unit_source_statistics_', 
    source_zone_name, '.nc')
uss_attr = list('discretized_sources_file' = discretized_sources_file,
    'parent_script_name' = parent_script_name())
write_table_to_netcdf(unit_source_statistics, 
    basename(unit_source_statistics_file),
    global_attributes_list = uss_attr,
    add_session_info_attribute=TRUE)

# Save earthquake events to csv
write.csv(all_eq_events,
    paste0('all_uniform_slip_earthquake_events_', source_zone_name, '.csv'),
    row.names=FALSE, quote=FALSE)
# Also save earthquake events to netcdf
eq_events_attr = c(
    list('unit_source_statistics_file' = unit_source_statistics_file,
        'slip_type' = 'uniform slip'), 
    uss_attr)
all_eq_events_nc_file = paste0('all_uniform_slip_earthquake_events_', 
    source_zone_name, '.nc')
write_table_to_netcdf(
    all_eq_events,
    all_eq_events_nc_file,
    global_attributes_list = eq_events_attr,
    add_session_info_attribute=TRUE)

#
# Stochastic slip events, and uniform-slip events with variable dimensions
#
for(variable_uniform_slip in c(FALSE, TRUE)){

    stochastic_events_store = vector(mode='list', 
        length=length(all_eq_events[,1]))

    for(i in 1:length(all_eq_events[,1])){

        eq_event = all_eq_events[i,]   
     
        # Find a centroid for the earthquake
        usi = get_unit_source_indices_in_event(eq_event)
        # The next steps assume the following condition holds
        stopifnot(all(unit_source_statistics$subfault_number[usi] == usi))

        # Find mean of unit-source locations, accounting for spherical
        # coordinates
        unit_source_locations = cbind(unit_source_statistics$lon_c[usi], 
            unit_source_statistics$lat_c[usi])
        if(length(usi) > 1){
            event_hypocentre = geomean(unit_source_locations, 
                w=rep(1, length(usi)) )
        }else{
            event_hypocentre = unit_source_locations
        }

        event_magnitude = eq_event$Mw

        # Determine the number of stochastic events -- ensuring there are
        # enough events in this magnitude category, and that there are a
        # minimum number of events for each uniform slip event
        number_of_uniform_events_with_same_magnitude = 
            sum(all_eq_events$Mw == event_magnitude)

        number_of_sffm = max(
            number_stochastic_events_for_each_uniform_event, 
            ceiling(minimum_number_stochastic_events_for_each_magnitude / 
                number_of_uniform_events_with_same_magnitude)
            )

        # Make stochastic events
        all_events = sffm_make_events_on_discretized_source(
            unit_source_statistics,    
            target_location = event_hypocentre,
            target_event_mw = event_magnitude,
            num_events = number_of_sffm,
            zero_low_slip_cells_fraction=0.0,
            sourcename = source_zone_name,
            uniform_slip = variable_uniform_slip)

        events_with_Mw = sffm_events_to_table(all_events, 
            slip_significant_figures=
                config_env$stochastic_slip_table_significant_figures)

        # Add additional variables we will need. Use an obviously
        # floating point number to get the netcdf output to be in double
        events_with_Mw$uniform_event_row = i
        events_with_Mw$rate_annual = config_env$null_double
        events_with_Mw$rate_annual_lower_ci = config_env$null_double
        events_with_Mw$rate_annual_upper_ci = config_env$null_double

        stochastic_events_store[[i]] = events_with_Mw
        
    }
    stochastic_events_table = do.call(rbind, stochastic_events_store)

    if(!variable_uniform_slip){
        #
        # Save stochastic events to csv
        #
        write.csv(stochastic_events_table,
            paste0('all_stochastic_slip_earthquake_events_', source_zone_name, 
                '.csv'),
            row.names=FALSE, quote=FALSE)
        # Also save to netcdf
        eq_events_attr = c(
            list('unit_source_statistics_file' = unit_source_statistics_file,
                 'corresponding_uniform_slip_events_file' = all_eq_events_nc_file, 
                'slip_type' = 'stochastic slip'), 
            uss_attr)
        write_table_to_netcdf(stochastic_events_table,
            paste0('all_stochastic_slip_earthquake_events_', source_zone_name, 
                '.nc'),
            global_attributes_list = eq_events_attr,
            add_session_info_attribute=TRUE)

    }else{
        #
        # Save variable_uniform slip events to csv
        #
        write.csv(stochastic_events_table,
            paste0('all_variable_uniform_slip_earthquake_events_', 
                source_zone_name, '.csv'),
            row.names=FALSE, quote=FALSE)
        # Also save to netcdf
        eq_events_attr = c(
            list('unit_source_statistics_file' = unit_source_statistics_file,
                 'corresponding_uniform_slip_events_file' = all_eq_events_nc_file, 
                'slip_type' = 'variable uniform slip'), 
            uss_attr)
        write_table_to_netcdf(stochastic_events_table,
            paste0('all_variable_uniform_slip_earthquake_events_', 
                source_zone_name, '.nc'),
            global_attributes_list = eq_events_attr,
            add_session_info_attribute=TRUE)

    }
}
