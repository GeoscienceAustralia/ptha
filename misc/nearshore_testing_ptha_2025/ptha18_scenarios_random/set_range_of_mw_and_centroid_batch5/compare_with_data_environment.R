#
# Find scenarios with magntiude and centroid location within specified ranges
#

library(rptha)

#' For a given source-zone, read the event metadata and return an environment
#' containing functions that can be used to select events with magnitude and centroid
#' location within specified ranges.
#'
#' @param tsunami_events_directory The PTHA18 directory containing the tsunami events
#' of interest. 
#' @return The function environment is returned, and this contains all the event data
#' as well as functions to do the event selection (documented inside the following).
#' @example
#' tsunami_events_directory = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/puysegur2/TSUNAMI_EVENTS/'
#' puysegur2_env = make_scenario_selection_environment(tsunami_events_directory)
#' 
#' # Get a few types of constant rigidity scenarios, with magnitude and location
#' # roughly matching the 2009 Puysegur earthquake.
#'
#' # Fixed-area uniform-slip
#' puysegur2_uniform_slip_scenarios = puysegur2_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
#'     7.65, 7.85, c(166.56, -45.76), c(166.56, -45.76))
#' # Heterogeneous-slip
#' puysegur2_heterogeneous_slip_scenarios = puysegur2_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
#'     7.65, 7.85, c(166.56, -45.76), c(166.56, -45.76), use_stochastic_slip_runs=TRUE)
#' # Variable-area uniform-slip
#' puysegur2_variable_area_uniform_slip_scenarios = puysegur2_env$find_events_with_magnitude_and_alongstrike_centroid_constraints(
#'     7.65, 7.85, c(166.56, -45.76), c(166.56, -45.76), use_variable_uniform_slip_runs=TRUE)
#'
make_scenario_selection_environment<-function(tsunami_events_directory){

    #
    # Data for this source zone
    #
    #source_name = basename(dirname(getwd()))
    source_name = basename(dirname(tsunami_events_directory))
    # NetCDF file with uniform slip earthquake events
    earthquake_events = read_table_from_netcdf(
        paste0(tsunami_events_directory, '/all_uniform_slip_earthquake_events_', 
               source_name, '.nc'))
    # NetCDF file with stochastic slip earthquake events
    earthquake_events_stochastic = read_table_from_netcdf(
        paste0(tsunami_events_directory, '/all_stochastic_slip_earthquake_events_', 
               source_name, '.nc'))
    # NetCDF file with variable-uniform slip earthquake events
    earthquake_events_variable_uniform = try(read_table_from_netcdf(
        paste0(tsunami_events_directory, '/all_variable_uniform_slip_earthquake_events_', 
               source_name, '.nc')))
    # NetCDF file with unit-source statistics
    unit_source_statistics = read_table_from_netcdf(
        paste0(tsunami_events_directory, '/unit_source_statistics_', 
               source_name, '.nc'))

    # Shapefile with unit-source grid geometry
    unit_source_geometry = readOGR(
        dsn=paste0(dirname(dirname(dirname(unit_source_statistics$initial_condition_file[1]))),
            '/unit_source_grid'),
        layer=source_name)


    #'
    #' Find earthquake events near a point close to a given magnitude
    #'
    find_events_with_magnitude_and_alongstrike_centroid_constraints<-function(
        event_magnitude_lower,
        event_magnitude_upper,
        event_point1,
        event_point2,
        use_stochastic_slip_runs = FALSE,
        use_variable_uniform_slip_runs = FALSE,
        fixed_mu = TRUE){

        # Select events based on the input slip type
        if(use_stochastic_slip_runs | use_variable_uniform_slip_runs){

            if(use_stochastic_slip_runs & use_variable_uniform_slip_runs){
                stop('Cannot have TRUE values for BOTH use_stochastic_slip_runs AND use_variable_uniform_slip_runs')
            }
            if(use_stochastic_slip_runs){
                events = earthquake_events_stochastic
            }else{
                events = earthquake_events_variable_uniform
            }

        }else{

            events = earthquake_events

        }

        # Filter by magnitude
        if(fixed_mu){
            k = which(events$Mw >= event_magnitude_lower & 
                      events$Mw <= event_magnitude_upper)
        }else{
            k = which(events$variable_mu_Mw >= event_magnitude_lower & 
                      events$variable_mu_Mw <= event_magnitude_upper)
        }
        if(length(k) == 0) stop('No events in specified magnitude range')
        events_with_Mw = events[k,]

        #
        # Filter by slip location.
        #

        # First get the unit-source alongstrike indices associated with
        # point1 and point2
        i1 = find_unit_source_index_containing_point(event_point1, unit_source_geometry, unit_source_statistics)
        i2 = find_unit_source_index_containing_point(event_point2, unit_source_geometry, unit_source_statistics)
        alongstrike_lower_ind = min(unit_source_statistics$alongstrike_number[c(i1, i2)])
        alongstrike_upper_ind = max(unit_source_statistics$alongstrike_number[c(i1, i2)])

        # Next get the alongstrike indices of each event centroid
        centroid_alongstrike_ind = rep(NA, nrow(events_with_Mw))
        for(i in 1:nrow(events_with_Mw)){
            ai = get_event_slip_weighted_centroid(events_with_Mw[i,], unit_source_statistics, as_subfault_number=TRUE)
            centroid_alongstrike_ind[i] = ai$alongstrike_number
        }

        # To the selection
        k = which(centroid_alongstrike_ind >= alongstrike_lower_ind & 
                  centroid_alongstrike_ind <= alongstrike_upper_ind)
        if(length(k) == 0) stop('No centroids within target range')
        events_with_Mw = events_with_Mw[k,]

        return(events_with_Mw)
    }

    return(environment())
}
