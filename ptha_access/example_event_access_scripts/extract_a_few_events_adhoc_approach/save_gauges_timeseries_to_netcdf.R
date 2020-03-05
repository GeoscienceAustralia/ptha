#
# What needs to change
# -- Need to include a global attribute 'reference point for stage exceedance rates'
# -- Include an attribute with the 'stage exceedance rate curves file'
# -- Include a global attribute with the date?
# -- Each event needs a variable 'stage exceedance rate at reference point'
# -- Need to extend e.g. max_stage to have dimensions [ntime, nevent, nstation],
#    remove other flow variables, and add uh, vh
# -- Can only store results from a single source-zone, assuming we put the unit-source information
#    inside the file. ALTERNATIVE: Remove unit-source statistics, but give them the
#    unit-source-statistics netcdf files AS WELL AS the associated shapefiles separately.
#    Should also give them the tsunami_stage_exceedance_rates for the reference point.
# -- Another alternative. Use 'ncks' to extract subsets we require. e.g 'ncks -d mydimension,required_index input_file output_file'


write_tsunami_gauges_to_netcdf<-function(
    file_scenario_name,
    all_eq_events,
    gauge_locations,
    gauge_obs_times,
    gauge_event_STAGE,
    gauge_event_UH,
    gauge_event_VH,
    peak_stage_at_reference_point,
    peak_stage_exceedance_rate_at_reference_point,
    proportion_rate_from_this_source,
    reference_point,
    reference_point_rp_curve_stages,
    reference_point_rp_curve_rates,
    reference_point_rp_curve_rates_lower_ci,
    reference_point_rp_curve_rates_upper_ci,
    desired_exceedance_rate,
    output_file_name
    ){

    library(ncdf4)

    # Station dimension -- make this one unlimited for fast 'single-station'
    # access to the data
    dim_station = ncdim_def(name='station', units='', 
        vals=1:length(gauge_locations[,1]), unlim=TRUE,
        longname='integer index corresponding to the gauge locations')
    
    dim_reference_point_rp_curve = ncdim_def(name='nrp', units='', 
        vals=1:length(reference_point_rp_curve_stages), unlim=FALSE,
        longname='integer index corresponding to each point on the discretized return period curve')

    # Event dimension
    dim_event = ncdim_def(name='event', units='', 
        vals=1:length(all_eq_events[,1]), unlim=FALSE,
        longname='integer index for each tsunami event')

    dim_time = ncdim_def(name='time', units='s',
        vals=gauge_obs_times, unlim=FALSE,
        longname='time at which gauges are observed')



    # Figure out a dimension size for a string
    charlen_sourcename = max(nchar(all_eq_events$sourcename))
    charlen_event_index_string = max(nchar(all_eq_events$event_index_string))
    charlen_event_slip_string = max(nchar(all_eq_events$event_slip_string))

    dim_char_size = max(c(charlen_sourcename, charlen_event_index_string, charlen_event_slip_string))

    dim_nchar = ncdim_def(name='max_nchar', units='', vals=1:dim_char_size, unlim=FALSE,
        longname='integer index corresponding to the maximum number of characters in strings')

    #
    # Create netcdf variable for gauges
    #
    gauge_lon_v = ncvar_def(name='lon', units='degrees_east', 
        dim=list(dim_station), missval=NA, longname='station_longitude', 
        prec='float')
    gauge_lat_v = ncvar_def(name='lat', units='degrees_north', 
        dim=list(dim_station), missval=NA, longname='station_latitude', 
        prec='float')
    gauge_elev_v = ncvar_def(name='elev', units='m', dim=list(dim_station), 
        missval=NA, longname='station_ground_elevation_above_mean_sea_level', 
        prec='float')
    gauge_id_v = ncvar_def(name='gaugeID', units='', dim=list(dim_station), 
        missval=NA, longname='real_ID_for_each_station', prec='float')

    # Keep a list of variables for when we make the ncdf
    all_nc_var = list(gauge_lon_v, gauge_lat_v, gauge_elev_v, gauge_id_v)

    #
    # Create a netcdf variable for the gauge summary statistics
    #
    gauge_event_stage_v = ncvar_def(name='stage', units='m', 
        dim=list(dim_time, dim_event, dim_station), missval=NA, 
        longname='stage (i.e. waterlevel) above MSL',
        prec='float')
    gauge_event_UH_v= ncvar_def(name='uh', units='m', 
        dim=list(dim_time, dim_event, dim_station), missval=NA, 
        longname='(easterly_velocity)x(depth)',
        prec='float')
    gauge_event_VH_v= ncvar_def(name='vh', units='m', 
        dim=list(dim_time, dim_event, dim_station), missval=NA, 
        longname='(northerly_velocity)x(depth)',
        prec='float')

    # Keep a list of variables for when we make the netcdf 
    all_nc_var = c(all_nc_var, list(gauge_event_stage_v, gauge_event_UH_v, gauge_event_VH_v))


    #
    # Create a netcdf variable for the event summary statistics
    #
    # These statistics should vary depending on whether we do stochastic_slip or
    # uniform_slip. But here we assume stochastic
    #
    # Assumes stochastic slip case
    #
    event_Mw_v = ncvar_def(name='event_Mw', units='', dim=list(dim_event),
        missval=NA, longname='moment_magnitude_of_earthquake_event', 
        prec='float')
    event_target_lon_v = ncvar_def(name='event_target_lon', 
        units='degrees_east', dim=list(dim_event),
        missval=NA, longname='longitude_near(ish)_earthquake_event', 
        prec='float')
    event_target_lat_v = ncvar_def(name='event_target_lat', 
        units='degrees_north', dim=list(dim_event),
        missval=NA, longname='latitude_near(ish)_earthquake_event', prec='float')

    event_peak_slip_downdip_ind_v = ncvar_def(
        name='event_peak_slip_downdip_index', units='', 
        dim=list(dim_event), missval=NULL, 
        longname='down-dip_index_of_the_unit_source_with_peak_slip',
        prec='integer')
    event_peak_slip_alongstrike_ind_v = ncvar_def(
        name='event_peak_slip_alongstrike_index', units='', 
        dim=list(dim_event), missval=NULL, 
        longname='along-strike_index_of_the_unit_source_with_peak_slip',
        prec='integer')
    event_uniform_event_row_v = ncvar_def(
        name='event_uniform_event_row', units='', 
        dim=list(dim_event), missval=NULL, 
        longname='row_index_of_the_uniform_slip_event_used_to_define_location_and_magnitude',
        prec='integer')
    
    event_sourcename_v = ncvar_def(name='event_sourcename', units='', 
        dim=list(dim_nchar, dim_event),
        missval=NULL, longname='source_zone_name_hosting_the_earthquake_event', 
        prec='char')

    event_index_string_v = ncvar_def(name='event_index_string', units='', 
        dim=list(dim_nchar, dim_event),
        missval=NULL, 
        longname='indices_of_unit_sources_included_in_earthquake_event_with_-_separator', 
        prec='char')
    event_slip_string_v = ncvar_def(name='event_slip_string', units='', 
        dim=list(dim_nchar, dim_event), missval=NULL, 
        longname='slip_(m)_of_unit_sources_included_in_earthquake_event_with_underscore_separator', 
        prec='char')

    event_wave_height_exceedance_rate_at_reference_point_v = ncvar_def(
        name='event_peak_stage_exceedance_rate_at_reference_point', units='events/year',
        dim=list(dim_event),
        missval=NULL,
        longname='exceedance rate of the maximum stage for this event at the reference point. This should be near to (but probably not equal) the desired value. Exact equality is usually not possible with a finite number of scenarios')
    event_peak_stage_at_reference_point_v = ncvar_def(
        name='event_peak_stage_at_reference_point', units='m',
        dim=list(dim_event),
        missval=NULL,
        longname='maximum stage for this event at the reference point. This should be near to (but probably not equal) the desired value. Exact equality is usually not possible with a finite number of scenarios')
    
    event_proportion_rate_from_this_source_v = ncvar_def(
        name='proportion_rate_from_this_source', units='',
        dim=list(dim_event),
        missval=NULL,
        longname='the proportion of "all events with this stage exceedance rate at the return period point" that will originate from the current source-zone. Values close to 1.0 mean that the current source-zone is a likely source of such events, while values close to 0.0 mean the current source zone is relatively unlikely to generate such an event')
    
    event_desired_exceedance_rate_v = ncvar_def(
        name='desired_exceedance_rate', units='',
        dim=list(dim_event),
        missval=NULL,
        longname='the stage exceedance rate value requested by the event selection script. This will be slightly different to the event_peak_stage_exceedance_rate_at_reference_point, because we only have a finite number of scenarios')

    all_nc_var = c(all_nc_var, list(event_Mw_v, event_target_lon_v, 
        event_target_lat_v, event_peak_slip_downdip_ind_v, 
        event_peak_slip_alongstrike_ind_v, event_uniform_event_row_v,
        event_sourcename_v, event_index_string_v, event_slip_string_v,
        event_wave_height_exceedance_rate_at_reference_point_v,
        event_peak_stage_at_reference_point_v,
        event_proportion_rate_from_this_source_v,
        event_desired_exceedance_rate_v))

    # Stage vs rate curve at reference point
    reference_point_rp_curve_stages_v = ncvar_def(name='rp_curve_stage',
        units='m', dim=list(dim_reference_point_rp_curve),
        missval=NULL,
        longname='stages for stage-vs-rate curve at reference point')
    reference_point_rp_curve_rates_v = ncvar_def(name='rp_curve_rate',
        units='m', dim=list(dim_reference_point_rp_curve),
        missval=NULL,
        longname='rates for stage-vs-rate curve at reference point')
    reference_point_rp_curve_rates_lower_ci_v = ncvar_def(name='rp_curve_rate_lower_ci',
        units='m', dim=list(dim_reference_point_rp_curve),
        missval=NULL,
        longname='rates (lower 95% credible interval) for stage-vs-rate curve at reference point')
    reference_point_rp_curve_rates_upper_ci_v = ncvar_def(name='rp_curve_rate_upper_ci',
        units='m', dim=list(dim_reference_point_rp_curve),
        missval=NULL,
        longname='rates (upper 95% credible interval) for stage-vs-rate curve at reference point')
       
    all_nc_var = c(all_nc_var, list(reference_point_rp_curve_stages_v,
        reference_point_rp_curve_rates_v, 
        reference_point_rp_curve_rates_lower_ci_v, 
        reference_point_rp_curve_rates_upper_ci_v
        )) 

    # Finished defining variables

    # Make the file
    output_nc_file = nc_create(output_file_name, vars=all_nc_var)

    #
    # Add global attributes
    #
    ncatt_put(output_nc_file, varid=0, attname='file_scenario_name',
        attval=file_scenario_name, prec='text')
    ncatt_put(output_nc_file, varid=0, attname='creation_date',
        attval=as.character(Sys.time()), prec='text')
    ncatt_put(output_nc_file, varid=0, 
        attname='reference_point_coordinates_where_stage_exceedance_rates_were_computed',
        attval=paste(reference_point, collapse=" , "), prec='text')

    #
    # Add nc variables
    #

    ncvar_put(output_nc_file, 'lon', gauge_locations$lon) ; gc()
    ncvar_put(output_nc_file, 'lat', gauge_locations$lat) ; gc()
    ncvar_put(output_nc_file, 'elev', gauge_locations$elev) ; gc()
    ncvar_put(output_nc_file, 'gaugeID', gauge_locations$gaugeID) ; gc()

    ncvar_put(output_nc_file, 'stage', gauge_event_STAGE); gc()
    ncvar_put(output_nc_file, 'uh', gauge_event_UH); gc()
    ncvar_put(output_nc_file, 'vh', gauge_event_VH); gc()

    # Stochastic slip
    ncvar_put(output_nc_file, 'event_Mw', all_eq_events$Mw); gc()
    ncvar_put(output_nc_file, 'event_target_lon', all_eq_events$target_lon); gc()
    ncvar_put(output_nc_file, 'event_target_lat', all_eq_events$target_lat); gc()
    ncvar_put(output_nc_file, 'event_peak_slip_downdip_index', all_eq_events$peak_slip_downdip_ind); gc()
    ncvar_put(output_nc_file, 'event_peak_slip_alongstrike_index', all_eq_events$peak_slip_alongstrike_ind); gc()
    ncvar_put(output_nc_file, 'event_sourcename', all_eq_events$sourcename); gc()
    ncvar_put(output_nc_file, 'event_uniform_event_row', all_eq_events$uniform_event_row); gc()
    ncvar_put(output_nc_file, 'event_index_string', all_eq_events$event_index_string); gc()
    ncvar_put(output_nc_file, 'event_slip_string', all_eq_events$event_slip_string); gc()
    ncvar_put(output_nc_file, 'event_peak_stage_exceedance_rate_at_reference_point', peak_stage_exceedance_rate_at_reference_point); gc()
    ncvar_put(output_nc_file, 'event_peak_stage_at_reference_point', peak_stage_at_reference_point); gc()
    ncvar_put(output_nc_file, 'proportion_rate_from_this_source', proportion_rate_from_this_source); gc()
    ncvar_put(output_nc_file, 'desired_exceedance_rate', desired_exceedance_rate); gc()

    ncvar_put(output_nc_file, 'rp_curve_stage', reference_point_rp_curve_stages); gc()
    ncvar_put(output_nc_file, 'rp_curve_rate', reference_point_rp_curve_rates); gc()
    ncvar_put(output_nc_file, 'rp_curve_rate_lower_ci', reference_point_rp_curve_rates_lower_ci); gc()
    ncvar_put(output_nc_file, 'rp_curve_rate_upper_ci', reference_point_rp_curve_rates_upper_ci); gc()

    nc_close(output_nc_file)

    gc()

    return(invisible(output_file_name))
}

