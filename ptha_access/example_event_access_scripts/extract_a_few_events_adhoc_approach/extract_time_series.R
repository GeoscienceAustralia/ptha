#
# This script takes a given point and a set of max-stage return periods, and uses
# an ad-hoc method to pick a few events from one or more source-zones.
#

#
# Inputs
#

ptha = new.env()
# EDIT THE FOLLOWING PATH AS REQUIRED
source('../../../../ptha/ptha_access/get_PTHA_results.R', local=ptha, chdir=TRUE)

# Hazard point here will be used to define return periods -- prefer points far offshore
return_period_point = c(154.6667, -26.6667)
# Only get wave time-series at points deeper than this (m). Prefer deep points (e.g. 1000+m)
lower_depth_limit = 200 

# We extract wave time series at hazard points in this polygon. It can have an
# arbitrary number of points (but the download will take longer)
# If you only want time-series at the return period point, then set this to NULL
# Note that if the return period point is not inside the polygon, will
# anyway include it 
point_extraction_polygon = rbind(
    c(154.5, -28.8),   # Lower left
    c(154.8, -28.8), # Lower right
    c(154.8, -24.5),  # Upper right
    c(154.5, -24.5)  # Upper left
    )

# Rates at which outputs are desired
all_desired_rates = c(1/500, 1/2500) #c(1/750, 1/3000, 1/10000)
# Alternative to specifying rates -- specify stages directly. In that case, the
# rates should be NA. 
all_desired_stages = c(NA, NA) # c(NA, NA, NA)

# Choose at most this many source-zones at each return period
max_number_of_source_zones = 2 # 4

# Specify the number of events that should be selected for each
# source-zone/rate combination
number_of_events_per_source_and_rate = 2 # 5

# Search for events that have max-stage within a given fraction of the desired_stage,
# or within the given absolute tolerance -- whichever is more generous.
max_stage_relative_tolerance = 0.05 # 0.1
max_stage_absolute_tolerance = 0.01 # 0.05

#
# End Inputs
#

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if(length(all_desired_stages) != length(all_desired_rates)){
    stop('length(all_desired_stages) does not equal length(all_desired_rates) -- but it should')
}

if(any(is.na(all_desired_stages) + is.na(all_desired_rates) != 1)){
    stop("Please specify EITHER all_desired_stages, OR all_desired_rates, with the other being a vector of NA's with the same length")
}

# Get stage-vs-exceedance at the return_period_point
return_period_info = ptha$get_stage_exceedance_rate_curve_at_hazard_point(
    target_point=return_period_point, make_plot=TRUE)

# Get indices of gauges at which we will extract waves
if(!all(is.null(point_extraction_polygon))){

    wave_target_indices = ptha$get_netcdf_gauge_indices_in_polygon(
        return_period_info$stage_exceedance_rate_curves_file,
        point_extraction_polygon)
    # FIXME: For a web app, consider limiting the number of
    # points here [e.g. to 20].
    if(!(return_period_info$target_index %in% wave_target_indices)){
        wave_target_indices = c(wave_target_indices, return_period_info$target_index)
    }

}else{
    wave_target_indices = return_period_info$target_index
}
# Get the coordinates and depth of the wave_target_indices
wave_target_points = ptha$get_netcdf_gauge_locations(
    return_period_info$stage_exceedance_rate_curves_file,
    indices_of_subset = wave_target_indices)
# Filter by depth -- but keep the return period point
k = which((wave_target_points$elev < ((-1)*lower_depth_limit)) |
          (wave_target_indices == return_period_info$target_index))
wave_target_indices = wave_target_indices[k]
wave_target_points = wave_target_points[k,]

#
#
# Using the above plot, the user can provide either a stage, or a return period
#
#

# Get return periods for all sources -- used repeatedly in the loop below we
# will repeatedly use this
all_source_return_periods = ptha$get_stage_exceedance_rate_curves_all_sources(
    target_index = return_period_info$target_index,
    only_mean_rate_curve=TRUE)

output_rdata_file = rep("", length(all_desired_rates)) 
for(ir in 1:length(all_desired_rates)){

    desired_rate = all_desired_rates[ir]
    desired_stage = all_desired_stages[ir]

    # If stage has not been provided, compute it
    if(is.na(desired_stage)){
        desired_stage = approx(return_period_info$stochastic_slip_rate, 
            return_period_info$stage, xout=desired_rate, ties=min)$y
    }
    # If desired_rate has not been provided, compute it
    if(is.na(desired_rate)){
        desired_rate = approx(
            return_period_info$stage,
            return_period_info$stochastic_slip_rate, 
            xout=desired_stage)$y
    }
    # Check everything is consistent
    if(is.na(desired_stage) | is.na(desired_rate)){
        stop('No information available for the desired stage/return_period')
    }

    # Find domainant source-zones, by finding the exceedance rate of
    # 'desired_stage' at the return_period_point on a source-by-source basis
    desired_stage_exceedance_rate_by_source = unlist(lapply(all_source_return_periods,
        f<-function(x) approx(x$stage, x$stochastic_slip_rate, xout=desired_stage)$y))
    # Order by decreasing rate
    desired_stage_exceedance_rate_by_source = sort(desired_stage_exceedance_rate_by_source, 
        decreasing=TRUE)

    # Check that the summed exceedance rate over all sources is indeed equal to the desired rate.
    # Otherwise something has gone wrong
    stopifnot(isTRUE(all.equal(sum(desired_stage_exceedance_rate_by_source), desired_rate, tol=1.0e-06)))

    # Choose source_zones with non-zero rate -- up to 4 here
    n = min(max(1, max(which(desired_stage_exceedance_rate_by_source > 0))), max_number_of_source_zones)
    source_zones_of_interest = names(desired_stage_exceedance_rate_by_source)[1:n]
    rates_by_source_zone_of_interest = desired_stage_exceedance_rate_by_source[1:n]

    #
    # Extract events at source-zones of interest that produce stage close to
    # the desired stage
    #
    sourcezone_event_stage_Mw = ptha$get_peak_stage_at_point_for_each_event(
        target_index=return_period_info$target_index,
        all_source_names = source_zones_of_interest)
    relative_tolerance = max_stage_relative_tolerance 
    absolute_tolerance = max_stage_absolute_tolerance 
    # Indices of events in the earthquake_events_table
    sourcezone_event_inds_near_desired_stage = lapply(
        sourcezone_event_stage_Mw,
        f<-function(x){
            which((abs(x$max_stage - desired_stage) < 
                max(relative_tolerance*desired_stage, absolute_tolerance))&
                (x$scenario_rate_is_positive))
        }
    )
    names(sourcezone_event_inds_near_desired_stage) = source_zones_of_interest

    # Check we did not have zero events. If we do, remove that source-zone from the list
    n_near_events = unlist(lapply(sourcezone_event_inds_near_desired_stage, length))
    if(any(n_near_events == 0)){
        keepers = which(n_near_events > 0)
        if(length(keepers) == 0) next  # Break out of overall loop

        # Remove the event with no points 
        source_zones_of_interest = source_zones_of_interest[keepers]
        rates_by_source_zone_of_interest = rates_by_source_zone_of_interest[keepers]
        sourcezone_event_stage_Mw = sourcezone_event_stage_Mw[keepers]
        sourcezone_event_inds_near_desired_stage = sourcezone_event_inds_near_desired_stage[keepers]
        n = length(keepers)
        n_near_events = n_near_events[keepers]
        print(c('Removing source zones without close events. Keeping ', source_zones_of_interest) )

    }


    # Extract the actual earthquake events
    sourcezone_events_near_desired_stage = vector(mode='list', length=n)
    names(sourcezone_events_near_desired_stage) = source_zones_of_interest
    for(i in 1:n){
        sourcezone_events_near_desired_stage[[i]] = ptha$get_source_zone_events_data(
            source_zones_of_interest[i], slip_type='stochastic', 
            desired_event_rows=sourcezone_event_inds_near_desired_stage[[i]])
    }

    #
    # Choose some events from the above selection
    #
    event_summaries = lapply(sourcezone_events_near_desired_stage, ptha$summarise_events)
    # Ad-hoc event selection method based on mahalanobis distance. This should give events
    # with properties that are more 'central'. However, it has little justification, need
    # to study event reduction
    chosen_events = lapply(event_summaries, f<-function(x){
        desired_nevents = number_of_events_per_source_and_rate
        n = min(desired_nevents, length(x$mh_distance)) # In case we do not have enough events
        output = order(x$mh_distance)[1:n]
        return(output) 
        }
    )


    # Get the tsunami
    sourcezone_waves = vector(mode='list', length=n)
    names(sourcezone_waves) = source_zones_of_interest
    for(i in 1:n){
        print(i)
        sourcezone_waves[[i]] = ptha$get_flow_time_series_at_hazard_point(
            sourcezone_events_near_desired_stage[[i]], 
            #event_ID=1:nrow(sourcezone_events_near_desired_stage[[i]]$events),
            event_ID = chosen_events[[i]],
            target_indices=wave_target_indices)
    }

    # TEST
    # 1) Reconstructed wave heights agree with the target range, and with the
    #   heights extracted from the summary info. This requires that the 'return_period_point'
    #   is also a wave_target_index
    for(i in 1:length(sourcezone_waves)){
        sw = sourcezone_waves[[i]]
        site = which.min(abs(c(sw$locations[,1]) - return_period_info$lon[1]) + abs(c(sw$locations[,2]) - return_period_info$lat[1]))
        flow = sw$flow[[site]]
        max_stage_summary = sourcezone_event_stage_Mw[[i]]$max_stage[sourcezone_event_inds_near_desired_stage[[i]]]
        for(j in 1:length(chosen_events[[i]])){
            stopifnot(isTRUE(all.equal(max(flow[j,,1]), max_stage_summary[chosen_events[[i]][j]], tolerance=1.0e-05)))
        }
        print(paste0('PASS peak wave height test', i))
    }

    # TEST
    # 2) No events have zero rate.
    for(i in 1:length(sourcezone_events_near_desired_stage)){
        stopifnot(all(sourcezone_events_near_desired_stage[[i]]$events$rate_annual > 0))
        print(paste0('PASS rate check', i))
    }


    output_rdata_file_local = paste0('time_series_extraction_image_', return_period_point[1], 
        '_', return_period_point[2],'_', signif(return_period_info$elev,4),
        '_', signif(desired_rate, 4), '_', signif(desired_stage, 4), 
        '_.Rdata')
    save.image(file=output_rdata_file_local)

    output_rdata_file[ir] = output_rdata_file_local

    # Save memory
    rm(sourcezone_waves, event_summaries, chosen_events, sourcezone_event_stage_Mw,
        sourcezone_event_inds_near_desired_stage, sourcezone_events_near_desired_stage)
    gc()
}


#
# netcdf output
#
source('save_gauges_timeseries_to_netcdf.R')

# remove unwritten files
k = which(output_rdata_file != '')
full_output_rdata_file = output_rdata_file[k]

# Load all the Rdata files into a list of R images
full_images = vector(mode='list', length=length(full_output_rdata_file))
for(i in 1:length(full_output_rdata_file)){
    full_images[[i]] = new.env()
    load(full_output_rdata_file[i], envir=full_images[[i]])
}

# Loop over return periods
for(i in 1:length(full_images)){

    im = full_images[[i]]

    # Loop over source-zones
    for(j in 1:length(im$sourcezone_waves)){

        sourcezone = names(im$sourcezone_waves)[j]

        # Loop over events, and write a single netcdf for each
        for(k in 1:nrow(im$sourcezone_waves[[j]]$events)){

            # Get the event metadata
            event = im$sourcezone_waves[[j]]$events[k,]
            # Gauge locations
            gauge_locations = im$sourcezone_waves[[j]]$locations

            # Pack the flow data into suitably dimensioned arrays
            # Since we only store 1 event with this output format, 
            # the final index of the packed data is 1. Although we can
            # store more than one event in the netcdf file, I wonder if
            # it will be confusing for users?
            gauge_obs_times = im$sourcezone_waves[[j]]$time
            gauge_event_STAGE = array(NA, 
                dim=c(length(gauge_obs_times), nrow(event), nrow(gauge_locations)))
            gauge_event_UH = array(NA, 
                dim=c(length(gauge_obs_times), nrow(event), nrow(gauge_locations)))
            gauge_event_VH = array(NA, 
                dim=c(length(gauge_obs_times), nrow(event), nrow(gauge_locations)))

            for(ll in 1:nrow(gauge_locations)){
                gauge_event_STAGE[,1,ll] = im$sourcezone_waves[[j]]$flow[[ll]][k,,1]
                gauge_event_UH[,1,ll] = im$sourcezone_waves[[j]]$flow[[ll]][k,,2]
                gauge_event_VH[,1,ll] = im$sourcezone_waves[[j]]$flow[[ll]][k,,3]
            }

            # Find the peak stage at the reference point for this event, so that we can 
            # report on the actual return period as well as the desired one (they should
            # be similar but not exactly equal, to allow us to extract a few events)
            event_index = im$sourcezone_event_inds_near_desired_stage[[j]][ im$chosen_events[[j]][k] ]
            peak_stage_at_reference_point = im$sourcezone_event_stage_Mw[[j]]$max_stage[event_index]
            peak_stage_exceedance_rate_at_reference_point = approx(
                im$return_period_info$stage, im$return_period_info$stochastic_slip_rate,
                xout=peak_stage_at_reference_point)$y
            proportion_rate_from_this_source = approx(
                im$all_source_return_periods[[sourcezone]]$stage, 
                im$all_source_return_periods[[sourcezone]]$stochastic_slip_rate,
                xout=peak_stage_at_reference_point)$y / peak_stage_exceedance_rate_at_reference_point

            # File name and scenario name
            dir.create('./event_time_series', showWarnings=FALSE)
            output_file_name = paste0('./event_time_series/', sourcezone, '_tsunami_event_Mw_', 
                round(event$Mw[1], 2), '_approx_stage_exceedance_rate_', signif(im$desired_rate,4), 
                '_relative_importance_of_this_sourcezone_', signif(proportion_rate_from_this_source, 3),
                '_index_', event_index, 
                '.nc')
            file_scenario_name = paste0(sourcezone, '_tsunami_event_Mw_', 
                round(event$Mw[1], 2), '_approx_peak_stage_exceedance_rate_', signif(im$desired_rate,5),
                '_at_reference_point_', round(im$return_period_point[1], 4), '_', 
                round(im$return_period_point[2], 4), '_depth_', -round(im$return_period_info$elev,2),
                '_relative_importance_of_this_sourcezone_', signif(proportion_rate_from_this_source, 3),
                '_index_', event_index)

            # write to netcdf
            write_tsunami_gauges_to_netcdf(
                file_scenario_name = file_scenario_name,
                all_eq_events = event,
                gauge_locations = gauge_locations,
                gauge_obs_times = gauge_obs_times,
                gauge_event_STAGE = gauge_event_STAGE,
                gauge_event_UH = gauge_event_UH,
                gauge_event_VH = gauge_event_VH,
                peak_stage_at_reference_point = peak_stage_at_reference_point,
                peak_stage_exceedance_rate_at_reference_point = peak_stage_exceedance_rate_at_reference_point,
                proportion_rate_from_this_source = proportion_rate_from_this_source,
                reference_point = im$return_period_point,
                reference_point_rp_curve_stages = im$return_period_info$stage,
                reference_point_rp_curve_rates = im$return_period_info$stochastic_slip_rate,
                reference_point_rp_curve_rates_lower_ci = im$return_period_info$stochastic_slip_rate_lower_ci,
                reference_point_rp_curve_rates_upper_ci = im$return_period_info$stochastic_slip_rate_upper_ci,
                desired_exceedance_rate = im$desired_rate,
                output_file_name = output_file_name)
        }
    }
}


