# Code to make earthquake events
#
# - Tests to help ensure results are sensible:
#   - All the events do have the desired magnitude (back-calculated)
#   - Event areas/lengths/widths are consistent with Mw + scaling relation
#   (make a plot and compare graphically) [Actually at the moment we check this in the function]
#   - Good to have a test that links the rupture size with the deformation

#' Given the summary statistics for a source zone, compute **all** uniform slip events
#' with moment magnitude Mw.
#'
#' The events all have the same number of unit sources along-strike and
#' down-dip. \code{Mw_2_rupture_size} is used to compute the target event
#' dimensions.
#'
#' @param Mw Moment magnitude.
#' @param unit_source_stats Output of
#' discretized_source_approximate_summary_statistics or similar.
#' @param mu Shear modulus in Pascals.
#' @param constant value of constant passed to \code{M0_2_Mw}
#' @param relation The scaling relation used  to determine width,
#' length and area. This is the argument 'relation' passed to
#' \code{Mw_2_rupture_size}
#' @return A list with information which can be used to create all events.
#' @export
get_all_earthquake_events_of_magnitude_Mw<-function(Mw, unit_source_stats, mu=3.0e+10, constant=9.05,
    relation='Strasser'){

    # Quick check that the input table is ok
    if(!all(unit_source_stats$subfault_number == 1:length(unit_source_stats[,1]))){
        msg = paste0('Problem with unit_source_stats table:\n',
            'subfault numbers should range in order from 1 to number_of_subfaults')
        stop(msg)
    }
    

    ## Get scaling-relation based area, width, length
    rupture_stats = Mw_2_rupture_size(Mw, relation=relation, detailed=TRUE, CI_sd = 2.0)
    desired_ALW = rupture_stats$values

    # Record M0 for later testing
    M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

    # FIXME: This area computation needs to be fixed
    subfault_areas = unit_source_stats[,'length']*unit_source_stats[,'width']
    mean_subfault_area = mean(unit_source_stats[,'length']*unit_source_stats[,'width'])
    mean_subfault_width = mean(unit_source_stats[,'width'])
    mean_subfault_length = mean(unit_source_stats[,'length'])

    ## Determine the number of subfaults along the length/width 
    ndip = max(unit_source_stats['downdip_number'])
    nstrike = max(unit_source_stats['alongstrike_number'])
    desired_subfault_count = max(round(desired_ALW['area']/mean_subfault_area), 1)

    if(desired_subfault_count > (ndip*nstrike)){
        print('Reached limit of rupture size')
    }
    

    ## Figure out the desired width/length in terms of number of cells
    nlength = ceiling(desired_ALW['length']/mean_subfault_length)
    nwidth = ceiling(desired_ALW['width']/mean_subfault_width)

    if(nlength > nstrike){
        # Length constrained
        nlength = nstrike
        nwidth = min(max(round(desired_subfault_count/nlength), 1), ndip)
    }else if(nwidth > ndip){
        # Width constrained
        nwidth = ndip
        nlength = min(max(round(desired_subfault_count/nwidth), 1), nstrike)
    }else{
        # Try a few different lengths/widths, and use the one with the best aspect ratio
        l0 = rep(NA, 4)
        w0 = rep(NA, 4)

        desired_aspect_ratio = desired_ALW['length']/desired_ALW['width']

        # Choose length first
        l0[1] = min( max(floor(desired_ALW['length']/mean_subfault_length), 1), nstrike)
        l0[2] = min( max(ceiling(desired_ALW['length']/mean_subfault_length), 1), nstrike)

        w0[1] = min( max(round(desired_subfault_count/l0[1]), 1), ndip)
        w0[2] = min( max(round(desired_subfault_count/l0[2]), 1), ndip)

        # Now choose width first
        w0[3] = min( max(floor(desired_ALW['width']/mean_subfault_width), 1), ndip)
        w0[4] = min( max(ceiling(desired_ALW['width']/mean_subfault_width), 1), ndip)

        l0[3] = min( max(round(desired_subfault_count/w0[3]), 1), nstrike)
        l0[4] = min( max(round(desired_subfault_count/w0[4]), 1), nstrike)

        # Get the 'best'
        # According to Strasser, log10(length/width) = log10(length) - log10(width)
        # is normally distributed, with mean = log10(desired_ALW['length']/desired_ALW['width']),
        # and variance = sum of variances [which can be obtained from the scaling relation
        # Therefore, the 'best' aspect ratio is the one which minimises
        # abs(log10(l0*mean_subfault_length/(w0*mean_subfault_width)) - log10(desired_aspect_ratio))

        aspect_error_ratios = abs(log10((l0*mean_subfault_length)/(w0*mean_subfault_width)) - 
            log10(desired_aspect_ratio))
        #area_error_ratios = abs(log10((l0*w0*mean_subfault_area)/desired_ALW['area']))

        chosen_ind = which.min(aspect_error_ratios)
        #chosen_ind = which.min(area_error_ratios)
        #chosen_ind = which.min(aspect_error_ratios + area_error_ratios)

        nwidth = w0[chosen_ind]
        nlength = l0[chosen_ind]
    }

    # We don't allow > subfault width than length
    #if(nwidth > nlength){
    #    nlength = round(sqrt(desired_subfault_count))
    #    nwidth = nlength
    #}

    #nwidth = round(desired_ALW['width']/mean_subfault_width)
    #nlength = round(desired_subfault_count/nwidth)

    if(nlength > nstrike) nlength = nstrike
    if(nwidth > ndip) nwidth = ndip

    ## Find all possible topleft corners
    topleft_corner_diprange = 1:(ndip - nwidth + 1)
    topleft_corner_strikerange = 1:(nstrike - nlength + 1)

    all_topleft_indices = expand.grid(dipIndex = topleft_corner_diprange, 
        strikeIndex=topleft_corner_strikerange)

    ## Other useful variables
    event_dim = c(nlength, nwidth)
    names(event_dim) = c('length', 'width')

    nevents = length(all_topleft_indices[,1])

    # Store event summary statistics here
    event_statistics = data.frame(area = rep(0.0, nevents), 
        mean_length = rep(0.0, nevents), mean_width = rep(0.0, nevents), 
        slip = rep(0.0, nevents), Mw = rep(Mw, nevents), 
        mean_depth = rep(0.0, nevents), max_depth = rep(0.0, nevents))
    # Store event indices here
    event_indices = list()

    ## Loop over all events and compute their summary statistics
    for(i in 1:nevents){
        topleft_index = all_topleft_indices[i,]

        # Get indices in unit_source_stats which contribute to this event
        event_indices[[i]] = which(
            (unit_source_stats$downdip_number %in% (topleft_index$dipIndex + (1:nwidth - 1)))&
            (unit_source_stats$alongstrike_number %in% (topleft_index$strikeIndex + (1:nlength - 1)))
            )

        # Logical test that we have the right number of subfaults
        if(!(length(event_indices[[i]]) == (nlength*nwidth))){
            print(paste0('Mw: ', Mw))
            print('Event indices: ')
            print(event_indices[[i]])
            print(paste0('length(event_indices[[i]]) = ', length(event_indices[[i]])))
            print(paste0('nlength: ', nlength, ' nwidth: ', nwidth, ' nlength*nwidth: ', nlength*nwidth))
            stop('Incorrect number of subfaults')
        }

        # Shorthand
        ll = unit_source_stats$length[event_indices[[i]]]
        ww = unit_source_stats$width[event_indices[[i]]]
        dd = unit_source_stats$depth[event_indices[[i]]]
        dd_M = unit_source_stats$max_depth[event_indices[[i]]]

        # Store key statistics in table.
        # Note that width*length != area although they will be similar
        event_statistics$area[i] = sum(ll*ww)
        event_statistics$mean_length[i] = sum(ll)/nwidth
        event_statistics$mean_width[i] = sum(ww)/nlength
        event_statistics$slip[i] = slip_from_Mw_area_mu(Mw, event_statistics$area[i], mu, constant=constant)
        event_statistics$mean_depth[i] = mean(dd)
        event_statistics$max_depth[i] = max(dd_M)

        # Logical test that the magnitude is correct
        local_M0 = event_statistics$slip[i] * (event_statistics$area[i] * 1e+06) * mu
        stopifnot(isTRUE(all.equal(local_M0, M0)))

        # Logical test that length/width/area are (very) reasonable
        aa = event_statistics$area[i]
        test1 = ((aa > rupture_stats$minus_CI['area']) &
                 (aa < rupture_stats$plus_CI['area']))

        aa = event_statistics$mean_width[i]
        test2 = ((aa > rupture_stats$minus_CI['width']) &
                 (aa < rupture_stats$plus_CI['width']))

        aa = event_statistics$mean_length[i]
        test3 = ((aa > rupture_stats$minus_CI['length']) &
                 (aa < rupture_stats$plus_CI['length']))

        if(!(test1 & test2 & test3)){
            #browser()
            print('Warning -- event length/width/area not within 2sd of the scaling relation')
        }
    }

    # Useful summary statistic for comparison
    mean_event_ALW = c(mean(event_statistics$area), 
        mean(event_statistics$mean_width),
        mean(event_statistics$mean_length))
    names(mean_event_ALW) = c('area', 'width', 'length')

    output = list(
        topleft_indices = all_topleft_indices,
        event_dim = event_dim,
        desired_subfault_count = desired_subfault_count,
        actual_subfault_count = prod(event_dim),
        desired_ALW = desired_ALW,
        mean_event_ALW = mean_event_ALW,
        event_statistics = event_statistics,
        event_indices = event_indices,
        unit_source_stats = unit_source_stats,
        Mw = Mw)

    return(output)
}


#
# Our current subfaults have a tendency to have lower width where the
# source-zone itself has lower width.
#
# This is forced by the logically rectangular layout.
#
# Suppose Mw is fixed, and each event has uniform slip on all subfaults. We can
# find collections of subfaults with length/width similar to the desired value. 
# A simple way to construct these would be to find the number of subfaults AS and DD
# that most closely matched the mean value.
#
# If we make these have Mw exactly = desired Mw, then events with lower area would have 
# higher slip than those with greater area.
#
# Let i index the subfault, and j index the event. We should have:
#     sum_{j in Mw}(event_probability_j) = probability of an event of size Mw.
# It would also be desirable to have:
#     sum_{j}(event_slip_ij * event_probability_j) = global slip (for every i)
#
# Supposing that all subfaults were involved in the same number of events, we could
# satisfy this relation by having probability_j ~ 1/slip_j
#
# However, we generally can't have all subfaults involved in the same number of
# events (because of boundary effects). Subfaults in the interior are more
# likely to slip than subfaults on the edges. Perhaps we should just live with
# this?
#
# Further, in reality the rates of slip and coupling might differ over the source-zone.
# So we might not really need to make the long-terms slip rate constant everywhere.
#

#' Compute all synthetic uniform slip earthquake events for a sourcezone
#'
#' Apply \code{get_all_earthquake_events_of_magnitude_Mw} to a range of Mw
#' values and combine the results
#'
#' @param discrete_source The discrete source for the source-zone, as from 
#' \code{discretized_source_from_source_contours}. If NULL, you must provide
#' unit_source_statistics (prefered approach)
#' @param unit_source_statistics Pre-computed unit-source-statistics for the
#' discrete-source. It is best to provide this, but if it is NULL, a local
#' computation is made with \code{discretized_source_summary_statistics},
#' assuming approx_dx=approx_dy=5000, and using default arguments for other
#' inputs.  It's normally best to compute this yourself beforehand to take
#' control of those parameters.
#' @param Mmin The minimum Mw value in the table
#' @param Mmax The maximum Mw value in the table
#' @param dMw The increment between the Mw values
#' @param mu Shear modulus in Pascals.
#' @param constant value of constant passed to \code{M0_2_Mw}
#' @param target_location c(longitude, latitude) or NULL. If not NULL, then
#' only return earthquake events which include the unit source 'closest' to the
#' target location (computed assuming a spherical earth with lon/lat coordinates)
#' @param source_zone_name name of the source-zone. This will be
#' put in the table. Even though the entire column is constant, it can be useful 
#' useful when working with tables that combine results from many sources.
#' @param relation Scaling relation type (e.g. 'Strasser'), passed as argument 'relation'
#' to \code{get_all_earthquake_events_of_magnitude_Mw}
#' @return A large table containing information on all earthquake events, and
#' the unit sources they involve
#' @export
#' @examples
#'
#' # Get source contours
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip)
#'
#' puysegur_unit_source_stats = discretized_source_summary_statistics(
#'    puysegur_discretized_source, approx_dx=10000, approx_dy=10000)
#' puysegur_events = get_all_earthquake_events(
#'     unit_source_statistics = puysegur_unit_source_stats,
#'     Mmin = 7.2, Mmax = 9.2, dMw = 0.1, 
#'     source_zone_name='puysegur')
#' print(head(puysegur_events))
#'
get_all_earthquake_events<-function(discrete_source = NULL, unit_source_statistics = NULL,
    Mmin=7.5, Mmax = 9.6, dMw = 0.1, mu=3.0e+10, constant=9.05, target_location = NULL,
    source_zone_name = NA, relation='Strasser'){

    if(is.null(unit_source_statistics)){

        if(is.null(discrete_source)){
            stop('Cannot have both "discrete_source" and "unit_source_statistics" being NULL')
        }

        unit_source_statistics = discretized_source_summary_statistics(discrete_source, 
            approx_dx = 5000, approx_dy = 5000)
    }

    relation = relation

    ## Get all earthquake events
    all_eq_events = lapply(as.list(seq(Mmin, Mmax, dMw)), 
        f<-function(x){ 
            get_all_earthquake_events_of_magnitude_Mw(x, unit_source_statistics, 
                mu=mu, constant=constant, relation=relation)}
        )

    source_zone_name = source_zone_name
    ## Convert to a single table which holds all the events + their subfaults
    add_unit_source_indices_to_event_table<-function(event){
        event_statistics = event$event_statistics

        # Make a character vector, where each entry is a string
        # with all unit source indices for that event, separated by '-'
        event_index_string = unlist(lapply(event$event_indices, 
            f<-function(x) paste0(x, "-", collapse="")))
        local_source_name = rep(source_zone_name, length(event_index_string))

        event_statistics = cbind(event_statistics, 
            data.frame(event_index_string = event_index_string,
                sourcename = source_zone_name))

        return(event_statistics)
    }

    all_eq_tables = lapply(all_eq_events, add_unit_source_indices_to_event_table)

    # big_eq_table + unit_source_statistics hold everything we need about the event geometry
    # (but not probability)
    big_eq_table = all_eq_tables[[1]]
    if(length(all_eq_tables) > 1){
        for(i in 2:length(all_eq_tables)) big_eq_table = rbind(big_eq_table, all_eq_tables[[i]])
    }

    if(!is.null(target_location)){
        # Find the row of unit_source_statistics with location closest to
        # target_location
        N = length(unit_source_statistics[,1])
        unit_source_nearest_target = which.min(distHaversine(
            as.matrix(unit_source_statistics[,c('lon_c', 'lat_c')]),
            cbind(rep(target_location[1], N), rep(target_location[2], N)))) 

        # Find events containing unit_source_nearest_target 
        target_eq_event_indices = which(
            grepl(paste0('-', unit_source_nearest_target, '-'), 
                paste0('-',as.character(big_eq_table$event_index_string))))
    
        big_eq_table = big_eq_table[target_eq_event_indices,]
    }

    return(big_eq_table)
}

#' Find unit sources contained in event
#'
#' The earthquake event table (generated by \code{get_all_earthquake_events})
#' contains a column 'event_index_string' containing the indices of
#' unit sources which contribute to the event, separated by '-'.
#' This convenience function converts those indices to a numeric vector of integers.
#' The latter correspond to the row indices of unit_source_statistics (which
#' is an argument to \code{get_all_earthquake_events}, see the help of that function
#' for details)
#' 
#' @param earthquake_event A single row of the earthquake event table made with
#' \code{get_all_earthquake_events}
#' @param also_return_slip logical - if TRUE the result is a list with 'inds'
#' and 'slip', otherwise we directly return 'inds'
#' @return Integer indices in earthquake_event$event_index_string. These correspond
#' to rows in the unit_source_summary_statistics for the source.
#' @export
get_unit_source_indices_in_event<-function(earthquake_event, also_return_slip=FALSE){

    if(nrow(earthquake_event) != 1){
        stop('only a single row of the earthquake event table can be passed as input')
    }

    unit_source_indices_text = gsub('-', ' ', as.character(earthquake_event$event_index_string))
    event_unit_source_inds = scan(text=unit_source_indices_text, quiet=TRUE)

    if(!also_return_slip){
        return(event_unit_source_inds)
    }else{
        # Get the slip as well
        if('event_slip_string' %in% names(earthquake_event)){
            slip = as.numeric(strsplit(earthquake_event$event_slip_string, '_')[[1]])
        }else{
            slip = rep(earthquake_event$slip, length(event_unit_source_inds))
        }

        output = list(inds=event_unit_source_inds, slip=slip)
        return(output)
    }
}

#' Earthquake centroid calculation, weighted by slip
#'
#' Compute either: 1) the weighted mean of the unit-source locations lon_c and lat_c 
#' (using the angular mean) or; 2) a weighted mean of the downdip and alongstrike 
#' unit_source indices. In this case we return the rounded values, and also values
#' prior to rounding (to give the user greater control in corner cases).
#'
#' @param event_row A single row of the event table. The code can work with either 
#' uniform_slip or variable/stochastic slip formats.
#' @param unit_source_statistics data.frame giving the unit-source summary statistics.
#' This could be the output from the function \code{discretized_source_summary_statistics}.
#' @param as_subfault_number If FALSE then return a lon/lat coordinate. If TRUE, then
#' compute the weighted downdip_number and alongstrike_number. These might not be integers,
#' so we return the unrounded values as well as the rounded values, and the subfault_number
#' corresponding to the rounded values (i.e. the corresponding row of the unit_source_statistics)
#' @return either a coordinate (subfault_number == FALSE) or a list with the downdip_number, 
#' alongstrike_number, subfault_number, mean_downdip_number (i.e. before rounding), 
#' mean_alongstrike_number (i.e. before rounding).
#' @export
#' @examples
#' puysegur = readOGR(system.file('extdata/puysegur.shp', package='rptha'), layer='puysegur')
#' # Get downdip lines
#' puysegur_downdip = readOGR(system.file('extdata/puysegur_downdip.shp', package='rptha'), 
#'    layer='puysegur_downdip')
#' # Make discretized_source with 50km x 50km unit-sources (approximately)
#' puysegur_discretized_source = discretized_source_from_source_contours(
#'     source_shapefile=puysegur,
#'    desired_subfault_length=50,
#'    desired_subfault_width=50,
#'    downdip_lines=puysegur_downdip)
#' 
#' puysegur_unit_source_stats = discretized_source_summary_statistics(
#'    puysegur_discretized_source, approx_dx=10000, approx_dy=10000)
#' puysegur_events = get_all_earthquake_events(
#'     unit_source_statistics = puysegur_unit_source_stats,
#'     Mmin = 7.2, Mmax = 9.2, dMw = 0.1, 
#'     source_zone_name='puysegur')
#' 
#' #
#' # Example earthquake with a single unit-source
#' #
#' 
#' for(event_ind in seq(1, nrow(puysegur_events))){
#'     pinds = as.numeric(strsplit(puysegur_events$event_index_string[event_ind], '-')[[1]])
#' 
#'     # Get the coordinates
#'     p = get_event_slip_weighted_centroid(puysegur_events[event_ind,], puysegur_unit_source_stats)
#'     # Check the coordinates, being careful about rounding issues
#'     stopifnot(isTRUE(all.equal(p[1], mean_angle(puysegur_unit_source_stats$lon_c[pinds]))) & 
#'               isTRUE(all.equal(p[2],  mean_angle(puysegur_unit_source_stats$lat_c[pinds]))))
#' 
#'     # Get the subfault
#'     p = get_event_slip_weighted_centroid(puysegur_events[event_ind,], puysegur_unit_source_stats, 
#'                                          as_subfault_number=TRUE)
#'     # To avoid issues with rounding and tiny floating point issues, here we check the unrounded
#'     # indices 
#'     stopifnot(isTRUE(all.equal(p$mean_downdip_number, 
#'                                mean(puysegur_unit_source_stats$downdip_number[pinds]))) &
#'               isTRUE(all.equal(p$mean_alongstrike_number, 
#'                                mean(puysegur_unit_source_stats$alongstrike_number[pinds]))))
#' }
#' 
get_event_slip_weighted_centroid<-function(event_row, unit_source_statistics, as_subfault_number=FALSE){

    inds_slip = get_unit_source_indices_in_event(event_row, also_return_slip=TRUE)
    inds = inds_slip$inds
    slip = inds_slip$slip

    if(!as_subfault_number){
        # Get the weighted mean angle of the unit sources
        lon = unit_source_statistics$lon_c[inds]
        lat = unit_source_statistics$lat_c[inds]
        mean_lon = mean_angle(lon, weights=slip)
        mean_lat = mean_angle(lat, weights=slip)
        return(c(mean_lon, mean_lat))

    }else{
        # Get the nearest along-strike and down-dip indices
        mean_downdip_number = weighted.mean(unit_source_statistics$downdip_number[inds], w=slip)
        dd = round(mean_downdip_number)
        mean_alongstrike_number = weighted.mean(unit_source_statistics$alongstrike_number[inds], w=slip)
        as = round(mean_alongstrike_number)
        nearest_subfault = which(unit_source_statistics$alongstrike_number == as &
                                 unit_source_statistics$downdip_number == dd)
        output = list(subfault_number=nearest_subfault, alongstrike_number=as, 
                      downdip_number=dd, mean_downdip_number=mean_downdip_number,
                      mean_alongstrike_number = mean_alongstrike_number)
        return(output)
    }

}


#' Graphical checks of earthquake events
#'
#' Plot some properties of the earthquake events, with comparisons to scaling
#' relations (solid lines) and +- 2sd (dashed lines)
#'
#' @param all_eq_events data.frame with all earthquake events. This could be
#' generated using the function \code{get_all_earthquake_events}
#' @param relation the scaling relation type used for theoretical comparison (see \code{Mw_2_rupture_size})
#' @return nothing, but make a plot
#' @export
plot_earthquake_event_properties<-function(all_eq_events, relation='Strasser'){
    
    par(mfrow=c(3,2))

    # Initial plot showing numbers of events
    plot(table(all_eq_events$Mw), main='Number of events by magnitude')

    # Plot of event Mw vs mean depth
    plot(all_eq_events$Mw, all_eq_events$mean_depth, 
        main='Mw vs mean rupture depth', xlab='Mw', ylab='depth (km)')

    # Compute scaling relations, which our events should roughly satisfy
    mw_seq = sort(unique(all_eq_events$Mw))
    mw_area_theory = mw_seq*NA
    mw_area_theory_p2sd = mw_seq*NA
    mw_area_theory_m2sd = mw_seq*NA
    mw_length_theory = mw_seq*NA
    mw_length_theory_p2sd = mw_seq*NA
    mw_length_theory_m2sd = mw_seq*NA
    mw_width_theory = mw_seq*NA
    mw_width_theory_p2sd = mw_seq*NA
    mw_width_theory_m2sd = mw_seq*NA

    for(i in 1:length(mw_seq)){

        scaling_results = Mw_2_rupture_size(mw_seq[i], relation=relation, detailed=TRUE, CI_sd=2)

        mw_area_theory[i] = scaling_results$values['area']
        mw_area_theory_p2sd[i] = scaling_results$plus_CI['area']
        mw_area_theory_m2sd[i] = scaling_results$minus_CI['area']

        mw_width_theory[i] = scaling_results$values['width']
        mw_width_theory_p2sd[i] = scaling_results$plus_CI['width']
        mw_width_theory_m2sd[i] = scaling_results$minus_CI['width']
        
        mw_length_theory[i] = scaling_results$values['length']
        mw_length_theory_p2sd[i] = scaling_results$plus_CI['length']
        mw_length_theory_m2sd[i] = scaling_results$minus_CI['length']

    }

    # Main plot
    plot(all_eq_events$Mw, all_eq_events$area, log='y', main='Mw-area', 
        xlab='Mw', ylab='area')
    points(mw_seq, mw_area_theory, t='o', col='red')
    points(mw_seq, mw_area_theory_p2sd, t='o', lty='dashed', col='red')
    points(mw_seq, mw_area_theory_m2sd, t='o', lty='dashed', col='red')

    plot(all_eq_events$Mw, all_eq_events$mean_length, log='y', 
        main='Mw-length', xlab='Mw', ylab='length')
    points(mw_seq, mw_length_theory, t='o', col='red')
    points(mw_seq, mw_length_theory_p2sd, t='o', lty='dashed', col='red')
    points(mw_seq, mw_length_theory_m2sd, t='o', lty='dashed', col='red')

    plot(all_eq_events$Mw, all_eq_events$mean_width, log='y', main='Mw-width', 
        xlab='Mw', ylab='width')
    points(mw_seq, mw_width_theory, t='o', col='red')
    points(mw_seq, mw_width_theory_p2sd, t='o', lty='dashed', col='red')
    points(mw_seq, mw_width_theory_m2sd, t='o', lty='dashed', col='red')

    plot(all_eq_events$Mw, all_eq_events$slip, log='y', main='Mw-slip', 
        xlab='Mw', ylab='slip')
    points(mw_seq, slip_from_Mw(mw_seq), t='o', col='red')

}
