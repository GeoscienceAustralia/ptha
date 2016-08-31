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
#' @param unit_source_stats Output of discretized_source_approximate_summary_statistics or similar.
#' @param mu Shear modulus in Pascals.
#' @param constant value of constant passed to \code{M0_2_Mw}
#' @return A list with information which can be used to create all events.
#' @export
get_all_earthquake_events_of_magnitude_Mw<-function(Mw, unit_source_stats, mu=3.0e+10, constant=9.05){

    ## Get scaling-relation based area, width, length
    rupture_stats = Mw_2_rupture_size(Mw, detailed=TRUE, CI_sd = 2.0)
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
#' \code{discrete_source_from_source_contours}. If NULL, you must provide
#' unit_source_statistics
#' @param unit_source_statistics Pre-computed unit-source-statistics for the
#' discrete-source. If NULL, a local computation is made with
#' \code{discretized_source_summary_statistics}, assuming
#' approx_dx=approx_dy=5000, and using default arguments for other inputs. 
#' It's normally best to compute this yourself beforehand to take control of
#' those parameters.
#' @param Mmin The minimum Mw value in the table
#' @param Mmax The maximum Mw value in the table
#' @param dMw The increment between the Mw values
#' @param mu Shear modulus in Pascals.
#' @param constant value of constant passed to \code{M0_2_Mw}
#' @param target_location c(longitude, latitude) or NULL. If not NULL, then
#' only return earthquake events which include the unit source 'closest' to the
#' target location (computed assuming a spherical earth with lon/lat coordinates)
#' @return A large table containing information on all earthquake events, and
#' the unit sources they involve
#' @export
get_all_earthquake_events<-function(discrete_source = NULL, unit_source_statistics = NULL,
    Mmin=7.5, Mmax = 9.6, dMw = 0.1, mu=3.0e+10, constant=9.05, target_location = NULL){

    if(is.null(unit_source_statistics)){

        if(is.null(discrete_source)){
            stop('Cannot have both "discrete_source" and "unit_source_statistics" being NULL')
        }

        unit_source_statistics = discretized_source_summary_statistics(discrete_source, 
            approx_dx = 5000, approx_dy = 5000)
    }

    ## Get all earthquake events
    all_eq_events = lapply(as.list(seq(Mmin, Mmax, dMw)), 
        f<-function(x){ 
            get_all_earthquake_events_of_magnitude_Mw(x, unit_source_statistics, 
                mu=mu, constant=constant)}
        )

    ## Convert to a single table which holds all the events + their subfaults
    add_unit_source_indices_to_event_table<-function(event){
        event_statistics = event$event_statistics

        # Make a character vector, where each entry is a string
        # with all unit source indices for that event, separated by '-'
        event_index_string = unlist(lapply(event$event_indices, 
            f<-function(x) paste0(x, sep="-", collapse="")))

        event_statistics = cbind(event_statistics, 
            data.frame(event_index_string = event_index_string))

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

