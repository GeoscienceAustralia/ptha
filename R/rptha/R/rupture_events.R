# Code to make earthquake events
#
# FIXME: Need to carefully review discretized_source_approximate_summary_statistics
# e.g. 
# - needs a 'depth_in_km' argument [also there is a magic number 1000 which
# needs to account for this]
# - the area computation should be done more carefully (actually it is a real
# area, and is equal to length*width, but check all this again)
# - Potentially the computations should employ integration over the sub-source
# points
#
#
# - Tests to help ensure results are sensible:
#   - All the events do have the desired magnitude (back-calculated)
#   - Event areas/lengths/widths are consistent with Mw + scaling relation
#   (make a plot and compare graphically)
#   - Good to have a test that links the rupture size with the deformation
#
# TODO #
# - Integrate changes in rupture_scaling to main code
# -  
#
# Mw = 8.0
#
#' Given the summary statistics for a source zone, compute **all** events
#' with moment magnitude Mw.
#'
#' The events all have the same number of unit sources along-strike and down-dip.
#'
#'
#' @param Mw Moment magnitude.
#' @param unit_source_stats Output of discretized_source_approximate_summary_statistics or similar.
#' @param mu Shear modulus in Pascals.
#' @return A list with information which can be used to create all events.
#' 
get_all_events_of_magnitude_Mw<-function(Mw, unit_source_stats, mu=3.0e+10){

    ## Get scaling-relation based area, width, length
    rupture_stats = Mw_2_rupture_size(Mw, detailed=TRUE, CI_sd = 2.0)
    desired_ALW = rupture_stats$values

    # Record M0 for later testing
    M0 = M0_2_Mw(Mw, inverse=TRUE)

    # FIXME: This area computation needs to be fixed
    subfault_areas = unit_source_stats[,'length']*unit_source_stats[,'width']
    mean_subfault_area = mean(unit_source_stats[,'length']*unit_source_stats[,'width'])
    mean_subfault_width = mean(unit_source_stats[,'width'])
    mean_subfault_length = mean(unit_source_stats[,'length'])

    ## Determine the number of subfaults along the length/width 
    ndip = max(unit_source_stats['downdip_number'])
    nstrike = max(unit_source_stats['alongstrike_number'])
    desired_subfault_count = round(desired_ALW['area']/mean_subfault_area)

    # Round length towards the longer side of the target
    nlength = round(desired_ALW['length']/mean_subfault_length)
    nwidth = round(desired_subfault_count/nlength)
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
        slip = rep(0.0, nevents), Mw = rep(Mw, nevents))
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
        stopifnot(length(event_indices[[i]]) == (nlength*nwidth))

        # Shorthand
        ll = unit_source_stats$length[event_indices[[i]]]
        ww = unit_source_stats$width[event_indices[[i]]]

        # Store key statistics in table.
        # Note that width*length != area although they will be similar
        event_statistics$area[i] = sum(ll*ww)
        event_statistics$mean_length[i] = sum(ll)/nwidth
        event_statistics$mean_width[i] = sum(ww)/nlength
        event_statistics$slip[i] = slip_from_Mw_area_mu(Mw, event_statistics$area[i], mu)

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
            browser()
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

