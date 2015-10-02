context('test_rupture_creation_and_probabilities')

test_that("test_rupture_creation_and_probabilities", {
    # Alaska source zone computations

    alaska_contour_shapefile = './testshp/alaska.shp'

    # Make discrete source
    discrete_source = discretized_source_from_source_contours(alaska_contour_shapefile,
        desired_subfault_length=100, desired_subfault_width=50, make_plot=FALSE)
      
    # Get summary statistics for unit sources 
    unit_source_summary_statistics = discretized_source_summary_statistics(discrete_source,
        approx_dx=5000, approx_dy = 5000) 

    # Get a table with all earthquake events between Mw = 7.5 and Mw = 9.4
    dMw = 0.1
    earthquake_event_table = get_all_earthquake_events(discrete_source, unit_source_summary_statistics,
        Mmin = 7.5, Mmax = 9.4, dMw = dMw)

    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        earthquake_event_table,
        conditional_probability_model = 'inverse_slip')

    # Recurrence parameters
    slip_rate = c(44.00, 49.50, 55.00)/1000 # m/year
    slip_rate_prob = c(0.6, 0.2, 0.2)  
    # The b parameters from GEM are indeed based on Gutenberg-Richter in base 10
    # (not e) -- according to the text in the GEM report
    #b = c(0.54, 0.40, 0.80) # Parameters from the 2013 report
    # Updated GEM parameters [Berryman et al 2015] -- big difference!
    ## FIXME: Go through all our parameters and compare with new GEM ones. There may
    ## have been mistakes in the early publication -- the 'b' values seem to have changed greatly
    b = c(0.95, 0.7, 1.2)
    b_prob = c(0.6, 0.2, 0.2)
    Mw_min = 7.5
    Mw_min_prob = 1
    Mw_max = c(9.40, 9.40, 9.00)
    Mw_max_prob = c(0.6, 0.3, 0.1)

    # Get total area of source-zone in km^2
    sourcezone_total_area = sum(unit_source_summary_statistics$length*unit_source_summary_statistics$width)

    rate_function = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate,
        slip_rate_prob = slip_rate_prob,
        b = b,
        b_prob = b_prob,
        Mw_min = Mw_min,
        Mw_min_prob = Mw_min_prob,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = earthquake_event_table,
        event_conditional_probabilities = event_conditional_probabilities)

    # Test the rate of 'Great earthquakes' is reasonable
    # These bounds may change but they work at present (Oct 2) and are
    # consistent with other estimates for the "Alaska" sourcezone
    freq_gt9 = 1/rate_function(9.0)
    expect_that( (freq_gt9 > 600) & (freq_gt9 < 700), is_true())

    # Back-calculate slip on each fault
    event_rate = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw-dMw/2) - 
        rate_function(earthquake_event_table$Mw + dMw/2))
    event_longterm_slip = earthquake_event_table$slip * event_rate
    unit_source_longterm_slip = rep(0, length(unit_source_summary_statistics[,1]))

    event_longterm_slip_area = event_longterm_slip * earthquake_event_table$area
    unit_source_longterm_slip_area = rep(0, length(unit_source_summary_statistics[,1]))

    # Loop over all events
    for(ee in 1:length(event_longterm_slip)){
        # Find the unit sources on this event
        unit_sources = strsplit(as.character(earthquake_event_table$event_index_string[ee]), '-')[[1]]
        unit_sources = as.numeric(unit_sources)

        # Add the slip to all the involved unit sources
        unit_source_longterm_slip[unit_sources] = 
            unit_source_longterm_slip[unit_sources] + event_longterm_slip[ee]

        # Add the seismic moment to all the involved sources
        unit_source_longterm_slip_area[unit_sources] = 
            unit_source_longterm_slip_area[unit_sources] + event_longterm_slip_area[ee]
    }
    m1 = matrix(unit_source_longterm_slip, ncol=6,nrow=7,byrow=TRUE)
    m2 = matrix(unit_source_longterm_slip_area, ncol=6,nrow=7,byrow=TRUE)

    # Check that the mean slip is close to the desired value
    # FIXME: Why is it not identical?
    theoretical_slip = sum(slip_rate*slip_rate_prob)
    expect_that(abs(mean(m1) - theoretical_slip)/theoretical_slip < 0.01, is_true())

})
