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

    # Get total area of source-zone in km^2
    sourcezone_total_area = sum(unit_source_summary_statistics$length *
        unit_source_summary_statistics$width)


    ###########################################################################
    #
    # A 'REALISTIC' TEST
    #
    ###########################################################################

    # Recurrence parameters
    slip_rate = c(44.00, 49.50, 55.00)/1000 # m/year
    slip_rate_prob = c(0.6, 0.2, 0.2)  
    b = c(0.95, 0.7, 1.2)
    b_prob = c(0.6, 0.2, 0.2)
    Mw_min = 7.5
    Mw_min_prob = 1
    Mw_max = c(9.40, 9.40, 9.00)
    Mw_max_prob = c(0.6, 0.3, 0.1)

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
    # These bounds may change but they work at present (2/10/2015) and are
    # consistent with other estimates for the "Alaska" sourcezone
    freq_gt9 = 1/rate_function(9.0)
    expect_that( (freq_gt9 > 600) & (freq_gt9 < 700), is_true() )
   
    # Check that if we account for the moment below mwmin, then the rate
    # of earthquakes is less. 
    rate_functionB = rate_of_earthquakes_greater_than_Mw_function(
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
        event_conditional_probabilities = event_conditional_probabilities,
        account_for_moment_below_mwmin=TRUE)

    freq_gt9B = 1/rate_functionB(9.0)
    expect_that( freq_gt9 < freq_gt9B, is_true())

    #
    # Check that we get pretty much the same as rate_functionB by setting Mw_min to 0,
    # and ignoring the moment below mwmin.
    #
    # To do this test, we need to hack events with Mw = 0 to Mw = 7.4 into to
    # the event table. Those events need the right seismic moment.
    fake_event_table = earthquake_event_table[1,]
    fake_event_table$Mw = 0
    fake_event_table$slip = M0_2_Mw(0, inverse=TRUE)/(fake_event_table$area *1e+06 * 3e+10)
    fake_cond_prob = 1
    for(mwtmp in seq(0.1, Mw_min, by=0.1)){
        tmp_event_table = earthquake_event_table[1,]
        tmp_event_table$Mw = mwtmp
        tmp_event_table$slip = M0_2_Mw(mwtmp, inverse=TRUE)/(tmp_event_table$area * 1e+06 * 3e+10)
        fake_event_table = rbind(fake_event_table, tmp_event_table)
        fake_cond_prob = c(fake_cond_prob, 1)
    }

    fake_event_table = rbind(fake_event_table, earthquake_event_table)
    fake_cond_prob = c(fake_cond_prob, event_conditional_probabilities)

    rate_functionC = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate,
        slip_rate_prob = slip_rate_prob,
        b = b,
        b_prob = b_prob,
        Mw_min = 0,
        Mw_min_prob = 1,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = fake_event_table,
        event_conditional_probabilities = fake_cond_prob,
        account_for_moment_below_mwmin=FALSE)

    freq_gt9C = 1/rate_functionC(9.0)
    expect_that(abs(freq_gt9C - freq_gt9B) < 0.1, is_true())

    #
    # Compute some quantiles representing uncertainty in the rate function
    #
    rate_gt9_quantiles = rate_function(9.0, quantiles=seq(0.1, 0.9, by=0.1))
    # The mean of these values should be close to the raw value of rate_function, which
    # is based on the weighted mean of the logic-tree rates
    expect_that( abs(mean(rate_gt9_quantiles) - rate_function(9.0)) < 1e-04, is_true())

    # Back-calculate slip on each fault
    event_rate = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw-dMw/2) - 
        rate_function(earthquake_event_table$Mw + dMw/2))
    event_longterm_slip = earthquake_event_table$slip * event_rate
    unit_source_longterm_slip = 
        rep(0, length(unit_source_summary_statistics[,1]))

    event_longterm_slip_area = event_longterm_slip * earthquake_event_table$area
    unit_source_longterm_slip_area = 
        rep(0, length(unit_source_summary_statistics[,1]))

    unit_source_area = unit_source_summary_statistics$length * 
        unit_source_summary_statistics$width

    # Loop over all events
    for(ee in 1:length(event_longterm_slip)){
        # Find the unit sources on this event
        unit_sources = strsplit(as.character(
            earthquake_event_table$event_index_string[ee]), '-')[[1]]
        unit_sources = as.numeric(unit_sources)

        # Add the slip to all the involved unit sources
        unit_source_longterm_slip[unit_sources] = 
            unit_source_longterm_slip[unit_sources] + event_longterm_slip[ee]

        # Add the seismic moment to all the involved sources
        unit_source_longterm_slip_area[unit_sources] = 
            unit_source_longterm_slip_area[unit_sources] + 
            event_longterm_slip_area[ee]
    }
    #m1 = matrix(unit_source_longterm_slip, ncol=6,nrow=7,byrow=TRUE)
    #m2 = matrix(unit_source_longterm_slip_area, ncol=6,nrow=7,byrow=TRUE)

    # Check that the mean slip is close to the desired value
    # Why is it not identical? -- because our relations ensure that:
    # (long-term-slip x area) = sum(earthquake_slip x area)
    #, which is not the same as ensuring the mean slip itself is identical
    # [unless the mean is weighted by area]
    theoretical_slip = sum(slip_rate*slip_rate_prob)
    back_calc_slip = mean(unit_source_longterm_slip)
    err = abs(back_calc_slip - theoretical_slip)/theoretical_slip
    expect_that(err < 0.01, is_true())

    weighted_slip = weighted.mean(unit_source_longterm_slip, unit_source_area)
    err = abs(sum(weighted_slip) - theoretical_slip)/theoretical_slip
    expect_that(err < 1.0e-08, is_true())


    ###########################################################################
    #
    # Test 2 -- check that the weighting is occurring properly
    #
    ###########################################################################   
    slip_rate = 55.00/1000 # m/year
    slip_rate_prob =  1.0
    b = c(0.7, 1.2)
    b_prob = c(0.3, 0.7)
    Mw_min = 7.5
    Mw_min_prob = 1
    Mw_max = 9.40
    Mw_max_prob = 1

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

    # Compute rate functions individually and confirm their weightings
    rate_function_b0.7 = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate,
        slip_rate_prob = slip_rate_prob,
        b = b[1],
        b_prob = 1.0,
        Mw_min = Mw_min,
        Mw_min_prob = Mw_min_prob,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = earthquake_event_table,
        event_conditional_probabilities = event_conditional_probabilities)

    rate_function_b1.2 = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate,
        slip_rate_prob = slip_rate_prob,
        b = b[2],
        b_prob = 1.0,
        Mw_min = Mw_min,
        Mw_min_prob = Mw_min_prob,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = earthquake_event_table,
        event_conditional_probabilities = event_conditional_probabilities)

    Mws = seq(7.5, 9.5, by=0.1)
    r1 = rate_function(Mws)
    r2 = rate_function_b0.7(Mws)
    r3 = rate_function_b1.2(Mws)

    # This should be zero if the weighing is being correctly done
    err = r1 - (b_prob[1]*r2 + b_prob[2]*r3)

    expect_that( isTRUE(all.equal(err, err*0.0)), is_true())


    #####################################################################
    #
    # Test 3: Check that weighting of logic-tree probabilities with data is ok
    #
    #####################################################################
    #
    # Idea: Provide strongly varying slip rates. Give data that only agrees with one or the other
    # Check that the final averaged rate curve reflects this

    slip_rate = c(0.1, 100)/1000 # m/year
    slip_rate_prob =  c(0.5, 0.5)
    b = 1
    b_prob = 1
    Mw_min = 7.5
    Mw_min_prob = 1
    Mw_max = 9.40
    Mw_max_prob = 1

    # The following data is much more consistent with 100mm of slip/year, vs
    # 0.1mm
    Mw_count_duration = c(7.6, 3, 50)

    # Slip rate = 0.1mm/year alone
    rate_function1 = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate[1],
        slip_rate_prob = 1,
        b = b,
        b_prob = b_prob,
        Mw_min = Mw_min,
        Mw_min_prob = Mw_min_prob,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = earthquake_event_table,
        event_conditional_probabilities = event_conditional_probabilities)
    
    # Slip rate = 100mm/year alone
    rate_function2 = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = slip_rate[2],
        slip_rate_prob = 1,
        b = b,
        b_prob = b_prob,
        Mw_min = Mw_min,
        Mw_min_prob = Mw_min_prob,
        Mw_max = Mw_max,
        Mw_max_prob = Mw_max_prob,
        sourcezone_total_area = sourcezone_total_area,
        event_table = earthquake_event_table,
        event_conditional_probabilities = event_conditional_probabilities)

    # Weighted combination of both models, adjusted for the data
    rate_function3 = rate_of_earthquakes_greater_than_Mw_function(
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
        event_conditional_probabilities = event_conditional_probabilities,
        update_logic_tree_weights_with_data=TRUE,
        Mw_count_duration=Mw_count_duration)

    # Given the data, rate_function3 should be basically equal to rate_function2
    mws = seq(7.5, 9.399, len=100)
    expect_that(
        isTRUE(all.equal(rate_function3(mws), rate_function2(mws), tol=1.0e-05)), 
        is_true())

})
