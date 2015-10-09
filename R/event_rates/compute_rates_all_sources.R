# Given an appropriately formatted table with the sourcezone parameters, and a
# discrete-source RDS file, this code computes all earthquake events on each
# sourcezone, and assigns probabilities to each based on the sourcezone
# parameters.

library(rptha)

# Make a function to do the main computations in serial
compute_rate_curve<-function(discrete_source, Mw_min, Mw_max, slip_rate,
    b, Mw_min_prob, Mw_max_prob, slip_rate_prob, b_prob, dMw = 0.1){

    # Trick so we get something even if it fails
    on.exit(return(environment()))

    # Get summary statistics for unit sources 
    unit_source_summary_statistics = 
        discretized_source_summary_statistics(discrete_source,
            approx_dx=5000, approx_dy = 5000)

    # Get a table with all earthquake events between min(Mw_min) and max(Mw_max)
    earthquake_event_table = get_all_earthquake_events(discrete_source, 
        unit_source_summary_statistics, Mmin = min(Mw_min), Mmax = max(Mw_max), 
        dMw = dMw)

    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        earthquake_event_table,
        conditional_probability_model = 'inverse_slip')

    # Get total area of source-zone in km^2. This gives the down-dip area.
    sourcezone_total_area = sum(unit_source_summary_statistics$length * 
        unit_source_summary_statistics$width)

    # Get the rate function, integrating out the logic tree
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


    # Back-calculate slip on each fault
    event_rate = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw-dMw/2) - rate_function(earthquake_event_table$Mw + dMw/2))
    event_longterm_slip = earthquake_event_table$slip * event_rate
    event_longterm_slip_area = event_longterm_slip * earthquake_event_table$area

    unit_source_longterm_slip = rep(0, length(unit_source_summary_statistics[,1]))
    unit_source_longterm_slip_area = rep(0, length(unit_source_summary_statistics[,1]))

    unit_source_area = unit_source_summary_statistics$length*unit_source_summary_statistics$width

    # Loop over all events
    for(ee in 1:length(event_longterm_slip)){
        # Find the unit sources on this event
        unit_sources = strsplit(as.character(earthquake_event_table$event_index_string[ee]), '-')[[1]]
        unit_sources = as.numeric(unit_sources)

        # Add the slip to all the involved unit sources
        unit_source_longterm_slip[unit_sources] = unit_source_longterm_slip[unit_sources] + 
            event_longterm_slip[ee]
        # Add the seismic moment to all the involved sources
        unit_source_longterm_slip_area[unit_sources] = unit_source_longterm_slip_area[unit_sources] + 
            event_longterm_slip[ee]*unit_source_area[unit_sources]
    }

    theoretical_slip = sum(slip_rate*slip_rate_prob)

    # The area-weighted-mean of the unit_source_longterm_slip should equal the theoretical slip,
    # by construction
    wm_slip = weighted.mean(unit_source_longterm_slip, unit_source_area)
    stopifnot(abs( wm_slip - theoretical_slip)/theoretical_slip < 0.0001)

    # The mean slip should be fairly close to this, if not exactly equal
    m_slip = mean(unit_source_longterm_slip)
    stopifnot(abs(m_slip - theoretical_slip)/theoretical_slip < 0.1)
}

##############################################################################
#
# MAIN CODE HERE
#
#

all_discrete_sourcezones = readRDS(
    '../SOURCE_CONTOURS_2_UNIT_SOURCES/all_discretized_sources.RDS')

## Read the key source-zone parameters
#sourcezone_par = read.csv('Parameters/GAR_revised_pars.txt', stringsAsFactors=FALSE)
sourcezone_par = read.csv('Parameters/GAR_revised_pars_B.txt', stringsAsFactors=FALSE)
names(sourcezone_par)[1] = 'source'

sourcezones = sourcezone_par$source

# Store output
sourcezone_output = list()

# Make a function to compute the sourcezone rates, for convenient execution in
# parallel (or serial)
parallel_f<-function(sourcezone){

    data_row = match(sourcezone, sourcezone_par[,1])

    if(is.na(data_row)){
        return(NA)
    }

    discrete_source = all_discrete_sourcezones[[sourcezone]]

    data_line = sourcezone_par[match(sourcezone, sourcezone_par[,1]), ]

    Mw_min = 7.5
    Mw_min_prob = 1.0

    Mw_max = as.numeric(data_line[c('Mmax.pref', 'Mmax.min', 'Mmax.max')])

    # Truncate Mw_max so 7.6 is the minimum Mw_max [so it is higher than our
    # 7.5 Mw_min]
    if(any(Mw_max < (Mw_min + 0.1))){
        print('Forcing Mw_max >= Mw_min + 0.1')
        Mw_max = pmax(Mw_max, Mw_min + 0.1)
    }

    Mw_max = round(Mw_max, 1) # Need spacing comparable to Mw
    Mw_max_prob = c(0.6, 0.2, 0.2) # Changed from earlier GAR

    b = as.numeric(data_line[c('b.pref', 'b.min', 'b.max')])
    b_prob = c(0.6, 0.2, 0.2) # Follows earlier GAR

    # long-term slip rate in m/year
    slip_rate = as.numeric(data_line[c('Slip.pref', 'Slip.min', 'Slip.max')])/1000
    slip_rate_prob = c(0.6, 0.2, 0.2) # Follows earlier GAR


    output = compute_rate_curve(discrete_source, 
        Mw_min, Mw_max, slip_rate, b, Mw_min_prob, Mw_max_prob, slip_rate_prob, 
        b_prob)
    return(output)
}


# Run in parallel, or serial, over every source zone
parallel_run = TRUE

if(!parallel_run){
    for(sourcezone in sourcezones){
        print(sourcezone)
        sourcezone_output[[sourcezone]] = parallel_f(sourcezone)
    }
}else{
    library(parallel)
    sourcezone_output = mclapply(as.list(sourcezones), parallel_f, mc.cores=12, mc.preschedule=FALSE)
}

# Ensure each entry in the list has a good name
names(sourcezone_output) = sourcezones

# Save the outputs
#saveRDS(sourcezone_output, 'sourcezone_events_and_rates.RDS')
saveRDS(sourcezone_output, 'sourcezone_events_and_rates_B.RDS')
