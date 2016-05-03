# Given an appropriately formatted table with the sourcezone parameters, and a
# discrete-source RDS file, this code computes all earthquake events on each
# sourcezone, and assigns probabilities to each based on the sourcezone
# parameters.

library(rptha)

# A function to do our computations on a single discrete_source
compute_rate_curve<-function(
    discrete_source, 
    Mw_min, 
    Mw_max, 
    slip_rate, 
    b, 
    Mw_freq_distribution_type, 
    Mw_min_prob, 
    Mw_max_prob, 
    slip_rate_prob, 
    b_prob, 
    Mw_freq_distribution_type_prob,
    conv_slip_factor, 
    dMw = 0.1, 
    update_logic_tree_weights_with_data=FALSE, 
    Mw_count_duration=c(NA,NA,NA), 
    discrete_source_data_dir=NULL,
    use_variable_slip=FALSE, 
    sourcezone=NULL){
    
    # Ensure we get some output even if the code fails
    success = FALSE
    on.exit(return(environment()))

    # FIXME: These lines of code should be optionally replaced by reads from
    # files (since the unit source summary statistics and earthquake event
    # table are required elsewhere, when combining unit sources).
    # OR, alternatively, we could keep them, produce the output, and then read
    # it in the 'combine events' code
    # The second option is probably better, because it is highly convenient to
    # adapt Mw_min and Mw_max and dMw in this routine while checking that the
    # eq probabilities are reasonable.
    if(is.null(discrete_source_data_dir)){
        # Get summary statistics for unit sources 
        unit_source_summary_statistics = 
            discretized_source_summary_statistics(discrete_source,
                approx_dx=4000, approx_dy = 4000)

        # Get a table with all earthquake events between min(Mw_min) and
        # max(Mw_max)
        earthquake_event_table = get_all_earthquake_events(
            discrete_source=discrete_source, 
            unit_source_statistics=unit_source_summary_statistics, 
            Mmin = min(Mw_min), Mmax = max(Mw_max), 
            dMw = dMw)
    }else{
        unit_source_summary_statistics = readRDS(
            paste0(discrete_source_data_dir, '/', 
                basename(discrete_source_data_dir), 
                '_unit_source_statistics.RDS'))

        earthquake_event_table = readRDS(
            paste0(discrete_source_data_dir, '/', 
                basename(discrete_source_data_dir),
                '_earthquake_events_table.RDS'))

    }

    if( (!is.null(sourcezone)) && (sourcezone != 'gloria')){
        # Rescale slip by 1/cos(dip) -- except for Gloria.
        # Just use a single dip.
        # Motivation -- our long term slip rates are measured in a horizontal
        # plane, and for slip up-dip to balance this, we need to divide by
        # cos(dip).
        mean_dip = weighted.mean(
            unit_source_summary_statistics$dip, 
            weights=(unit_source_summary_statistics$length * 
                unit_source_summary_statistics$width))
        slip_rescale = 1/cos(mean_dip/180 * pi)
        msg = 'Warning: Rescaling long term slip rate by 1/cos(av_dip) = '
        print(paste0(msg, round(slip_rescale,4)))
        slip_rate = slip_rate * slip_rescale
    }


    if(use_variable_slip){

        # Make function to set the conditional probabilities based on the
        # nearby convergent slip and area. The function must take a subset of
        # the event table where all Mw are identical

        cond_prob<-function(event_table_subset){
            # For each event in the table, find the area x long-term slip
            # perpendicular to the trench.
            # Make the conditional probability of the event proportional to
            # this.
            event_slip_areas = rep(NA, length(event_table_subset[,1]))

            bird_name_match = which(bird_slip$src_nms == sourcezone)
            bird_slip_sourcezone = bird_slip[bird_name_match,]

            for(i in 1:length(event_slip_areas)){
                # Get unit sources included
                us_string = event_table_subset$event_index_string[i]
                us_included = as.numeric(
                    strsplit(as.character(us_string),'-')[[1]])
               
                # Get their area 
                area = unit_source_summary_statistics$length[us_included] *
                    unit_source_summary_statistics$width[us_included]

                # Get their long-term slip rates, based on the strike index,
                # and the shapefile data that maps Bird's slip rates to our
                # along-trench unit sources.
                strike_indices = 
                    unit_source_summary_statistics$alongstrike_number[us_included]
                match_strike_indices = 
                    match(strike_indices, bird_slip_sourcezone$strk_nd)

                slip_rates = bird_slip_sourcezone$Div_vel[match_strike_indices]
                # Ensure convergent slip is treated as positive, and divergent
                # slip gives 0 activity
                slip_rates = pmax(-slip_rates, 0)

                event_slip_areas[i] = sum(area*slip_rates)
            }

            conditional_prob = event_slip_areas/sum(event_slip_areas)

            return(conditional_prob)
        }

        # Compute the conditional probabilities with the above function
        event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
            earthquake_event_table,
            conditional_probability_model = cond_prob)

    }else{
        event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
            earthquake_event_table,
            conditional_probability_model = 'inverse_slip')
    }

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
        event_conditional_probabilities=event_conditional_probabilities,
        Mw_frequency_distribution=Mw_freq_distribution_type,
        Mw_frequency_distribution_prob=Mw_freq_distribution_type_prob,
        update_logic_tree_weights_with_data=update_logic_tree_weights_with_data,
        Mw_count_duration=Mw_count_duration,
        account_for_moment_below_mwmin=TRUE)

    # Do computations with deterministic Mw-max, so we can look at the impact
    # of logic-tree Mw-max
    extra_rate_functions = list()
    for(i in 1:3){
        
        if(i==1) local_mw_max = min(Mw_max[1:3])
        if(i==2) local_mw_max = median(Mw_max[1:3])
        if(i==3) local_mw_max = max(Mw_max[1:3])

        local_rate_function = rate_of_earthquakes_greater_than_Mw_function(
            slip_rate = slip_rate,
            slip_rate_prob = slip_rate_prob,
            b = b,
            b_prob = b_prob,
            Mw_min = Mw_min,
            Mw_min_prob = Mw_min_prob,
            Mw_max = local_mw_max,
            Mw_max_prob = 1,
            sourcezone_total_area = sourcezone_total_area,
            event_table = earthquake_event_table,
            event_conditional_probabilities = event_conditional_probabilities,
            Mw_frequency_distribution=Mw_freq_distribution_type,
            Mw_frequency_distribution_prob=Mw_freq_distribution_type_prob,
            update_logic_tree_weights_with_data=update_logic_tree_weights_with_data,
            Mw_count_duration=Mw_count_duration,
            account_for_moment_below_mwmin=TRUE)
        extra_rate_functions[[i]] = local_rate_function
    }

    # Back-calculate slip on each fault
    # Assign the rate of events between [Mw -dMw/2, Mw + dMw/2] to the Mw events
    event_rate = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw - dMw/2) - 
         rate_function(earthquake_event_table$Mw + dMw/2))
    
    event_rate_975 = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw - dMw/2, quantiles=0.975) - 
         rate_function(earthquake_event_table$Mw + dMw/2, quantiles=0.975))
    
    event_rate_025 = event_conditional_probabilities * 
        (rate_function(earthquake_event_table$Mw - dMw/2, quantiles=0.025) - 
         rate_function(earthquake_event_table$Mw + dMw/2, quantiles=0.025))

    event_longterm_slip = earthquake_event_table$slip * event_rate
    event_longterm_slip_area = event_longterm_slip * earthquake_event_table$area

    unit_source_longterm_slip = rep(0, length(unit_source_summary_statistics[,1]))
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
        unit_source_longterm_slip[unit_sources] = unit_source_longterm_slip[unit_sources] + 
            event_longterm_slip[ee]
        # Add the seismic moment to all the involved sources
        unit_source_longterm_slip_area[unit_sources] = 
            unit_source_longterm_slip_area[unit_sources] + 
            event_longterm_slip[ee]*unit_source_area[unit_sources]
    }

    # Get the updated slip rate probabilities. These may have changed from inputs
    # if Bayesian re-weighting of logic-tree branches is applied
    # theoretical_slip = sum(slip_rate * slip_rate_prob) # Original weights
    tmp = rate_function(NA, return_all_logic_tree_branches=TRUE)
    theoretical_slip = sum(tmp$all_par$slip_rate * tmp$all_par_prob)

    # The area-weighted-mean of the unit_source_longterm_slip should equal
    # the theoretical slip
    wm_slip = weighted.mean(unit_source_longterm_slip, unit_source_area)
    stopifnot(abs( wm_slip - theoretical_slip)/theoretical_slip < 0.0001)

    # The mean slip should be fairly close to this, if not exactly equal
    # because of our weighting by unit source area x local long-term slip
    m_slip = mean(unit_source_longterm_slip)
    stopifnot(abs(m_slip - theoretical_slip)/theoretical_slip < 0.1)

    # Keep this to easily check if there were no errors
    success = TRUE
}

##############################################################################
#
# MAIN CODE HERE
#
#

if(interactive()){
    if(!exists('parameter_file')){
        stop('variable parameter_file is needed to run this script')
    }
}else{
    parameter_file = commandArgs(trailingOnly=TRUE)[1]
}
output_dir = gsub('.txt', '', basename(parameter_file))

sourcezone_par = read.csv(parameter_file, stringsAsFactors=FALSE)
names(sourcezone_par)[1] = 'source'

all_discrete_sourcezones = readRDS(
    '../SOURCE_CONTOURS_2_UNIT_SOURCES/all_discretized_sources.RDS')

#counts_mw75_all_sourcezones = read.csv('count_events_gt_Mw75.csv', stringsAsFactors=FALSE)
counts_mw75_all_sourcezones = read.csv('count_thrust_events_gt_Mw75.csv', stringsAsFactors=FALSE)

# Get a shapefile containing info on long-term convergent slip rates at the top
# edge of each 'along-trench' unit source
bird_slip = readOGR('../SOURCE_CONTOURS_2_UNIT_SOURCES/topedge2bird', layer='topedge2bird')

sourcezones = sourcezone_par$source

# Make an output directory
dir.create(output_dir, showWarnings=FALSE)

# Store output
sourcezone_output = list()

# Make a function to compute the sourcezone rates, for convenient execution in
# parallel (or serial)
parallel_f<-function(sourcezone, read_source_statistics_and_earthquake_events=FALSE){

    data_row = match(sourcezone, sourcezone_par[,1])

    if(is.na(data_row)){
        return(NA)
    }

    discrete_source = all_discrete_sourcezones[[sourcezone]]

    data_line = sourcezone_par[match(sourcezone, sourcezone_par[,1]), ]

    Mw_min = 7.5
    Mw_min_prob = 1.0

    # We append Mmax.huge to Mw_max, but give it 0 probability. This allows us
    # to have events covering a large range (convenient if we need to extend
    # the range later)

    if(sum(grepl('Strasser', names(data_line)))>0){
        Mw_max = as.numeric(data_line[c('Mmax.Strasser',  'Mmax.McCaffrey', 
            'Mmax.Strasser1sd', 'Mmax.huge')])
        Mw_max_prob = c(0.4, 0.4, 0.2, 0.) # Changed from earlier GAR as explained above. 
    }else{
        Mw_max = as.numeric(data_line[c('Mmax.pref',  'Mmax.min', 'Mmax.max', 'Mmax.huge')])
        ##Mw_max_prob = c(0.6, 0.2, 0.2, 0.) # Following earlier gar
        Mw_max_prob = c(0.45, 0.1, 0.45, 0) # Revised
    }

    # Truncate Mw_max so 7.6 is the minimum Mw_max [so it is higher than our
    # 7.5 Mw_min]
    if(any(Mw_max < (Mw_min + 0.1))){
        print('Forcing Mw_max >= Mw_min + 0.1')
        Mw_max = pmax(Mw_max, Mw_min + 0.1)
    }

    Mw_max = round(Mw_max, 1) # Need spacing comparable to Mw

    if(FALSE){
        print('WARNING: Converting Mw_max to sequence, with equal weights')
        Mw_max = c(seq(min(Mw_max[1:3]), max(Mw_max[1:3]), len=11), Mw_max[4])
        Mw_max_prob = c( rep(1, length(Mw_max)-1)/(length(Mw_max) - 1), 0)
    }

    b = as.numeric(data_line[c('b.pref', 'b.min', 'b.max')])
    #b_prob = c(0.6, 0.2, 0.2) # Follows earlier GAR
    #b_prob = c(1.0, 0.0, 0.0) # Experimental
    b_prob = c(1/3, 1/3, 1/3) # Experimental
    #b_prob = c(0.2, 0.1, 0.7)
    #b_prob = c(0.1, 0.8, 0.1)

    if(FALSE){
        print('WARNING: Converting b to a sequence with equal weights')
        b = seq(min(b), max(b), len=11)
        b_prob = rep(1,11)/11
    }

    # long-term slip rate in m/year
    slip_rate = as.numeric(data_line[c('Slip.pref', 'Slip.min', 'Slip.max')])/1000 #* 0.5
    #slip_rate_prob = c(0.6, 0.2, 0.2) # Follows earlier GAR
    #slip_rate_prob = c(1.0, 0.0, 0.0) # Experimental 
    slip_rate_prob = c(1/3, 1/3, 1/3) # Follows earlier GAR
    #slip_rate_prob = c(0.7, 0.2, 0.1) # Follows earlier GAR

    if(FALSE){
        print('WARNING: Converting slip_rate to a sequence with equal weights')
        slip_rate = seq(min(slip_rate), max(slip_rate), len=11)
        slip_rate_prob = rep(1,11)/11
    }

    conv_slip_factor = as.numeric(data_line['conv.slip.factor'])
    print('WARNING: Adjusting slip to represent perpendicular component')
    slip_rate = slip_rate*(-conv_slip_factor) 

    use_variable_slip = as.logical(as.numeric(data_line['Use_spatially_variable_slip']))


    Mw_freq_distribution_type = c('characteristic_gutenberg_richter',
        'truncated_gutenberg_richter')
    Mw_freq_distribution_type_prob = c(0.5, 0.5)

    dMw = 0.1
    # Get Mw-count-data
    Mw_count_duration_index = which(counts_mw75_all_sourcezones[,1] == sourcezone)
    if(length(Mw_count_duration_index) != 1){
        stop(paste0('No count-data found for sourcezone ', sourcezone))
    }
    Mw_count_duration = counts_mw75_all_sourcezones[Mw_count_duration_index,2:4]
    # Must be 'numeric' to avoid some data_type issues
    Mw_count_duration = as.numeric(Mw_count_duration)

    if(read_source_statistics_and_earthquake_events){
        # Read statistics from RDS in this directory
        discrete_source_data_dir = paste0(
            '../COMBINE_TSUNAMI_SOURCES/NCI_OUTPUTS2/merged_outputs/', 
            sourcezone)

    }else{
        discrete_source_data_dir=NULL
    }

    output = compute_rate_curve(discrete_source, 
        Mw_min, Mw_max, slip_rate, b, Mw_freq_distribution_type, 
        Mw_min_prob, Mw_max_prob, slip_rate_prob, b_prob, 
        Mw_freq_distribution_type_prob,
        conv_slip_factor, dMw = dMw, 
        update_logic_tree_weights_with_data=TRUE,
        Mw_count_duration=Mw_count_duration,
        discrete_source_data_dir=discrete_source_data_dir,
        use_variable_slip=use_variable_slip,
        sourcezone=sourcezone)
    return(output)
}

parallel_run = TRUE


if(!parallel_run){
    for(sourcezone in sourcezones){
        print(sourcezone)
        sourcezone_output[[sourcezone]] = parallel_f(sourcezone, 
            read_source_statistics_and_earthquake_events=TRUE)
    }
}else{
    library(parallel)
    sourcezone_output = mclapply(as.list(sourcezones), parallel_f,
        read_source_statistics_and_earthquake_events=TRUE,
        mc.cores=12, mc.preschedule=FALSE)
}

names(sourcezone_output) = sourcezones
#saveRDS(sourcezone_output, 'sourcezone_events_and_rates.RDS')
saveRDS(sourcezone_output, paste0(output_dir, '/sourcezone_events_and_rates_B.RDS'))

## HACK -- forgot to write the unit source summary statistics for NCI outputs,
## so I reproduce them here.
#for(sourcename in sourcezones){
#    output_RDS = paste0('../COMBINE_TSUNAMI_SOURCES/NCI_OUTPUTS/outputs/', 
#        sourcename, '/', sourcename, '_unit_source_statistics.RDS')
#
#    try(
#    saveRDS(sourcezone_output[[sourcename]]$unit_source_summary_statistics,
#        output_RDS)
#    )
#}
