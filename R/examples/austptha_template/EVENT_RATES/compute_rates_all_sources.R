# Get environment we can use to provide spatially variable convergence information
source('make_spatially_variable_source_zone_convergence_rates.R')
bird2003_env = event_conditional_probability_bird2003_factory(
    return_environment=TRUE)

sourcezone_parameter_file = '../DATA/SOURCEZONE_PARAMETERS/sourcezone_parameters.csv'
sourcezone_parameters = read.csv(sourcezone_parameter_file, stringsAsFactors=FALSE)

source_names = sourcezone_parameters$sourcename

# Never allow Mw_max to be greater than this
MAXIMUM_ALLOWED_MW_MAX = 9.6
MINIMUM_ALLOWED_MW_MAX = 7.6

# Only assign non-zero probability to earthquakes with Mw greater than this
MW_MIN = 7.5

# Increment between Mw values in the earthquake_events table. We will check
# that the table holds the same value
dMw = 0.1

# Truncated or 'characteristic' Gutenberg Richter model
Mw_frequency_dists = c('truncated_gutenberg_richter', 'characteristic_gutenberg_richter')
Mw_frequency_dists_p = c(0.5, 0.5)


#
# Function to evaluate the rates for a given source-zone. This function returns
# it's environment, so we have easy access to key variables
source_rate_environment_fun<-function(source_name, i){

    #
    # Coupling
    #
    source_coupling   = sourcezone_parameters[i, c('cmin', 'cpref', 'cmax')]
    source_coupling_p = sourcezone_parameters[i, c('cmin_p', 'cpref_p', 'cmax_p')]

    #
    # Gutenberg Richter b-value
    #
    source_b   = sourcezone_parameters[i, c('bmin', 'bpref', 'bmax')]
    source_b_p = sourcezone_parameters[i, c('bmin_p', 'bpref_p', 'bmax_p')]

    #
    # Mw_max
    #
    unit_source_areas = bird2003_env$unit_source_tables[[source_name]]$length * 
        bird2003_env$unit_source_tables[[source_name]]$width 

    source_area = sum(unit_source_areas)

    source_Mw_max = c(
        # Largest observed,
        sourcezone_parameters$mw_max_observed[i],
        # Lower mw
        Mw_2_rupture_size_inverse(source_area, CI_sd=1 ),
        # Middle mw
        Mw_2_rupture_size_inverse(source_area, CI_sd=0  ),
        # Upper mw
        Mw_2_rupture_size_inverse(source_area, CI_sd=-1 ))

    # Check it is correctly ordered (of course!)
    stopifnot(all(diff(source_Mw_max) > 0))
    # Ensure all Mw meet out constraints
    source_Mw_max = pmax(source_Mw_max, sourcezone_parameters$mw_max_observed[i])
    source_Mw_max = pmax(source_Mw_max, MINIMUM_ALLOWED_MW_MAX)
    source_Mw_max = pmin(source_Mw_max, MAXIMUM_ALLOWED_MW_MAX)
    # Assign equal probabilities to all
    source_Mw_max_p = rep(1, length(source_Mw_max) ) / length(source_Mw_max)

    #
    # Tectonic convergence rate
    #
    mean_dip = mean(bird2003_env$unit_source_tables[[source_name]]$dip)
    cos_dip = cos(2*pi*mean_dip/180)

    if(sourcezone_parameters$use_bird_convergence[i] == 1){


        source_slip = weighted.mean(
            # Convergent slip
            x= pmax(0, -bird2003_env$unit_source_tables[[source_name]]$bird_vel_div), 
            # Weighted by area
            w = unit_source_areas)

    }else{
    
        source_slip = sourcezone_parameters$tectonic_slip[i]

    }

    # Account for non-zero dip, and convert from mm/year to m/year
    source_slip = source_slip/cos_dip * 1/1000

    #
    # Event table
    #
    event_table_file = paste0('../SOURCE_ZONES/', source_name, 
        '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_', 
        source_name, '.nc')
    event_table = read_table_from_netcdf(event_table_file)

    # Round away any finite-precision issues in the netcdf file
    # Our increments are by dMw
    event_table$Mw = round(event_table$Mw, 3)
    stopifnot( all(event_table$Mw == sort(event_table$Mw)) )
    # Check we didn't destroy the table by rouding!
    stopifnot( all(diff(event_table$Mw) == 0) | (abs(diff(event_table$Mw) - dMw) < 1.0e-12) )

    #
    # Event conditional probabilities
    #
    if(sourcezone_parameters$use_bird_convergence[i] == 1){
        conditional_probability_model = 
            bird2003_env$make_conditional_probability_function_uniform_slip(source_name)
    }else{
        conditional_probability_model = 'inverse_slip'
    }

    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        event_table, conditional_probability_model = conditional_probability_model)    

    #
    # Build rate function
    #
    mw_rate_function = rate_of_earthquakes_greater_than_Mw_function(
        slip_rate = as.numeric(source_slip * source_coupling),
        slip_rate_prob = as.numeric(source_coupling_p),
        b = as.numeric(source_b),
        b_prob = source_b_p,
        Mw_min = MW_MIN,
        Mw_min_prob = 1,
        Mw_max = as.numeric(source_Mw_max),
        Mw_max_prob = source_Mw_max_p,
        sourcezone_total_area = source_area,
        event_table = event_table,
        event_conditional_probabilities = event_conditional_probabilities,
        computational_increment = 0.01,
        Mw_frequency_distribution = Mw_frequency_dists,
        Mw_frequency_distribution_prob = Mw_frequency_dists_p,
        account_for_moment_below_mwmin = TRUE
        )


    #
    # Nominal rates on all events
    #

    event_rates = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw -dMw/2) - 
        mw_rate_function(event_table$Mw + dMw/2) )

    event_rates_upper = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw -dMw/2, quantiles=0.975) - 
        mw_rate_function(event_table$Mw + dMw/2, quantiles=0.975) )

    event_rates_lower = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw -dMw/2, quantiles=0.025) - 
        mw_rate_function(event_table$Mw + dMw/2, quantiles=0.025) )

    return(environment())

}
    

source_envs = vector(mode='list', length=length(source_names))
names(source_envs) = source_names
for(i in 1:length(source_names)){

    source_name = source_names[i]
    source_envs[[i]] = source_rate_environment_fun(source_name, i)

}



