# Compute the rates for all source-zones, and assign nominal 'individual event' rates
# to each tsunami event in the synthetic catalogue. 
#
# Note the 'individual event' rates are not by themselves particularly
# meaningful [they simply involve distributing the rate for each Mw value over
# all events with that Mw]. For example, the nominal 'individual event' rates
# are inversely related to the number of events in the magnitude category, which is
# somewhat arbitrary [especially for stochastic slip]. 
library(rptha)
config = new.env()
source('config.R', local=config)

#
# INPUTS
# 

sourcezone_parameter_file = config$sourcezone_parameter_file 
sourcezone_parameters = read.csv(sourcezone_parameter_file, stringsAsFactors=FALSE)

source_names = sourcezone_parameters$sourcename

# Never allow Mw_max to be greater than this
MAXIMUM_ALLOWED_MW_MAX = config$MAXIMUM_ALLOWED_MW_MAX  # 9.8
MINIMUM_ALLOWED_MW_MAX = config$MINIMUM_ALLOWED_MW_MAX  # 7.65

# Increment between Mw values in the earthquake_events table. We will check
# that the table holds the same value
dMw = config$dMw # 0.1

# Only assign non-zero probability to earthquakes with Mw greater than this
MW_MIN = config$MW_MIN # 7.45

# We ensure that (Mw_max >= maximum_observed_mw + mw_observed_perturbation)
# This ensures that no logic-tree curve assigns zero probability to the largest
# observed event
mw_observed_perturbation = config$mw_observed_perturbation #0.05

# Truncated or 'characteristic' Gutenberg Richter model
Mw_frequency_dists = config$Mw_frequency_distribution_types #c('truncated_gutenberg_richter', 'characteristic_gutenberg_richter')
Mw_frequency_dists_p = config$Mw_frequency_distribution_weights #c(0.7, 0.3)

# Interpolate logic-tree parameter variation over this many values, all assumed to have
# equal rate. For example, if we provide source_coupling = c(0.1, 0.2, 0.7), then the
# actual coupling values will be "approx(source_coupling, n=nbins)$y"
nbins = config$logic_tree_parameter_subsampling_factor # 9

#
# END INPUTS
#


# Get environment we can use to provide spatially variable convergence information
source('make_spatially_variable_source_zone_convergence_rates.R')
bird2003_env = event_conditional_probability_bird2003_factory(
    return_environment=TRUE)


#
# Function to evaluate the rates for a given source-zone. This function returns
# it's environment, so we have easy access to key variables
#
source_rate_environment_fun<-function(source_name, i, write_rates_to_event_table=FALSE){

    #
    # Coupling
    #
    source_coupling   = sourcezone_parameters[i, c('cmin', 'cpref', 'cmax')]
    source_coupling = approx(as.numeric(source_coupling), n=nbins)$y
    source_coupling_p = rep(1, length(source_coupling))/length(source_coupling)

    #
    # Gutenberg Richter b-value
    #
    source_b   = sourcezone_parameters[i, c('bmin', 'bpref', 'bmax')]
    source_b = approx(as.numeric(source_b), n=nbins)$y
    source_b_p = rep(1, length(source_b))/length(source_b)

    #
    # Mw_max
    #
    unit_source_areas = bird2003_env$unit_source_tables[[source_name]]$length * 
        bird2003_env$unit_source_tables[[source_name]]$width 

    source_area = sum(unit_source_areas)

    min_mw_max = max(sourcezone_parameters$mw_max_observed[i] + mw_observed_perturbation,
        MINIMUM_ALLOWED_MW_MAX)

    source_Mw_max = c(
        ## Largest observed plus a small value,
        #min_mw_max,
        # Middle Mw
        #0.5*(Mw_2_rupture_size_inverse(source_area, CI_sd=0) + 
        #    min_mw_max),
        # Another middle Mw
        Mw_2_rupture_size_inverse(source_area/2, CI_sd=0),
        # Another middle Mw
        Mw_2_rupture_size_inverse(source_area, CI_sd=0),
        # Upper mw [Strasser + 1SD]
        Mw_2_rupture_size_inverse(source_area, CI_sd=-1 ) )
    #
    # Simple test -- all weight on full source rupture area 
    #source_Mw_max = c(Mw_2_rupture_size_inverse(source_area, CI_sd=0), Mw_2_rupture_size_inverse(source_area, CI_sd=0)+0.01)
    #
 
    # Ensure ordered
    source_Mw_max = sort(source_Mw_max)

    # Check it is correctly ordered (of course!)
    stopifnot(all(diff(source_Mw_max) > 0))
    # Ensure all Mw meet out constraints
    source_Mw_max = pmax(source_Mw_max, min_mw_max)
    source_Mw_max = pmin(source_Mw_max, MAXIMUM_ALLOWED_MW_MAX)

    # Interpolate
    source_Mw_max = approx(source_Mw_max, n=nbins)$y

    # Assign equal probabilities to all
    source_Mw_max_p = rep(1, length(source_Mw_max) )/length(source_Mw_max)

    #
    # Tectonic convergence rate
    #

    if(sourcezone_parameters$use_bird_convergence[i] == 1){

        #
        # Idea: If plate convergence vector is between
        # "+-config$rake_deviation_thrust_events" of pure thrust,
        # then use the raw vector. Otherwise, project it onto the nearest
        # within that range
        #
        # This is consistent with our use of data, which extracts earthquakes
        # with rake being within some deviation of pure thrust.
        #
        
        # Shorthand divergent and right-lateral velocity 
        div_vec = pmax(0, -bird2003_env$unit_source_tables[[source_name]]$bird_vel_div)
        bvrl =  bird2003_env$unit_source_tables[[source_name]]$bird_vel_rl

        # Limit lateral component that we consider, based on the permitted rake
        # deviation from pure thrust
        allowed_rake_deviation_radians = config$rake_deviation_thrust_events / 180 * pi
        rl_vec = sign(bvrl) *pmin(abs(bvrl), div_vec*allowed_rake_deviation_radians)

        source_slip = weighted.mean(
            # Convergent slip
            #x= pmax(0, -bird2003_env$unit_source_tables[[source_name]]$bird_vel_div), 
            x = sqrt(div_vec**2 + rl_vec**2),
            # Weighted by area
            w = unit_source_areas)

    }else{
    
        source_slip = sourcezone_parameters$tectonic_slip[i]

    }

    # Account for non-zero dip, and convert from mm/year to m/year
    mean_dip = mean_angle(bird2003_env$unit_source_tables[[source_name]]$dip)
    cos_dip = cos(pi*mean_dip/180)
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
    stopifnot( all( (diff(event_table$Mw) == 0) | ( abs(diff(event_table$Mw) - dMw) < 1.0e-12 ) ) )

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


    if(write_rates_to_event_table){

        # Onto the event table nc file
        fid = nc_open(event_table_file, readunlim=FALSE, write=TRUE)
        ncvar_put(fid, 'rate_annual', event_rates)
        ncvar_put(fid, 'rate_annual_upper_ci', event_rates_upper)
        ncvar_put(fid, 'rate_annual_lower_ci', event_rates_lower)
        nc_close(fid)

        # Onto the uniform slip table nc file
        event_table_fileB = paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc')
        fid = nc_open(event_table_fileB, readunlim=FALSE, write=TRUE)
        ncvar_put(fid, 'event_rate_annual', event_rates)
        ncvar_put(fid, 'event_rate_annual_upper_ci', event_rates_upper)
        ncvar_put(fid, 'event_rate_annual_lower_ci', event_rates_lower)
        nc_close(fid)


        # Onto the stochastic slip table nc file
        event_table_fileC = paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_',
            source_name, '.nc')
        fid = nc_open(event_table_fileC, readunlim=FALSE, write=TRUE)

        # Index corresponding to uniform slip row
        event_uniform_event_row = ncvar_get(fid, 'event_uniform_event_row') 
        # Number of events corresponding to event row
        nevents = table(event_uniform_event_row)
        names_nevents = as.numeric(names(nevents))
        stopifnot(all(names_nevents == 1:length(event_rates)))
        # Make an array giving the number of events matching the uniform_event_row, for every stochastic event
        nevents_broad = nevents[match(event_uniform_event_row, names_nevents )]

        ncvar_put(fid, 'event_rate_annual', event_rates[event_uniform_event_row]/nevents_broad)
        ncvar_put(fid, 'event_rate_annual_upper_ci', event_rates_upper[event_uniform_event_row]/nevents_broad)
        ncvar_put(fid, 'event_rate_annual_lower_ci', event_rates_lower[event_uniform_event_row]/nevents_broad)
        nc_close(fid)
    }


    return(environment())

}
    

source_envs = vector(mode='list', length=length(source_names))
names(source_envs) = source_names

for(i in 1:length(source_names)){

    source_name = source_names[i]
    source_envs[[i]] = source_rate_environment_fun(source_name, i, 
        write_rates_to_event_table=TRUE)

}



