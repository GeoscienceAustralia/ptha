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

# Get environment we can use to provide spatially variable convergence
# information
source('make_spatially_variable_source_zone_convergence_rates.R')
bird2003_env = event_conditional_probability_bird2003_factory(
    return_environment=TRUE)


#' Function to evaluate the rates for a given source-zone. This function returns
#' it's environment, so we have easy access to key variables
#'
#' @param sourcezone_parameters_row row from sourcezone_parameters table
#' @return the function environment
#'
source_rate_environment_fun<-function(sourcezone_parameters_row){

    # Store key parameters for easy write-out later
    sourcepar = list()

    sourcepar$sourcezone_parameters_row = sourcezone_parameters_row
    source_name = sourcezone_parameters_row$sourcename
    sourcepar$name = source_name

    # Source area 
    unit_source_areas = bird2003_env$unit_source_tables[[source_name]]$length * 
        bird2003_env$unit_source_tables[[source_name]]$width 

    sourcepar$area = sum(unit_source_areas)

    #
    # Gutenberg Richter b-value
    #
    sourcepar$b = sourcezone_parameters_row[1, c('bmin', 'bpref', 'bmax')]
    sourcepar$b = approx(as.numeric(sourcepar$b), n=nbins)$y
    sourcepar$b_p = rep(1, length(sourcepar$b))/length(sourcepar$b)

    #
    # Coupling
    #
    sourcepar$coupling   = sourcezone_parameters_row[1, c('cmin', 'cpref', 'cmax')]
    sourcepar$coupling = approx(as.numeric(sourcepar$coupling), n=nbins)$y
    sourcepar$coupling_p = rep(1, length(sourcepar$coupling))/length(sourcepar$coupling)

    #
    # Mw_max
    #

    min_mw_max = max(
        sourcezone_parameters_row$mw_max_observed + mw_observed_perturbation,
        MINIMUM_ALLOWED_MW_MAX)

    sourcepar$Mw_max = c(
        ## Largest observed plus a small value,
        #min_mw_max,
        # Middle Mw
        #0.5*(Mw_2_rupture_size_inverse(sourcepar$area, CI_sd=0) + 
        #    min_mw_max),
        # Another middle Mw
        Mw_2_rupture_size_inverse(sourcepar$area/2, CI_sd=0),
        # Another middle Mw
        Mw_2_rupture_size_inverse(sourcepar$area, CI_sd=0),
        # Upper mw [Strasser + 1SD]
        Mw_2_rupture_size_inverse(sourcepar$area, CI_sd=-1 ) )
    #
    # Simple test -- all weight on full source rupture area 
    #sourcepar$Mw_max = c(Mw_2_rupture_size_inverse(sourcepar$area, CI_sd=0),
    #    Mw_2_rupture_size_inverse(sourcepar$area, CI_sd=0)+0.01)
    #

    # Ensure ordered
    sourcepar$Mw_max = sort(sourcepar$Mw_max)

    # Check it is correctly ordered (of course!)
    stopifnot(all(diff(sourcepar$Mw_max) > 0))
    # Ensure all Mw meet out constraints
    sourcepar$Mw_max = pmax(sourcepar$Mw_max, min_mw_max)
    sourcepar$Mw_max = pmin(sourcepar$Mw_max, MAXIMUM_ALLOWED_MW_MAX)

    # Interpolate
    sourcepar$Mw_max = approx(sourcepar$Mw_max, n=nbins)$y

    # Assign equal probabilities to all
    sourcepar$Mw_max_p = rep(1, length(sourcepar$Mw_max) )/length(sourcepar$Mw_max)

    #
    # Tectonic convergence rate
    #

    if(sourcezone_parameters_row$use_bird_convergence == 1){

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
        deg2rad = pi/180
        allowed_rake_deviation_radians = config$rake_deviation_thrust_events * deg2rad
        rl_vec = sign(bvrl) *pmin(abs(bvrl), div_vec*allowed_rake_deviation_radians)

        sourcepar$slip = weighted.mean(
            # Convergent slip
            #x= pmax(0, -bird2003_env$unit_source_tables[[source_name]]$bird_vel_div), 
            x = sqrt(div_vec**2 + rl_vec**2),
            # Weighted by area
            w = unit_source_areas)

    }else{

        sourcepar$slip = sourcezone_parameters_row$tectonic_slip

    }

    # Account for non-zero dip, and convert from mm/year to m/year
    mean_dip = mean_angle(bird2003_env$unit_source_tables[[source_name]]$dip)
    sourcepar$mean_dip = mean_dip
    deg2rad = pi/180
    cos_dip = cos(mean_dip*deg2rad)
    sourcepar$cos_dip = cos_dip

    sourcepar$slip = sourcepar$slip/cos_dip * 1/1000

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
    if(sourcezone_parameters_row$use_bird_convergence == 1){
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
        slip_rate = as.numeric(sourcepar$slip * sourcepar$coupling),
        slip_rate_prob = as.numeric(sourcepar$coupling_p),
        b = as.numeric(sourcepar$b),
        b_prob = sourcepar$b_p,
        Mw_min = MW_MIN,
        Mw_min_prob = 1,
        Mw_max = as.numeric(sourcepar$Mw_max),
        Mw_max_prob = sourcepar$Mw_max_p,
        sourcezone_total_area = sourcepar$area,
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
        (mw_rate_function(event_table$Mw - dMw/2) - 
         mw_rate_function(event_table$Mw + dMw/2) )

    event_rates_upper = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw -dMw/2, quantiles=0.975) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.975) )

    event_rates_lower = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw -dMw/2, quantiles=0.025) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.025) )


    return(environment())

}


#' Put the event rates onto a netcdf file
#'
#' @param source_env result of source_rate_environment_fun
#' @param scale_rate multiple the rate by this constant before putting on file
#' @param add_rate add the (scaled) rate to the values already existing in the
#' file
#'
write_rates_to_event_table<-function(source_env, scale_rate=1.0, 
    add_rate=FALSE){

    # Add function arguments to environment source_env 
    with(source_env, {
        scale_rate = scale_rate
        add_rate = add_rate
    })

    # Write to netcdf file. Easiest from within source_env
    with(source_env,{

        #
        # Function to do slightly modified 'add to netcdf' computation
        #
        # It can put a 'scaled' version of the rate onto the file, and optionally
        # add to the variable already on the file. 
        #
        # This would allow e.g. weighted treatments of segmented and unsegmented
        # models of the source-zone.
        ncvar_put_extra<-function(fid, varname, event_rate, add_rate=add_rate, 
            scale_rate=scale_rate){

            if(add_rate){
                extra_rate = ncvar_get(fid, varname)
            }else{
                extra_rate = event_rate * 0
            }

            ncvar_put(fid, varname, event_rate*scale_rate + extra_rate)

            return(invisible())
        }

        # Put rates onto the event table nc file
        fid = nc_open(event_table_file, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'rate_annual', event_rates)
        ncvar_put_extra(fid, 'rate_annual_upper_ci', event_rates_upper)
        ncvar_put_extra(fid, 'rate_annual_lower_ci', event_rates_lower)
        nc_close(fid)

        # Put rates onto the uniform slip table nc file
        event_table_fileB = paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc')
        fid = nc_open(event_table_fileB, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'event_rate_annual', event_rates)
        ncvar_put_extra(fid, 'event_rate_annual_upper_ci', event_rates_upper)
        ncvar_put_extra(fid, 'event_rate_annual_lower_ci', event_rates_lower)
        nc_close(fid)

        # Put rates onto the stochastic slip table nc file
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
        # Make an array giving the number of events matching the
        # uniform_event_row, for every stochastic event
        nevents_broad = nevents[match(event_uniform_event_row, names_nevents )]

        ncvar_put_extra(fid, 'event_rate_annual', 
            event_rates[event_uniform_event_row]/nevents_broad)
        ncvar_put_extra(fid, 'event_rate_annual_upper_ci', 
            event_rates_upper[event_uniform_event_row]/nevents_broad)
        ncvar_put_extra(fid, 'event_rate_annual_lower_ci', 
            event_rates_lower[event_uniform_event_row]/nevents_broad)
        nc_close(fid)
    })

}


# Do computations
source_envs = vector(mode='list', length=length(source_names))
names(source_envs) = source_names

source_log_dir = config$sourcezone_log_directory
dir.create(source_log_dir, showWarnings=FALSE)

for(i in 1:length(source_names)){

    source_envs[[i]] = source_rate_environment_fun(
        sourcezone_parameters[i,])

    # Write parameters to a log for later checks
    log_filename = paste0(source_log_dir, '/', names(source_envs)[i], '_', i, '.log')
    capture.output(source_envs[[i]]$sourcepar, log_filename)
}

# Write rates to netcdf
for(i in 1:length(source_names)){
    write_rates_to_event_table(source_envs[[i]])
}


