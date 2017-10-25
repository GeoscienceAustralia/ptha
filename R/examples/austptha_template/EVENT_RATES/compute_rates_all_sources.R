# Compute the rates for all source-zones, and assign nominal 'individual event' rates
# to each tsunami event in the synthetic catalogue. 
#
# Note the 'individual event' rates are not by themselves particularly
# meaningful [they simply involve distributing the rate for each Mw value over
# all events with that Mw]. For example, the nominal 'individual event' rates
# are inversely related to the number of events in the magnitude category, which is
# somewhat arbitrary [especially for stochastic slip]. 
#

library(rptha)
config = new.env()
source('config.R', local=config)

gcmt_access = new.env()
source('gcmt_subsetter.R', local=gcmt_access)
# FIXME: Edit GCMT data selection using strike of source-zone.

#
# INPUTS
# 

sourcezone_parameter_file = config$sourcezone_parameter_file 
sourcezone_parameters = read.csv(sourcezone_parameter_file, stringsAsFactors=FALSE)

# Get the source-name. Segmented cases have an extra 'segment name' that
# distinguishes them
source_segment_names = paste0(sourcezone_parameters$sourcename, 
    sourcezone_parameters$segment_name)

# Never allow Mw_max to be greater than this
MAXIMUM_ALLOWED_MW_MAX = config$MAXIMUM_ALLOWED_MW_MAX  # 9.8
# Never allow Mw_max to be less than this
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

# Currently the code only treats pure normal or pure thrust events
stopifnot(all(sourcezone_parameters$rake[!is.na(sourcezone_parameters$rake)] %in% c(-90, 90)))

# Get environment we can use to provide spatially variable convergence
# information
source('make_spatially_variable_source_zone_convergence_rates.R')
# FIXME: Mappings to some Bird source-zones require manual edits.
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

    # We might be on a specific segment 
    segment_name = sourcezone_parameters_row$segment_name
    source_segment_name = paste0(source_name, segment_name)

    # Find the lower/upper alongstrike numbers for this segment. If missing, we
    # assume the 'segment' is actually the entire source-zone
    alongstrike_lower = as.numeric(
        sourcezone_parameters_row$segment_boundary_alongstrike_index_lower)
    if(is.na(alongstrike_lower)) alongstrike_lower = 1

    alongstrike_upper = as.numeric(
        sourcezone_parameters_row$segment_boundary_alongstrike_index_upper)
    if(is.na(alongstrike_upper)) alongstrike_upper = Inf

    stopifnot(alongstrike_lower < alongstrike_upper)

    # Get the CMT data in this segment
    target_rake = bird2003_env$unit_source_tables[[source_name]]$rake[1]
    stopifnot(target_rake %in% c(-90, 90))
    stopifnot(all(bird2003_env$unit_source_tables[[source_name]]$rake == target_rake))

    # Get GCMT data in this source-zone, if it is 'near-enough' to pure thrust or normal
    # Many default parameters controlling events we select are define in config.R.
    gcmt_data = gcmt_access$get_gcmt_events_in_poly(source_name, 
        alongstrike_index_min=alongstrike_lower,
        alongstrike_index_max=alongstrike_upper,
        target_rake_value=target_rake,
        filter_by_strike = (target_rake != -90), # For normal fault sources, allow the strike to go either way (like outer rise events)
        unit_source_table = bird2003_env$unit_source_tables[[source_name]])


    # Get a vector which is true/false depending on whether each unit-source is
    # inside this particular segment
    is_in_segment = 
        ((bird2003_env$unit_source_tables[[source_name]]$alongstrike_number >= alongstrike_lower) &
         (bird2003_env$unit_source_tables[[source_name]]$alongstrike_number <= alongstrike_upper))

    # Source area 
    unit_source_areas = bird2003_env$unit_source_tables[[source_name]]$length * 
        bird2003_env$unit_source_tables[[source_name]]$width 

    sourcepar$area = sum(unit_source_areas)
    sourcepar$area_in_segment = sum(unit_source_areas*is_in_segment)

    #
    # Gutenberg Richter b-value
    #
    sourcepar$b = sourcezone_parameters_row[1, c('bmin', 'bpref', 'bmax')]
    sourcepar$b = approx(as.numeric(sourcepar$b), n=nbins)$y
    sourcepar$b_p = rep(1, length(sourcepar$b))/length(sourcepar$b)

    #
    # Coupling
    #
    #sourcepar$coupling   = c(0.05, 0.5, 0.95)  # Agnostic approach
    sourcepar$coupling = sourcezone_parameters_row[1, c('cmin', 'cpref', 'cmax')]
    sourcepar$coupling = approx(as.numeric(sourcepar$coupling), n=nbins)$y
    sourcepar$coupling_p = rep(1, length(sourcepar$coupling))/length(sourcepar$coupling)

    #
    # Mw_max
    #

    min_mw_max = max(
        sourcezone_parameters_row$mw_max_observed + mw_observed_perturbation,
        MINIMUM_ALLOWED_MW_MAX)

    max_mw_max_strasser = Mw_2_rupture_size_inverse(sourcepar$area_in_segment, relation = sourcezone_parameters_row$scaling_relation, CI_sd=-1) 

    sourcepar$Mw_max = c(
        ## Largest observed plus a small value,
        min_mw_max,
        # Middle Mw
        0.5*(max_mw_max_strasser + min_mw_max),
        # Another middle Mw
        #Mw_2_rupture_size_inverse(sourcepar$area_in_segment/2, relation = sourcezone_parameters_row$scaling_relation, CI_sd=0),
        # Another middle Mw
        #Mw_2_rupture_size_inverse(sourcepar$area_in_segment, relation = sourcezone_parameters_row$scaling_relation, CI_sd=0),
        # Upper mw [Strasser + 1SD]
        max_mw_max_strasser)
    #
    # Simple test -- all weight on full source rupture area 
    #sourcepar$Mw_max = c(Mw_2_rupture_size_inverse(sourcepar$area_in_segment, relation = sourcezone_parameters_row$scaling_relation, CI_sd=0),
    #    Mw_2_rupture_size_inverse(sourcepar$area_in_segment, relation = sourcezone_parameters_row$scaling_relation, CI_sd=0)+0.01)
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
        # "+-config$rake_deviation" of pure thrust (or normal),
        # then use the raw vector. Otherwise, project it onto the nearest
        # within that range
        #
        # This is consistent with our use of data, which extracts earthquakes
        # with rake being within some deviation of pure thrust (or normal)
        #
        #
        # FIXME: The logic below only works for pure thrust (rake = 90), or pure
        # normal (rake = -90)
        stopifnot(all(bird2003_env$unit_source_tables[[source_name]]$rake %in% c(-90, 90)))

        if(bird2003_env$unit_source_tables[[source_name]]$rake[1] == -90){
            # Normal -- in this case, positive 'vel_div' values contribute to the seismic moment
            div_vec = pmax(0, bird2003_env$unit_source_tables[[source_name]]$bird_vel_div) * is_in_segment
        }else{
            # Thrust -- in this case, negative 'vel_div' values contribute to the seismic moment
            div_vec = pmax(0, -bird2003_env$unit_source_tables[[source_name]]$bird_vel_div) * is_in_segment
        }
        if(max(div_vec) <= 0) stop(paste0('No tectonic moment on source ', source_name, '-- suggests an input bug'))

        # Shorthand divergent (or convergent) and right-lateral velocity 
        # NOTE: We zero velocities that are not on this segment. Also using abs(divergent_velocity),
        # we compute positive tectonic moment for both normal and thrust events. 
        bvrl =  bird2003_env$unit_source_tables[[source_name]]$bird_vel_rl * is_in_segment

        # Limit lateral component of motion that we consider, based on the permitted rake
        # deviation from pure thrust
        deg2rad = pi/180
        allowed_rake_deviation_radians = config$rake_deviation * deg2rad
        rl_vec = sign(bvrl) * pmin(abs(bvrl), div_vec*tan(allowed_rake_deviation_radians))

        # NOTE: If we have segmentation, then this source-zone averaged slip
        # value will be a localised value. 
        sourcepar$slip = weighted.mean(
            # Convergent slip
            x = sqrt(div_vec**2 + rl_vec**2),
            # Weighted by area of unit-sources in the segment
            w = unit_source_areas * is_in_segment)

    }else{

        sourcepar$slip = as.numeric(sourcezone_parameters_row$tectonic_slip)*
            as.numeric(sourcezone_parameters_row$convergent_fraction)

    }

    # Account for non-zero dip, and convert from mm/year to m/year -- only
    # averaging over unit-sources in the segment.
    mean_dip = mean_angle(bird2003_env$unit_source_tables[[source_name]]$dip[which(is_in_segment)])
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
        # Make a function which uses Bird's spatially variable convergence to
        # assign conditional probabilities
        conditional_probability_model = 
            bird2003_env$make_conditional_probability_function_uniform_slip(
                source_name, is_in_segment)

    }else{
        ## conditional_probability_model = 'inverse_slip'

        # Make a conditional probability model like 'inverse_slip', but which
        # considers whether the event is in the segment
        conditional_probability_model<-function(event_table_with_fixed_Mw){

            slip_inv = 1/event_table_with_fixed_Mw$slip

            # Find the fraction of the unit sources in each event which are 
            # inside our segment
            fraction_in_segment = slip_inv * 0.0
            for(i in 1:nrow(event_table_with_fixed_Mw)){
                unit_source_indices = get_unit_source_indices_in_event(
                    event_table_with_fixed_Mw[i,])
                fraction_in_segment[i] = mean(is_in_segment[unit_source_indices])
            }
    
            # Conditional probability will entirely be weighted on events that
            # touch the segment. Events that are only partially contained will
            # be downweighted accordingly
            output = slip_inv * fraction_in_segment
            output = output/sum(output)
            return(output)
        }

    }

    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        event_table, 
        conditional_probability_model = conditional_probability_model)    

    # Get data for rate function update
    if(nrow(gcmt_data) > 0){
        gcmt_data_for_rate_function = list(
            # Magnitude
            Mw = gcmt_data$Mw,
            # Time since the start of observation, in years
            t = (gcmt_data$julianDay1900 - 
                gcmt_access$cmt_start_time_julianDay1900)/gcmt_access$days_in_year
        )
    }else{
        gcmt_data_for_rate_function = list(Mw = NULL, t = NULL)
    }

    MW_MIN = MW_MIN

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
        # Need to pass the local segment area here 
        sourcezone_total_area = sourcepar$area_in_segment, 
        event_table = event_table,
        event_conditional_probabilities = event_conditional_probabilities,
        computational_increment = 0.01,
        Mw_frequency_distribution = Mw_frequency_dists,
        Mw_frequency_distribution_prob = Mw_frequency_dists_p,
        update_logic_tree_weights_with_data=TRUE,
        Mw_count_duration = c(gcmt_access$mw_threshold, nrow(gcmt_data), 
            gcmt_access$cmt_duration_years),
        Mw_obs_data = gcmt_data_for_rate_function,
        account_for_moment_below_mwmin = TRUE
        )


    #
    # Nominal rates on all events
    #

    event_rates = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2) - 
         mw_rate_function(event_table$Mw + dMw/2) )

    # Upper credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    event_rates_upper = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$upper_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$upper_ci_inv_quantile) )
        )

    # Lower credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    event_rates_lower = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$lower_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$lower_ci_inv_quantile) )
        )

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
    source_env$scale_rate = scale_rate 
    source_env$add_rate = add_rate

    # Write to netcdf file. Easiest from within source_env
    with(source_env,{

        # Pull this into the local environment with a new name, to prevent that
        # 'promise already under evaluation' error
        scale_rate_local = scale_rate
        add_rate_local = add_rate

        #
        # Function to do slightly modified 'add to netcdf' computation
        #
        # It can put a 'scaled' version of the rate onto the file, and optionally
        # add to the variable already on the file. 
        #
        # This would allow e.g. weighted treatments of segmented and unsegmented
        # models of the source-zone.
        #
        ncvar_put_extra<-function(fid, varname, event_rate, add_rate=add_rate_local, 
            scale_rate=scale_rate_local){

            if(add_rate){
                # Need to add the current rate to the existing value
                extra_rate = ncvar_get(fid, varname)
            }else{
                # Replace the current rate with the existing value
                extra_rate = event_rate * 0
            }

            output_var = event_rate * scale_rate + extra_rate

            if(any(output_var < 0.0)) stop('negative rate')

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

        # Put rates onto the stochastic and variable uniform slip table nc file
        for(slip_type in c('stochastic', 'variable_uniform')){

            event_table_fileC = paste0(
                '../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_', slip_type, 
                '_slip_earthquake_events_tsunami_',
                source_name, '.nc')

            fid = nc_open(event_table_fileC, readunlim=FALSE, write=TRUE)

            # Index corresponding to uniform slip row
            event_uniform_event_row = ncvar_get(fid, 'event_uniform_event_row')
            # Number of events corresponding to event row
            nevents = table(event_uniform_event_row)
            names_nevents = as.numeric(names(nevents))
            stopifnot(all(names_nevents == 1:length(event_rates)))
            # Make an array giving the number of events matching the
            # uniform_event_row, for every stochastic/variable_uniform event
            nevents_broad = nevents[match(event_uniform_event_row, names_nevents )]

            ncvar_put_extra(fid, 'event_rate_annual', 
                event_rates[event_uniform_event_row]/nevents_broad)
            ncvar_put_extra(fid, 'event_rate_annual_upper_ci', 
                event_rates_upper[event_uniform_event_row]/nevents_broad)
            ncvar_put_extra(fid, 'event_rate_annual_lower_ci', 
                event_rates_lower[event_uniform_event_row]/nevents_broad)
            nc_close(fid)
        }

    })

}


# Do computations
source_envs = vector(mode='list', length=length(source_segment_names))
names(source_envs) = source_segment_names

source_log_dir = config$sourcezone_log_directory
dir.create(source_log_dir, showWarnings=FALSE)

for(i in 1:length(source_segment_names)){

    source_envs[[i]] = source_rate_environment_fun(
        sourcezone_parameters[i,])

    # Write parameters to a log for later checks
    log_filename = paste0(source_log_dir, '/', names(source_envs)[i], 
        '_', i, '.log')
    capture.output(source_envs[[i]]$sourcepar, file=log_filename)
}

#
# Plot all the rate curves to a single pdf
#

xlim = c(7.0, 9.7)
ylim = c(1.0e-06, 10)
pdf('rate_curves_on_source_zones.pdf', width=9, height=7)
for(i in 1:length(source_segment_names)){
     
    # Get all the information
    all_rate_curves = source_envs[[i]]$mw_rate_function(NA, return_all_logic_tree_branches=TRUE)

    mw = all_rate_curves$Mw_seq
    plot(xlim, ylim, log='y', col=0, xlab='Mw', ylab='Exceedance Rate')

    # We will plot quantiles 0, 0.1, 0.2, ... 0.8, 0.9, (1.0-eps) 
    # The 'eps' is used because our quantile evaluation function does not work
    # at the extreme end-point. 
    qntls = seq(0, 1, by=0.1)
    qntls[length(qntls)] = 1 - 1.0e-08 # Must be just < 1
    for(j in 1:length(qntls)){
        curve = as.numeric(source_envs[[i]]$mw_rate_function(mw, quantile=qntls[j]))
        # Remove zero values so we can see the line 'drop' on log axes
        curve = pmax(curve, 1e-100) 
        points(mw, curve, t='l', col='grey')
    }

    points(mw, source_envs[[i]]$mw_rate_function(mw), t='o', col='black', pch=17, cex=0.5)

    # Add empirical Mw-vs-rate for GCMT data
    gcmt_data = source_envs[[i]]$gcmt_data
    if(nrow(gcmt_data) > 0){
        rnk = rank(gcmt_data$Mw)
        N = nrow(gcmt_data)

        # Empirical rate 
        aep = (N+1 - rnk) / gcmt_access$cmt_duration_years

        ordr = order(gcmt_data$Mw)

        points(gcmt_data$Mw[ordr], aep[ordr], col='red', pch=19, t='o')

    }
    

    title(names(source_envs)[i])
    grid(col='orange')
    abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), col='orange', lty='dotted')
}

# Globally integrated rates
# Get all the information
all_rate_curves = source_envs[[1]]$mw_rate_function(NA, return_all_logic_tree_branches=TRUE)
mw = all_rate_curves$Mw_seq
rate_vals = mw*0
gcmt_global = data.frame()
for(i in 1:length(source_segment_names)){
    # Sum the rate value 
    rate_vals = rate_vals + (
        source_envs[[i]]$mw_rate_function(mw) * 
        as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight))
    # Get the GCMT data ONLY for the unsegmented models, so we avoid double-counting
    if(source_envs[[i]]$segment_name == ''){
        gcmt_global = rbind(gcmt_global, source_envs[[i]]$gcmt_data)
    }

}

# Make the globally integrated plot
plot(mw, rate_vals, t='o', log='y', xlab='Mw', ylab='Exceedance Rate (events/year)', 
    main='All sources integrated rate', xlim=xlim, ylim=ylim)
if(nrow(gcmt_global) > 0){
    rnk = rank(gcmt_global$Mw)
    N = nrow(gcmt_global)

    # Empirical rate 
    aep = (N+1 - rnk) / gcmt_access$cmt_duration_years

    ordr = order(gcmt_global$Mw)

    points(gcmt_global$Mw[ordr], aep[ordr], col='red', pch=19, t='o')

}
grid(col='orange')
abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), col='orange', lty='dotted')

dev.off()

#
# OUTPUT RATES TO NETCDF FILES BELOW HERE
#

# Zero rates in netcdf files by setting scale_rate to zero
for(i in 1:length(source_segment_names)){
    write_rates_to_event_table(source_envs[[i]], scale_rate = 0.0, add_rate=FALSE)
}
#
# Now add rates to netcdf files, scaled by the row_weight to allow multiple
# weighted models of the source-zone segmentation
#
for(i in 1:length(source_segment_names)){
    rate_scale = as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight)
    write_rates_to_event_table(source_envs[[i]], scale_rate = rate_scale, add_rate=TRUE)
}


#
# Insert a test here, to confirm that rates are written to all files?
#


#
#
#
save.image('compute_rates_all_sources_session.RData')

## for(i in 1:length(source_envs)){
##     print('')
##     print('##################')
##     print(names(source_envs)[i])  
##     xx = source_envs[[i]]$mw_rate_function(NA, return_all_logic_tree_branches=TRUE)
##     print('    mean coupling')
##     print(weighted.mean(xx$all_par$slip_rate, xx$all_par_prob)/source_envs[[i]]$sourcepar$slip)
##     print('    mean b')
##     print(weighted.mean(xx$all_par$b, xx$all_par_prob))
##     print('    mean Mw_max')
##     print(weighted.mean(xx$all_par$Mw_max, xx$all_par_prob))
##     print('    mean Mw_max -- change')
##     print(weighted.mean(xx$all_par$Mw_max, xx$all_par_prob) - weighted.mean(xx$all_par$Mw_max, xx$all_par_prob_prior))
## 
## }

#
# Various other one-off computations that were required
#
if(FALSE){

    #
    # Moment rate from GCMT on sources with outer-rise
    # 
    thrust_sources = c('sunda', 'kermadectonga', 'puysegur', 'newhebrides', 'timor', 'solomon')
    outerrise_sources = c('outerrisesunda', 'outerrise_kermadectonga', 'outerrise_puysegur', 
        'outerrisenewhebrides', 'outer_rise_timor', 'outerrisesolomon')

    thrust_moment_rate = unlist(lapply(source_envs[thrust_sources], f<-function(x) sum(M0_2_Mw(x$gcmt_data$Mw, inverse=TRUE))))
    normal_moment_rate = unlist(lapply(source_envs[outerrise_sources], f<-function(x) sum(M0_2_Mw(x$gcmt_data$Mw, inverse=TRUE))))

    sum(normal_moment_rate)/sum(thrust_moment_rate)

    #
    # Cross check values
    #
    sources_2 = c('kurilsjapan', 'izumariana', 'mexico', 'southamerica', 'alaskaaleutians')
    thrust_moment_rate2 = rep(NA, length(sources_2))
    normal_moment_rate2 = rep(NA, length(sources_2))
    for(i in 1:length(sources_2)){
        # Thrust events
        gcmt_subset = gcmt_access$get_gcmt_events_in_poly(sources_2[i], target_rake_value=90, 
            filter_by_strike=TRUE, unit_source_table=bird2003_env$unit_source_tables[[sources_2[i]]])
        thrust_moment_rate2[i] = sum(M0_2_Mw(gcmt_subset$Mw, inverse=TRUE))
        # Normal events
        gcmt_subset = gcmt_access$get_gcmt_events_in_poly(sources_2[i], target_rake_value=-90, 
            filter_by_strike=FALSE)
        normal_moment_rate2[i] = sum(M0_2_Mw(gcmt_subset$Mw, inverse=TRUE))
    }

    sum(normal_moment_rate2)/sum(thrust_moment_rate2)

    ## Combined

    sum(c(normal_moment_rate2, normal_moment_rate))/sum(c(thrust_moment_rate2, thrust_moment_rate))


    #
    # Moment rate on sources with outer-rise. This was used to correct the originally estimated
    # moment rates on outer-rise sources (which were too low)
    thrust_sources = c('sunda', 'kermadectonga', 'puysegur', 'newhebrides', 'timor', 'solomon')
    for(i in 1:length(thrust_sources)){
        uss = bird2003_env$unit_source_tables[[thrust_sources[i]]]
        slip_x_area = sum(uss$length*uss$width*1e+06*pmax(0, -uss$bird_vel_div/1000))
        print(c(thrust_sources[i], slip_x_area))
    }
    outerrise_sources = c('outerrisesunda', 'outerrise_kermadectonga', 'outerrise_puysegur', 
        'outerrisenewhebrides', 'outer_rise_timor', 'outerrisesolomon')
    for(i in 1:length(outerrise_sources)){
        uss = bird2003_env$unit_source_tables[[outerrise_sources[i]]]
        slip_x_area = sum(uss$length*uss$width*1e+06*pmax(0, uss$bird_vel_div/1000))
        print(c(outerrise_sources[i], slip_x_area))
    }

}
