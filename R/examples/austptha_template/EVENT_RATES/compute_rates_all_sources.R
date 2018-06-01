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


if(config$edge_correct_event_rates){
    edge_correct = new.env()
    source('back_calculate_convergence.R', local=edge_correct)
}

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
MAXIMUM_ALLOWED_MW_MAX_NORMAL = config$MAXIMUM_ALLOWED_MW_MAX_NORMAL
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
bird2003_env = event_conditional_probability_bird2003_factory(
    return_environment=TRUE)


#' Function to evaluate the rates for a given source-zone. This function returns
#' it's environment, so we have easy access to key variables
#'
#' @param sourcezone_parameters_row row from sourcezone_parameters table
#' @param unsegmented_edge_rate_multiplier Factor by which we increase the
#'   conditional probability of ruptures on the edge of the source-zone. Used to
#'   make the spatial distribution of moment release on the source better
#'   approximate the desired value. If NULL, it is computed, or ignored (see
#'   also config$edge_correct_event_rates)
#' @return the function environment
#'
source_rate_environment_fun<-function(sourcezone_parameters_row, unsegmented_edge_rate_multiplier=NULL){

    #
    # PART 1
    # Read parameters, basic datasets, etc
    #

    # Store key parameters for easy write-out later
    sourcepar = list()

    sourcepar$sourcezone_parameters_row = sourcezone_parameters_row
    source_name = sourcezone_parameters_row$sourcename
    sourcepar$name = source_name

    # We might be on a specific segment 
    segment_name = sourcezone_parameters_row$segment_name
    source_segment_name = paste0(source_name, segment_name)

    # If the input segment_name is blank, then we are on a 'full-zource-zone'. For some
    # operations we need to be careful if this is not the case
    is_a_segment = (sourcezone_parameters_row$segment_name != '')
    # Check this is right    
    if(is_a_segment){
        if(sum(sourcezone_parameters$sourcename == sourcezone_parameters_row$sourcename)<=1){
            msg = paste0('Segment name provided on source with only one entry on ', 
                        'sourcezone_parameters. \n In this case the segment name should be blank,',
                        ' or there should be another segment in the table. Check ', source_segment_name)
            stop(msg)
        }
    }else{
        if(sum(sourcezone_parameters$sourcename == sourcezone_parameters_row$sourcename & 
                (sourcezone_parameters$segment_name == '')) != 1){
            msg = paste0('Source appears more than once with blank segment names. Check', 
                source_segment_name)
            stop(msg)
        }
    }

    # Find the lower/upper alongstrike numbers for this segment. If missing, we
    # assume the 'segment' is actually the entire source-zone
    alongstrike_lower = as.numeric(
        sourcezone_parameters_row$segment_boundary_alongstrike_index_lower)
    if(is.na(alongstrike_lower)){
        alongstrike_lower = 1
        stopifnot(!is_a_segment)
    }

    alongstrike_upper = as.numeric(
        sourcezone_parameters_row$segment_boundary_alongstrike_index_upper)
    if(is.na(alongstrike_upper)){
        alongstrike_upper = Inf
        stopifnot(!is_a_segment)
    }

    stopifnot(alongstrike_lower < alongstrike_upper)

    # Get the CMT data in this segment
    target_rake = bird2003_env$unit_source_tables[[source_name]]$rake[1]
    stopifnot(target_rake %in% c(-90, 90))
    stopifnot(all(bird2003_env$unit_source_tables[[source_name]]$rake == target_rake))

    # Get GCMT data in this source-zone, if it is 'near-enough' to pure thrust or normal
    # Many default parameters controlling events we select are defined in config.R.
    gcmt_data = gcmt_access$get_gcmt_events_in_poly(source_name, 
        alongstrike_index_min=alongstrike_lower,
        alongstrike_index_max=alongstrike_upper,
        target_rake_value=target_rake,
        filter_by_strike = (target_rake != -90), # For normal fault sources, allow the strike to go either way (like outer rise events)
        unit_source_table = bird2003_env$unit_source_tables[[source_name]])


    # Get a vector which is true/false depending on whether each unit-source is
    # inside this particular segment. (Not to be confused with 'is_a_segment', which is a scaler logical,
    # telling us if the current source is just a segment of a source-zone)
    is_in_segment = 
        ((bird2003_env$unit_source_tables[[source_name]]$alongstrike_number >= alongstrike_lower) &
         (bird2003_env$unit_source_tables[[source_name]]$alongstrike_number <= alongstrike_upper))
    which_is_in_segment = which(is_in_segment)

    # Source area 
    unit_source_areas = bird2003_env$unit_source_tables[[source_name]]$length * 
        bird2003_env$unit_source_tables[[source_name]]$width 

    sourcepar$area = sum(unit_source_areas)
    sourcepar$area_in_segment = sum(unit_source_areas*is_in_segment)

    #
    # Gutenberg Richter b-value
    #
    sourcepar$b = sourcezone_parameters_row[1, c('bmin', 'bpref', 'bmax')]
    sourcepar$b = approx(as.numeric(sourcepar$b), n=nbins*config$b_subsampling_increase)$y
    sourcepar$b_p = rep(1, length(sourcepar$b))/length(sourcepar$b)

    #
    # Coupling
    #
    if(config$use_uniform_coupling_prior){
        sourcepar$coupling  = config$uniform_coupling_prior 
    }else{
        sourcepar$coupling = sourcezone_parameters_row[1, c('cmin', 'cpref', 'cmax')]
    }

    if(any(sourcepar$coupling == 0)) stop('Cannot treat coupling of zero, instead set prob_Mmax_below_Mmin > 0')
    # Log spacing to equally resolve all coupling values in a relative sense
    # However, we assign probabilities consistent with a uniform distribution
    # (the log spacing is to improve numerical aspects of the discretization,
    # we are not claiming that log(coupling) is uniformly distributed)
    sourcepar$coupling = exp(approx(log(as.numeric(sourcepar$coupling)), 
        n=nbins*config$coupling_subsampling_increase)$y)
    # Discretize the probabilities as a uniform distribution (NOT log-uniform!)
    bin_size = diff(c(min(sourcepar$coupling), sourcepar$coupling, max(sourcepar$coupling)), lag=2)/2
    sourcepar$coupling_p = bin_size/sum(bin_size)

    # We may have put some weight on 'aseismic' behaviour in the magnitude range of interest
    sourcepar$prob_Mmax_below_Mmin = as.numeric(sourcezone_parameters_row$prob_Mmax_below_Mmin)
    if(sourcepar$prob_Mmax_below_Mmin > 0){
        # We have a certain probability of '0' coupling.
        # Append that as an additional bin
        p0 = sourcepar$prob_Mmax_below_Mmin
        stopifnot(p0 > 0 & p0 <= 1)
        sourcepar$coupling = c(0, sourcepar$coupling)
        sourcepar$coupling_p = c(p0, sourcepar$coupling_p * (1-p0))
    }

    #
    # Mw_max
    #

    min_mw_max = max(
        sourcezone_parameters_row$mw_max_observed + mw_observed_perturbation,
        MINIMUM_ALLOWED_MW_MAX)

    # Upper bound mw-max --> scaling relation, with area at -1 standard deviation
    max_mw_max_strasser_area = Mw_2_rupture_size_inverse(sourcepar$area_in_segment, 
        relation = sourcezone_parameters_row$scaling_relation, CI_sd=-1) 
    
    # Alternative upper bound mw-max --> scaling relation, with width at -2 standard deviation
    # This is a good idea from a practical perspective, because I generate variable-uniform
    # and stochastic slip events by simulating width/length within +-2 SD limits. If width
    # cannot accomodate this, then it is truncated (and in 50% of cases, length
    # is expanded commensurately, while in the other 50% of cases, length is
    # unchanged). If I allow Mw_max such that the source-zone width is < 2 standard deviations
    # below the scaling value, then I will end up having very small area earthquakes in the
    # cases where width is truncated, but length is not expanded. On some
    # source-zones, this can lead to 100s of metres of mean slip even for e.g.
    # Mw 9.1 earthquakes (say on Manus which is only 1-unit-source down-dip, but long enough
    # to have an area-based Mw_max of 9.15).
    sourcezone_widths = aggregate(
        bird2003_env$unit_source_tables[[source_name]]$width[which_is_in_segment],
        by=list(bird2003_env$unit_source_tables[[source_name]]$alongstrike_number[which_is_in_segment]),
        sum)
    mean_sourcezone_width = mean(sourcezone_widths$x)

    max_mw_max_strasser_width = uniroot(
        f<-function(x){ 
            output = Mw_2_rupture_size(x, relation=sourcezone_parameters_row$scaling_relation,
                detailed=TRUE, CI_sd=2)$minus_CI['width'] - mean_sourcezone_width
            return(output)
            },
        interval=c(2, 20), # Mw must be between these bounds!
        tol=1e-10
        )$root

    max_mw_max_strasser = min(max_mw_max_strasser_area, max_mw_max_strasser_width)

    if(max_mw_max_strasser > min_mw_max){

        sourcepar$Mw_max = c(
            # Largest observed plus a small value,
            min_mw_max,
            # Upper mw [Strasser - 1SD area, OR Strasser -2SD width]
            max_mw_max_strasser)
    }else{

        stop(paste0('Scaling relation mw_max < max observed @', source_segment_name))
    }

    # Ensure ordered
    sourcepar$Mw_max = sort(sourcepar$Mw_max)

    # Check it is correctly ordered (of course!)
    stopifnot(all(diff(sourcepar$Mw_max) > 0))
    # Ensure all Mw meet out constraints
    sourcepar$Mw_max = pmax(sourcepar$Mw_max, min_mw_max)

    # Interpolate
    sourcepar$Mw_max = approx(sourcepar$Mw_max, n=nbins*config$mwmax_subsampling_increase)$y
    # Assign equal probabilities to all
    sourcepar$Mw_max_p = rep(1, length(sourcepar$Mw_max) )/length(sourcepar$Mw_max)
   
    # Clip AFTER interpolation 
    if(target_rake == -90){
        sourcepar$Mw_max = pmin(sourcepar$Mw_max, MAXIMUM_ALLOWED_MW_MAX_NORMAL)
    }else{
        sourcepar$Mw_max = pmin(sourcepar$Mw_max, MAXIMUM_ALLOWED_MW_MAX)
    }


    #
    # PART 2.
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
        # Using a constant convergence rate everywhere

        sourcepar$slip = as.numeric(sourcezone_parameters_row$tectonic_slip)*
            as.numeric(sourcezone_parameters_row$convergent_fraction)

        div_vec = rep(sourcepar$slip, length(unit_source_areas)) * is_in_segment

    }

    # Account for non-zero dip, and convert from mm/year to m/year -- only
    # averaging over unit-sources in the segment.
    mean_dip = mean_angle(bird2003_env$unit_source_tables[[source_name]]$dip[which(is_in_segment == 1)])
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

    if(target_rake == 90){
        #
        # Cumulative distribution function of the difference between the 'variable shear modulus'
        # magnitude, and the 'fixed shear modulus' magnitude. By treating these differences as
        # 'errors in observed magnitudes', we get a convenient method of treating variable shear 
        # modulus.
        # This will differ depending on whether we compute the function using
        # uniform/variable-uniform/heterogeneous slip scenarios. However, experimentation
        # suggests the variation is not great, although uniform has clearer 'discretization'
        # artefacts due to reduced variability. So we use stochastic
        #

        unit_source_depth = bird2003_env$unit_source_tables[[source_name]]$depth
        unit_source_mu_variable = shear_modulus_depth(unit_source_depth, type='default')
    
        stochastic_slip_event_table_file = paste0('../SOURCE_ZONES/', source_name, 
        '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_', source_name, '.nc')

        # Get slip and indices
        stochastic_slip_fid = nc_open(stochastic_slip_event_table_file, readunlim=FALSE)
        stoc_event_slip_string = ncvar_get(stochastic_slip_fid, 'event_slip_string')
        stoc_event_index_string = ncvar_get(stochastic_slip_fid, 'event_index_string')
        stoc_mw_mu_constant = ncvar_get(stochastic_slip_fid, 'Mw')
        stoc_mw_mu_constant = round(stoc_mw_mu_constant, 3) # Address imperfect floating point storage in netcdf
        nc_close(stochastic_slip_fid); rm(stochastic_slip_fid)

        stoc_ess = lapply(stoc_event_slip_string, f<-function(x) as.numeric(strsplit(x, '_')[[1]]))
        stoc_eis = lapply(stoc_event_index_string, f<-function(x) as.numeric(strsplit(x, '-')[[1]]))

        stoc_mw_mu_variable = rep(NA, length(stoc_ess)) 
        for(ei in 1:length(stoc_mw_mu_variable)){
            inds = stoc_eis[[ei]]
            slp = stoc_ess[[ei]]
            stoc_moment = sum(unit_source_areas[inds] * 1e+06 * unit_source_mu_variable[inds] * slp)
            stoc_mw_mu_variable[ei] = M0_2_Mw(stoc_moment)
        }
        # Save memory
        rm(stoc_ess, stoc_eis, stoc_event_slip_string, stoc_event_index_string)

        # Difference between variable slip and uniform slip shear modulus
        mw_obs_deviation = stoc_mw_mu_variable - stoc_mw_mu_constant
        # Sanity check (specific to our case)
        stopifnot( (min(mw_obs_deviation) > -0.32) & (max(mw_obs_deviation) < 0.25) )
   
        # Conditional empirical CDF of difference between 'Mw observation' and
        # 'Mw with constant shear modulus' 
        mw_deviation_cdf_variable_shear_modulus = make_conditional_ecdf(mw_obs_deviation, stoc_mw_mu_constant)

        # Reduce memory
        rm(mw_obs_deviation, stoc_mw_mu_constant); gc()

    }else{
        # For normal fault events, ignore shear modulus variations
        mw_deviation_cdf_variable_shear_modulus = NULL
    }
    

    #
    # PART 3
    # Get a preliminary event conditional probability model; use it to compute
    # the mw_rate_function, then update the event conditional probability model
    # (this will not affect the mw_rate_function)
    #

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

    # NOTE: These conditional probabilities may be updated to deal with 'edge-effects'
    # further in the code. The aim is to make the integrated slip more consistent with
    # moment conservation -- whereas the current approach tends to concentrate moment
    # release more towards the centre of the rupture.
    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        event_table, 
        conditional_probability_model = conditional_probability_model)    

    # Get data for rate function update
    if(nrow(gcmt_data) > 0){
        gcmt_data_for_rate_function = list(
            # Magnitude
            Mw = gcmt_data$Mw, 
            # Time since the start of observation, in years
            #t = (gcmt_data$julianDay1900 - 
            #    gcmt_access$cmt_start_time_julianDay1900)/gcmt_access$days_in_year
            t = NULL # Censored likelihood is biased for rates, better to use poisson count approach.
        )

        if(min(sourcepar$coupling) == 0) stop('Cannot have zero coupling logic tree branch when GCMT data is present')
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
        computational_increment = 0.02,
        Mw_frequency_distribution = Mw_frequency_dists,
        Mw_frequency_distribution_prob = Mw_frequency_dists_p,
        update_logic_tree_weights_with_data=config$update_logic_tree_weights_with_data,
        Mw_count_duration = c(gcmt_access$mw_threshold, nrow(gcmt_data), 
            gcmt_access$cmt_duration_years),
        Mw_obs_data = gcmt_data_for_rate_function,
        account_for_moment_below_mwmin = TRUE,
        mw_observation_error_cdf = mw_deviation_cdf_variable_shear_modulus
        )

    # Nearly finished .....
    #
    # Optionally adjust the conditional probability of the events, to better satisfy spatial
    # variations in seismic moment conservation
    #
    if(config$edge_correct_event_rates){
        #
        # Derive new event conditional probabilities which better represent
        # spatial variations in tectonic slip. We do this by inflating the
        # conditional probabilities of 'edge_events' (i.e. events which
        # occur on the boundary of the unsegmented source-zone) by a factor 'edge_multiplier'.
        #
        # Here, we find an 'edge_multiplier' that well represents spatial
        # variations in convergence.
        #
        fun_to_optimize<-function(edge_multiplier, return_conditional_probabilities=FALSE){

            # Rates based on the 'edge_multiplier = 0' case
            event_rates = event_conditional_probabilities * 
                (mw_rate_function(event_table$Mw - dMw/2) - 
                 mw_rate_function(event_table$Mw + dMw/2) )

            # Compute spatially variable integrated slip, given the edge_multiplier argument
            back_calculate_convergence_env = new.env()
            back_calculate_convergence_env = edge_correct$back_calculate_convergence(
                sourcename=source_name,
                slip_type='uniform',
                edge_multiplier=edge_multiplier,
                uss = bird2003_env$unit_source_tables[[source_name]],
                event_rates = event_rates,
                event_Mw = event_table$Mw,
                event_index_string = event_table$event_index_string,
                slip = event_table$slip)
            model = back_calculate_convergence_env$output$integrated_slip

            #
            # We want the shape of the convergent slip to be like 'div_vec'.
            # Coupling might rescale this. Return a goodness of fit value that
            # reflects that, unless we explicitly request to return the above
            # function environment
            #
            if(!return_conditional_probabilities){ 
                output = sum( is_in_segment*( model/sum(model) - div_vec/sum(div_vec) )**2 )
            }else{
                # This is useful once we know the best edge_multiplier
                output = back_calculate_convergence_env
            }
            return(output)
        }


        if(!is_a_segment){

            # Tests to weed out problematic cases
            f1 = fun_to_optimize(0.0)
            f2 = fun_to_optimize(10.0)

            if(isTRUE(all.equal(f1, f2))){
                # This should only happen if events on the boundary all have
                # a-prior weight of zero, which should only happen for 'is_a_segment'
                # cases
                msg = paste0('Use of an edge multiplier is having no impact on source ', 
                    source_segment_name)
                stop(msg)
            }else{
                # Find the optimal edge_multiplier
                best_edge_mult = optimize(fun_to_optimize, lower=0, upper=30)
                # Check that results are not hitting these boundaries
                if(best_edge_mult$minimum < 1.0e-03 | best_edge_mult$minimum > 29.999){
                    print(paste0('Check edge_multiplier on ', source_segment_name))
                }
                # Get the conditional probabilities from the 'best' edge_multiplier
                best_conv_env = fun_to_optimize(best_edge_mult$minimum, return_conditional_probabilities=TRUE)
                event_conditional_probabilities = best_conv_env$new_conditional_probability
                rm(best_conv_env)
                gc()
            }

        }else{
            #
            # We are on a segment. It might be hard to stably estimate the edge_multiplier, because
            # the edge events might have little impact on the convergence on this source. Therefore,
            # it seems better to use the edge_multiplier from the unsegmented source
            if( is.null(unsegmented_edge_rate_multiplier) ){
                msg = 'Must provide an unsegmented_edge_rate_multiplier on segments'
                stop(msg)
            }
            best_edge_mult = list()
            best_edge_mult$minimum = unsegmented_edge_rate_multiplier
            if(best_edge_mult$minimum < 1.0e-03 | best_edge_mult$minimum > 29.999){
                print(paste0('Check edge_multiplier on ', source_segment_name))
            }
            # Get the conditional probabilities from the 'best' edge_multiplier
            best_conv_env = fun_to_optimize(best_edge_mult$minimum, return_conditional_probabilities=TRUE)
            event_conditional_probabilities = best_conv_env$new_conditional_probability
            rm(best_conv_env)
            gc()
        }
        sourcepar$best_edge_multiplier = best_edge_mult
    }



    #
    # Nominal rates on all events
    #

    event_rates = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2) - 
         mw_rate_function(event_table$Mw + dMw/2) )

    event_rates_mu_vary = event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, account_for_mw_obs_error=TRUE) )

    # Upper credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    # NOTE: Later we sum over [row-weights x event_rates_upper] on segmented
    # (or optionally segmented) source-zones. This is valid if we assume
    # that ALL representations of the source behave like separate, real sources, with
    # the given fraction of the moment rate, AND furthermore that epistemic uncertainties
    # in all the sources are co-monotonic [which sounds a sensible assumption, because it
    # prevents us 'averaging away' the risk of rare events].
    # It is clearly reasonable for separate segments to be like separate, real sources. It is
    # less obvious that the 'full-source-zone' model should be as well, but we can think of 
    # this as us representing the source as having 'some tendency for rupture-segment-sized events,
    # but also some tendency to behave as a full source-zone.
    event_rates_upper = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$upper_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$upper_ci_inv_quantile) )
        )

    event_rates_upper_mu_vary = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$upper_ci_inv_quantile, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$upper_ci_inv_quantile, account_for_mw_obs_error=TRUE) )
        )

    # Lower credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    # NOTE: Later we sum over [row-weights x event_rates_upper] on segmented
    # (or optionally segmented) source-zones. This is valid if we assume
    # that ALL representations of the source behave like separate, real sources, with
    # the given fraction of the moment rate, AND furthermore that epistemic uncertainties
    # in all the sources are co-monotonic [which sounds a sensible assumption, because it
    # prevents us 'averaging away' the risk of rare events]
    # It is clearly reasonable for separate segments to be like separate, real sources. It is
    # less obvious that the 'full-source-zone' model should be as well, but we can think of 
    # this as us representing the source as having 'some tendency for rupture-segment-sized events,
    # but also some tendency to behave as a full source-zone.
    event_rates_lower = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$lower_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$lower_ci_inv_quantile) )
        )
    event_rates_lower_mu_vary = as.numeric(
        event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, 
            quantiles=config$lower_ci_inv_quantile, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, 
            quantiles=config$lower_ci_inv_quantile, account_for_mw_obs_error=TRUE) )
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

            ncvar_put(fid, varname, output_var)

            return(invisible())
        }

        # Put rates onto the event table nc file
        fid = nc_open(event_table_file, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'rate_annual', event_rates)
        # Summation of credible intervals OK for co-monotonic epistemic uncertainties
        ncvar_put_extra(fid, 'rate_annual_upper_ci', event_rates_upper)
        ncvar_put_extra(fid, 'rate_annual_lower_ci', event_rates_lower)
        nc_close(fid)

        # Put rates onto the uniform slip table nc file
        event_table_fileB = paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc')
        fid = nc_open(event_table_fileB, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'event_rate_annual', event_rates)
        # Summation of credible intervals OK for co-monotonic epistemic uncertainties
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
            # Summation of credible intervals OK for co-monotonic epistemic uncertainties
            ncvar_put_extra(fid, 'event_rate_annual_upper_ci', 
                event_rates_upper[event_uniform_event_row]/nevents_broad)
            ncvar_put_extra(fid, 'event_rate_annual_lower_ci', 
                event_rates_lower[event_uniform_event_row]/nevents_broad)
            nc_close(fid)

            #
            # Now to the same, for the file that contains only the earthquakes
            # 
            event_table_fileD = paste0(
                '../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_', slip_type, 
                '_slip_earthquake_events_',
                source_name, '.nc')

            fid = nc_open(event_table_fileD, readunlim=FALSE, write=TRUE)

            # Index corresponding to uniform slip row
            event_uniform_event_row = ncvar_get(fid, 'uniform_event_row')
            # Number of events corresponding to event row
            nevents = table(event_uniform_event_row)
            names_nevents = as.numeric(names(nevents))
            stopifnot(all(names_nevents == 1:length(event_rates)))
            # Make an array giving the number of events matching the
            # uniform_event_row, for every stochastic/variable_uniform event
            nevents_broad = nevents[match(event_uniform_event_row, names_nevents )]

            ncvar_put_extra(fid, 'rate_annual', 
                event_rates[event_uniform_event_row]/nevents_broad)
            # Summation of credible intervals OK for co-monotonic epistemic uncertainties
            ncvar_put_extra(fid, 'rate_annual_upper_ci', 
                event_rates_upper[event_uniform_event_row]/nevents_broad)
            ncvar_put_extra(fid, 'rate_annual_lower_ci', 
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

# Function to run in parallel over all source zones
parfun<-function(i, unsegmented_edge_rate_multiplier=NULL){

    output = try( source_rate_environment_fun(sourcezone_parameters[i,],
        unsegmented_edge_rate_multiplier = unsegmented_edge_rate_multiplier))

    if(!is.environment(output)){
        # Return an environment if it fails. This should make debugging easier
        failed_environment = as.environment(list(
            source_rate_environment_fun_return_value = as.list(output),
            i = i,
            sourcezone_parameters_i = sourcezone_parameters[i,],
            unsegmented_edge_rate_multiplier = unsegmented_edge_rate_multiplier))

        output = failed_environment
    }else{
        # It seems to have worked.
        # Write parameters to a log for later checks
        log_filename = paste0(source_log_dir, '/', source_segment_names[i], 
            '_', i, '.log')
        capture.output(output$sourcepar, file=log_filename)
    }
    return(output)
}

#
# Run for all source zones
#
# Do all unsegmented cases first. This is done, because we will use edge_multipliers
# from unsegmented segments to treat segmented cases (because stably and consistently 
# estimating edge_multipliers may be tricky with segmented cases, given that the multipliers
# may have little impact on some segments)
unseg = which(sourcezone_parameters$segment_name == '')
if(config$MC_CORES > 1){
    # Parallel run
    library(parallel)
    if(length(unseg) > 0){
        source_envs[unseg] = mclapply(as.list(1:length(source_segment_names))[unseg], parfun, 
            mc.cores=config$MC_CORES)
    }
}else{
    # Serial run
    if(length(unseg) > 0){
        source_envs[unseg] = lapply(as.list(1:length(source_segment_names))[unseg],
            parfun)
    }
}

# Get the edge_multiplier for the segmented models, and copy to the unsegmented models
unsegmented_edge_rate_multiplier = vector(mode='list', length=length(source_segment_names))
if(config$edge_correct_event_rates){
    for(i in 1:length(source_segment_names)){
        if(i %in% unseg){
            k = which(sourcezone_parameters$sourcename == sourcezone_parameters$sourcename[i])
            # Copy the edge_multiplier from the unsegmented model to all other models on the same source
            for(j in k){
                unsegmented_edge_rate_multiplier[[j]] = source_envs[[i]]$best_edge_mult$minimum
            }
        }
    }
}

#
# Finally run segmented models
#
seg = which(sourcezone_parameters$segment_name != '')
if(config$MC_CORES > 1){
    # Parallel run
    library(parallel)

    if(length(seg) > 0){
        source_envs[seg] = mcmapply(parfun, 
            i = as.list(1:length(source_segment_names))[seg], 
            unsegmented_edge_rate_multiplier=unsegmented_edge_rate_multiplier[seg], 
            SIMPLIFY=FALSE,
            mc.cores=config$MC_CORES)
    }

}else{
    # Serial run
    if(length(seg) > 0){
        source_envs[seg] = mapply(parfun, 
            i=as.list(1:length(source_segment_names))[seg],
            unsegmented_edge_rate_multiplier=unsegmented_edge_rate_multiplier[seg], 
            SIMPLIFY=FALSE)
    }
}

#
# Quick check for order mistakes
#
for(i in 1:length(source_envs)){
    if(names(source_envs)[i] != source_envs[[i]]$source_segment_name){
        print(paste0('Problem with ', names(source_envs)[i]))
    }
}


#
# OUTPUT RATES TO NETCDF FILES BELOW HERE
#
if(config$write_to_netcdf){

    # Zero rates in netcdf files by setting scale_rate to zero
    for(i in 1:length(source_segment_names)){
        write_rates_to_event_table(source_envs[[i]], scale_rate = 0.0, add_rate=FALSE)
    }
    #
    # Now add rates to netcdf files, scaled by the row_weight to allow multiple
    # weighted models of the source-zone segmentation. Note we are adding epistemic
    # credible intervals of rates here too. This is ok if the epistemic uncertainties
    # are co-monotonic [i.e. 95% quantile of segmentA ==> 95% quantile of segmentB
    # and segment C and the full source-zone]. That avoids risk-reduction through averaging,
    # which seems a reasonable approach [given it's hard to specify 'correlations' between
    # our epistemic uncertainties on different sources, but we could well imagine the 
    # correlations not being zero]
    #
    for(i in 1:length(source_segment_names)){
        rate_scale = as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight)
        write_rates_to_event_table(source_envs[[i]], scale_rate = rate_scale, add_rate=TRUE)
    }
}

#
# Plot all the rate curves to a single pdf
#

xlim = c(7.0, 9.7)
ylim = c(1.0e-06, 10)
pdf('rate_curves_on_source_zones.pdf', width=9, height=7)
for(i in 1:length(source_segment_names)){
     
    # Get all the information
    all_rate_curves = source_envs[[i]]$mw_rate_function(NA, 
        return_all_logic_tree_branches=TRUE)

    mw = all_rate_curves$Mw_seq
    par(mfrow=c(1,1))
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

    points(mw, source_envs[[i]]$mw_rate_function(mw), t='o', col='black', 
        pch=17, cex=0.5)

    # Plot the by-earthquake-event rates, integrated. Should be 'the same' (up
    # to mw discretization) as the fitted curve
    empirical_mean_curve = sapply(mw, 
        f<-function(x) sum(source_envs[[i]]$event_rates * (source_envs[[i]]$event_table$Mw >= x)))
    points(mw, empirical_mean_curve, pch=19, cex=0.2, col='green')

    # Also read rates that were just written to file -- this will be for the whole source-zone
    fid = nc_open(source_envs[[i]]$event_table_file, readunlim=FALSE)
    event_rates_file = ncvar_get(fid, 'rate_annual')
    event_Mw_file = round(ncvar_get(fid, 'Mw'), 3)
    nc_close(fid)
    empirical_mean_curve = sapply(mw, f<-function(x) sum(event_rates_file * (event_Mw_file >= x)))
    points(mw, empirical_mean_curve, pch=19, cex=0.2, col='pink')

    # Mean prior curve
    #mean_prior_curve = colMeans(all_rate_curves$all_rate_matrix)
    mean_prior_curve = apply(all_rate_curves$all_rate_matrix, 2, 
        f<-function(x) weighted.mean(x, w=all_rate_curves$all_par_prob_prior))
    points(mw, mean_prior_curve, t='l', lwd=2, col='orange', lty='dashed')

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

    # Posterior of slip, Mw_max, and b
    num_diff<-function(x, y){
        # Quick numerical derivative
        N = length(x)
        c( (y[2]-y[1])/(x[2]-x[1]), 
           (y[3:N] - y[1:(N-2)])/(x[3:N] - x[1:(N-2)]), 
           (y[N] - y[N-1])/(x[N] - x[N-1]) )
    }

    plot_derivs <-function(var, prior=FALSE){
        vars = sort(unique(all_rate_curves$all_par[[var]]))
        vars_cdf = sapply(vars, f<-function(x){
            sum(all_rate_curves$all_par_prob*(all_rate_curves$all_par[[var]] <= x))} )
        vars_cdf_prior = sapply(vars, f<-function(x){
            sum(all_rate_curves$all_par_prob_prior*(all_rate_curves$all_par[[var]] <= x))})
        vars_dens = num_diff(vars, vars_cdf)
        vars_dens_prior = num_diff(vars, vars_cdf_prior)
        var_info=list(x=vars, density=vars_dens, prior_density=vars_dens_prior)

        plot(vars, vars_dens, t='h', lend=1, ylim=c(0, max(vars_dens)), 
            main=paste0(var, ' density'), xlab=var, ylab='density')
        points(vars, vars_dens_prior, t='l', col='red', lty='dotted')

        return(invisible(var_info))
    }

    par(mfrow=c(2,2))
    
    plot_derivs('slip_rate')
    plot_derivs('Mw_max')
    plot_derivs('b')
    # Dummy plot to fill the page
    plot(0:1, col='white')
}

# Globally integrated rates
mw = seq(MW_MIN, MAXIMUM_ALLOWED_MW_MAX, by=dMw/2) #all_rate_curves$Mw_seq
rate_vals = mw*0
gcmt_global = data.frame()
for(i in 1:length(source_segment_names)){
    # Sum the rate value 
    rate_vals = rate_vals + (
        source_envs[[i]]$mw_rate_function(mw) * 
        as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight))
    # Get the GCMT data ONLY for the unsegmented models, so we avoid
    # double-counting
    if(source_envs[[i]]$segment_name == ''){
        gcmt_global = rbind(gcmt_global, source_envs[[i]]$gcmt_data)
    }

}

# Make the globally integrated plot
par(mfrow=c(1,1))
plot(mw, rate_vals, t='o', log='y', xlab='Mw', ylab='Exceedance Rate (events/year)', 
    main='All sources integrated rate', xlim=xlim, ylim=ylim)
if(nrow(gcmt_global) > 0){
    rnk = rank(gcmt_global$Mw, ties='random')
    N = nrow(gcmt_global)

    # Empirical rate 
    aep = (N+1 - rnk) / gcmt_access$cmt_duration_years

    ordr = order(gcmt_global$Mw)

    points(gcmt_global$Mw[ordr], aep[ordr], col='red', pch=19, t='o')

    # Add poisson type confidence intervals
    upper_ci = ordr*0
    lower_ci = ordr*0
    for(ii in ordr){
        p1 = poisson.test(N+1-rnk[ii], T=gcmt_access$cmt_duration_years, conf.level=0.95)
        upper_ci[ii] = p1$conf.int[2]
        lower_ci[ii] = p1$conf.int[1]
    }
    points(gcmt_global$Mw[ordr], lower_ci[ordr], t='l', col='red', lty='dashed', lwd=2)
    points(gcmt_global$Mw[ordr], upper_ci[ordr], t='l', col='red', lty='dashed', lwd=2)

}
grid(col='orange')
abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), col='orange', lty='dotted')

dev.off()

## #
## # Write to netcdf with a sourcezone specific approach
## #
## sourcenames = sourcezone_parameters$sourcename
## unique_sourcenames = unique(sourcenames)
## for(i in 1:length(unique_sourcenames)){
##     # When ruptures include segmentation, we need to integrate multiple source_envs
##     k = which(sourcenames == unique_sourcenames[[i]])
## 
##     # Extract key data for all segments
##     segment_altb = list() # all logic tree branches
##     segment_event_cond_prob = list()
##     segment_event_Mw = list()
##     segment_row_weight = list()
##     for(j in 1:length(k)){
##         kj = k[j]
##         segment_altb[[j]] = source_envs[[kj]]$mw_rate_function(NA, 
##             return_all_logic_tree_branches=TRUE)
##         segment_event_cond_prob[[j]] = source_envs[[kj]]$event_conditional_probabilities
##         segment_event_Mw[[j]] = source_envs[[kj]]$event_table$Mw
##         segment_row_weight[[j]] = source_envs[[kj]]$sourcezone_parameters_row$row_weight
## 
##         # Mw values should be consistent
##         stopifnot( all(segment_event_Mw[[j]] == segment_event_Mw[[1]]) )
##         stopifnot( length(segment_event_cond_prob[[j]]) == length(segment_event_cond_prob[[1]]) )
##     }
## 
##     mw_events = sort(unique(segment_event_Mw[[1]]))
##     mw_bin_lower = mw_events - config$dMw/2
##     mw_bin_upper = mw_events + config$dMw/2
## 
##     source_event_rates_mean = segment_event_cond_prob[[1]]*0
##     source_event_rates_upper_CI = segment_event_cond_prob[[1]]*0
##     source_event_rates_lower_CI = segment_event_cond_prob[[1]]*0
## 
##     # Loop over Mw bins defining event magnitudes
##     for(m in 1:length(mw_events)){
##         # Indices of events with Mw == mw_events[j] (identical for all segments, as checked above)
##         evnts = which(segment_event_Mw[[1]] == mw_events[m])
## 
##         # Loop over every segment and get rate-vs-probability curve for the m'th mw_bin
##         rate_vs_prob_sorted = list()
##         for(j in 1:length(segment_altb)){
##             # For every logic tree branch, get the rate of events in [
##             # mw_bin_lower[m], mw_bin_upper[m] ]
##             logic_tree_rate_in_bin = apply(segment_altb[[j]]$all_rate_matrix, 1, 
##                 f<-function(x){
##                      diff(approx(segment_altb[[j]]$Mw_seq, x, xout=c(mw_bin_upper[m], mw_bin_lower[m]))$y)
##                      }
##             )
## 
##             # Probability of each logic tree branch (assuming the segment is 'true')
##             logic_tree_prob = segment_altb[[j]]$all_par_prob
## 
##             rate_order = order(logic_tree_rate_in_bin)
##             # Make a rate-vs-prob curve
##             r1 = logic_tree_rate_in_bin[rate_order]
##             p1 = logic_tree_prob[rate_order]
##             p1_cumsum = cumsum(p1)
##             rate_vs_prob_sorted[[j]] = list(r1=r1, p1=p1, p1_cumsum=p1_cumsum)
##         }
## 
##         for(e in evnts){
##             # Conditional probability of the event on this segment (given
##             # that an event of the same Mw occurred)
##
##             # Idea: Compute a bound on the upper and lower probability quantiles for the summed rates,
##             # given we don't know the dependence structure
##             robust_quantile<-function(alphas, in_segment = c(2,3), e = 54, alpha=0.95, type='upper'){
##                
##                 la = length(in_segment) #end-start+1
##
##                 # Compute the last alpha value from the provided alpha values,
##                 # to ensure sum( 1 - all_alphas) = (1 - alpha)
##
##                 stopifnot(length(alphas) == (la-1))
##
##                 scaler = 1
##
##                 all_alphas = c(alphas, 0)
##                 if(type=='upper'){
##                     all_alphas[la] = 1 - ( (1-alpha) - sum( 1-alphas[1:(la-1)]))
##                 }else if(type == 'lower'){
##                     all_alphas[la] = alpha - sum(alphas[1:(la-1)])
##                 }else{
##                     stop('Incorrect value of "type"')
##                 }
##                
##                 if(any(all_alphas > 1 | all_alphas < 0)) return(Inf)
##
##                 ks = rep(NA, la)
##                 indiv_rates = rep(NA, la)
##                 for(i in 1:la){
##                     si = in_segment[i]
##                     p1_cs = rate_vs_prob_sorted[[si]]$p1_cumsum
##                     ks[i] = sum(p1_cs < all_alphas[i])
##                     if(ks[i] > 0){
##                         secp = segment_event_cond_prob[[si]][e]
##                         indiv_rates[i] = rate_vs_prob_sorted[[si]]$r1[ks[i]] * secp
##                     }else{
##                         indiv_rates[i] = 0
##                     }
##                 } 
##
##                 return(sum(indiv_rates))
##
##             }
##             ##
##             ## Bound on 0.95 quantile
##             # optimize(robust_quantile, lower=0.95, upper=1, alpha=0.95, in_segment=c(2,3), type='upper', maximum=FALSE)
##             ## Bound on 0.025 quantile
##             # optimize(robust_quantile, lower=0., upper=0.025, alpha=0.025, in_segment=c(2,3), type='lower', maximum=TRUE)
##         } 
##     }
## }


#
# Insert a test here, to confirm that rates are written to all files?
#


#
#
#
save.image('compute_rates_all_sources_session.RData')

if(config$edge_correct_event_rates){
    edge_correct$make_plot()
}


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
