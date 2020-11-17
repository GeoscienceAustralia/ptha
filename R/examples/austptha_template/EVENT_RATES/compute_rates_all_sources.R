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
# Quick check on the input file
source('check_sourcezone_parameters.R', local=TRUE)
check_sourcezone_parameter_row_weights(sourcezone_parameters)

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
MW_MIN = config$MW_MIN # 

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
#' an environment, so we have easy access to key variables. 
#'
#' This is a 'wrapper' function which uses a low-memory-demand "mock" if the
#' sourcezone_parameters_row$row_weight == 0. This saves time/memory when source-zones
#' have been replaced with others using the 'row_weight=0' approach.
#'
#' @param sourcezone_parameters_row data.frame a single row from sourcezone_parameters table
#' @param unsegmented_edge_rate_multiplier a constant, by which we increase the
#'   conditional probability of ruptures on the edge of the source-zone. Used to
#'   make the spatial distribution of moment release on the source better
#'   approximate the desired value. If NULL, it is computed, or ignored (see
#'   also config$edge_correct_event_rates)
#' @return the function environment
#'
source_rate_environment_fun<-function(sourcezone_parameters_row, unsegmented_edge_rate_multiplier=NULL){

    sourcezone_parameters_row = sourcezone_parameters_row
    unsegmented_edge_rate_multiplier = unsegmented_edge_rate_multiplier

    
    if(as.numeric(sourcezone_parameters_row$row_weight) == 0){
        # Shortcut for the case where the rate will always be zero!
        out_env = source_rate_environment_fun_row_weight_zero(sourcezone_parameters_row, unsegmented_edge_rate_multiplier)
    }else{
        # Regular case
        out_env = source_rate_environment_fun_standard(sourcezone_parameters_row, unsegmented_edge_rate_multiplier)
    }

    return(out_env)

}

#' Function to evaluate the rates for a given source-zone. This function returns
#' it's environment, so we have easy access to key variables.
#'
#' This is the main workhorse function, used when (row_weight != 0)
#'
#' @param sourcezone_parameters_row data.frame a single row from sourcezone_parameters table
#' @param unsegmented_edge_rate_multiplier constant by which we increase the
#'   conditional probability of ruptures on the edge of the source-zone. Used to
#'   make the spatial distribution of moment release on the source better
#'   approximate the desired value. If NULL, it is computed, or ignored (see
#'   also config$edge_correct_event_rates)
#' @return the function environment
#'
source_rate_environment_fun_standard<-function(sourcezone_parameters_row, unsegmented_edge_rate_multiplier=NULL){
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
    if(config$coupling_prior_type == 'uniform'){

        sourcepar$coupling  = config$uniform_coupling_prior 

        if(any(sourcepar$coupling == 0)) stop('Cannot treat coupling of zero, instead set prob_Mmax_below_Mmin > 0')

        # Log spacing to equally resolve all coupling values in a relative sense
        # However, we assign probabilities consistent with a uniform distribution
        # (the log spacing is to improve numerical aspects of the discretization,
        # we are not claiming that log(coupling) is uniformly distributed!)
        sourcepar$coupling = exp(approx(log(as.numeric(sourcepar$coupling)), 
            n=nbins*config$coupling_subsampling_increase)$y)
        # Discretize the probabilities as a uniform distribution (NOT log-uniform!)
        # This choice will ensure the numerical derivatives are constant (i.e. uniform density)
        # On the other hand, the discretization will lead to a small numerical difference between
        # the mean coupling, and the mean of a uniform distribution -- although this 'converges away' in
        # the limit of nbins-->Inf
        bin_size = sourcepar$coupling
        sourcepar$coupling_p = bin_size/sum(bin_size)

    }else if(config$coupling_prior_type == 'spreadsheet'){
        # Just use the spreadsheet values
        sourcepar$coupling = sourcezone_parameters_row[1, c('cmin', 'cpref', 'cmax')]
        sourcepar$coupling_p = rep(1, length(sourcepar$coupling))
        sourcepar$coupling_p = sourcepar$coupling_p / sum(sourcepar$coupling_p)

    }else if(config$coupling_prior_type == 'spreadsheet_and_uniform_50_50'){
        # 50% on spreadsheet values, 50% on uniform prior

        pr1 = config$uniform_coupling_prior
        pr2 = sourcezone_parameters_row[1, c('cmin', 'cpref', 'cmax')] 
    
        # Empirical CDF for the uniform coupling prior
        uniform_ecdf     = approxfun(pr1, seq(0, 1, len=length(pr1)), rule=2)
        spreadsheet_ecdf = approxfun(pr2, seq(0, 1, len=length(pr2)), rule=2)

        final_ecdf<-function(x) 0.5*(uniform_ecdf(x) + spreadsheet_ecdf(x))
    
        # Closer together when the values are small 
        range_c = range(c(pr1, pr2)) 
        coupling_vals = exp(approx(log(range_c), n=nbins*config$coupling_subsampling_increase)$y)

        ll = length(coupling_vals)
        # Use numerical derivative to get the density
        coupling_forward = 0.5*(coupling_vals + c(coupling_vals[-1], 2*coupling_vals[ll] - coupling_vals[ll-1]))
        coupling_backward = 0.5*(coupling_vals + c(2*coupling_vals[1] - coupling_vals[2], coupling_vals[1:(ll-1)]))
        coupling_p = final_ecdf(coupling_forward) - final_ecdf(coupling_backward)

        coupling_p = coupling_p / sum(coupling_p)
    
        sourcepar$coupling = coupling_vals
        sourcepar$coupling_p = coupling_p

    }else{
        stop('unrecognized coupling prior type')
    }


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
        # normal (rake = -90). This is OK for PTHA18 but might need to be generalised in future
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
    # Get the event table. We need to pass this to the rate function (so it can
    # moment balance)
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
    # Treatment of variable shear modulus (mathematically, it's like a 'magnitude-observation-error')
    #
    include_variable_shear_modulus = (target_rake == 90)
    if(include_variable_shear_modulus){
        # Treatment of variable shear modulus 
        #
        # Here we compute the conditional cumulative distribution function of
        # the magnitude-difference (i.e. between the 'variable shear modulus'
        # magnitude, and the 'fixed shear modulus' magnitude). It is
        # conditional on the 'fixed shear modulus' magnitude. By treating these
        # differences as 'errors in observed magnitudes', we get a convenient
        # method of treating variable shear modulus in the analysis.
        #
        # WE USE HETEROGENEOUS SLIP MODELS TO COMPUTE THE magnitude-difference
        # CDF IN ALL CASES! 
        # Why? The empirical CDF of magnitude-difference will actually differ
        # depending on whether we compute the function using
        # uniform/variable-uniform/heterogeneous slip scenarios.
        # However, experimentation suggests the impact of this variation on our
        # rate models is small. Heuristically, this is because the CDF affects
        # the 'observed magnitude rates' via a convolution integral, so
        # deviations are 'averaged out' in that process.
        # It is much simpler in terms of our workflow to use a single model, so
        # that is implemented here
        # The empirical CDF based on uniform-slip-with-fixed-size has
        # clearer 'discretization' artefacts due to reduced variability. So we
        # use stochastic slip to derive the empirical CDF, as the large number
        # of events + natural variability leads to a nicely behaved function.
        #
        stochastic_slip_event_table_file = paste0('../SOURCE_ZONES/', source_name, 
        '/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_', source_name, '.nc')
        stochastic_slip_fid = nc_open(stochastic_slip_event_table_file)
        stoc_mw_mu_constant = ncvar_get(stochastic_slip_fid, 'Mw')
        stoc_mw_mu_constant = round(stoc_mw_mu_constant, 3) # Address imperfect floating point storage in netcdf
        stoc_mw_mu_variable = ncvar_get(stochastic_slip_fid, 'variable_mu_Mw')

        # Difference between variable vs uniform shear modulus
        mw_obs_deviation = stoc_mw_mu_variable - stoc_mw_mu_constant
        # Sanity check (specific to our case)
        stopifnot( (min(mw_obs_deviation) > -0.32) & (max(mw_obs_deviation) < 0.25) )

        # If we have segmentation, then we should only look at events which touch the segment
        to_keep = rep(TRUE, length(stoc_mw_mu_variable))
        if(is_a_segment){
            # Read the event indices, and identify events that do not touch the current segment
            stoc_eis = ncvar_get(stochastic_slip_fid, 'event_index_string')
            stoc_eis = sapply(stoc_eis, f<-function(x) as.numeric(strsplit(x, '-')[[1]]), simplify=FALSE)
            for(ei in 1:length(stoc_mw_mu_variable)){
                inds = stoc_eis[[ei]]
                # Record events that are not in the current segment.
                if(!any(is_in_segment[inds])) to_keep[ei] = FALSE
            }
            # Save memory
            rm(stoc_eis)
        }
        nc_close(stochastic_slip_fid); rm(stochastic_slip_fid)
        # Restrict CDF based on events in the segment 
        keepers = which(to_keep)
        if(length(keepers) == 0) stop('Error: No events in segment. This suggests a bug')
        mw_obs_deviation = mw_obs_deviation[keepers]
        stoc_mw_mu_constant = stoc_mw_mu_constant[keepers] 
        # Conditional empirical CDF of difference between 'Mw observation' and
        # 'Mw with constant shear modulus' 
        mw_deviation_cdf_variable_shear_modulus = make_conditional_ecdf(mw_obs_deviation, 
            stoc_mw_mu_constant)

        # Reduce memory
        rm(mw_obs_deviation, stoc_mw_mu_constant, stoc_mw_mu_variable, keepers, to_keep); gc()

    }else{
        # For normal fault events, ignore shear modulus variations
        # Often we only have one 
        mw_deviation_cdf_variable_shear_modulus = NULL
    }
    

    #
    # PART 3
    # Get a preliminary event conditional probability model; use it to compute
    # the mw_rate_function, then update the event conditional probability model
    # to deal with edge-effects (this will not affect the mw_rate_function)
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

            # Note for uniform-slip models with fixed rigidity, inverse-slip is
            # proportional to area (with magnitude fixed).
            # That is the case treated here (notwithstanding perturbations to
            # deal with variable slip or rigidity).
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
    # release more towards the centre of the rupture, because by construction we have less
    # events touching unit-sources at the source-zone along-strike extremes
    # than we have events touching the middle of the source-zone.
    event_conditional_probabilities = get_event_probabilities_conditional_on_Mw(
        event_table, 
        conditional_probability_model = conditional_probability_model)    

    #
    # Get earthquake magnitude data for rate function update
    #
    if(nrow(gcmt_data) > 0){
        gcmt_data_for_rate_function = list(
            # Magnitude
            Mw = gcmt_data$Mw, 
            # Time since the start of observation, in years
            #t = (gcmt_data$julianDay1900 - 
            #    gcmt_access$cmt_start_time_julianDay1900)/gcmt_access$days_in_year
            t = NULL # Censored likelihood is biased for rates, better to use poisson count approach.
        )

        if(min(sourcepar$coupling) == 0){
            stop('Cannot have zero coupling logic tree branch when GCMT data is present')
        }
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
        # Note we pass the 'full' event table, but in segments some of these
        # events will have conditional probability of zero
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

    #
    # Get the fraction of the logic-tree which places non-zero weight on the
    # event being possible. This is nice for heuristically describing epistemic
    # uncertainties.
    #
    # Note that since the event rates are evaluated as GR(Mw-dMw/2) - GR(Mw+dMw/2),
    # it makes most sense to evaluate this term at Mw-dMw/2
    #
    # For segments, we need to be careful to not 'double-count' events which
    # are in 2 segments. In that case the conditional probability model is
    # partially weighted in each, and we should make sure the
    # weight_with_nonzero_rate is as well
    #
    if(is_a_segment){
        fraction_in_segment = sapply(event_table$event_index_string, f<-function(x){
            inds = as.numeric(strsplit(x, '-')[[1]])
            output = mean(is_in_segment[inds])
            return(output)
        })
    }else{
        fraction_in_segment = 1
    }
    event_weight_with_nonzero_rate = (event_conditional_probabilities > 0) * 
        fraction_in_segment *
        mw_rate_function(event_table$Mw-dMw/2, epistemic_nonzero_weight=TRUE)
    # As above with variable mu
    event_weight_with_nonzero_rate_mu_vary = (event_conditional_probabilities > 0) * 
        fraction_in_segment *
        mw_rate_function(event_table$Mw-dMw/2, epistemic_nonzero_weight=TRUE, account_for_mw_obs_error=TRUE)


    #
    # Quantiles below here
    #
    # Note that on source-zones with optional segmentation, we will have to change
    # these quantiles later. 
    #
    # The idea is that "for the segments", the uncertainties
    # should be co-monotinic (i.e. if 16th percentile is true on segmentA, then
    # 16th percentile is also true on segmentB, etc). So the 16th percentile exceedance
    # rate for 'all segments combined' can be obtained by adding their 16th percentiles. 
    #
    # However, we will want to compute quantiles resulting from the "weighted sum of
    # the unsegmented treatment, and the segmented treatment". For example, suppose
    # we give 50% weight to the segmented option (all segments, co-monotonic), and 
    # 50% weight to the unsegmented option. In this case, the 16th percentile for the
    # combination might not be based on combining the 16th percentile for each option
    # separately. Instead, we might take the 32-percentile from the segmented sources,
    # and no weight from the unsegmented source (this would be correct if the 32-percentile
    # in the segmented case was smaller than the smallest unsegmented percentile). Or we might
    # take 20% from the segmented models, and 12% from unsegmented. The actual
    # number will depend on the details of the source and segments. 
    # 
    #



    # Upper credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    event_rates_upper = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=config$upper_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=config$upper_ci_inv_quantile) )
        )

    event_rates_upper_mu_vary = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=config$upper_ci_inv_quantile, 
            account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=config$upper_ci_inv_quantile, 
            account_for_mw_obs_error=TRUE) )
        )

    # Lower credible interval bound. Wrap in as.numeric to avoid having a 1
    # column matrix as output
    event_rates_lower = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=config$lower_ci_inv_quantile) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=config$lower_ci_inv_quantile) )
        )
    event_rates_lower_mu_vary = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=config$lower_ci_inv_quantile, 
            account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=config$lower_ci_inv_quantile, 
            account_for_mw_obs_error=TRUE) )
        )


    #
    # Median
    #
    event_rates_median = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.5) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.5) )
        )
    event_rates_median_mu_vary = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.5, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.5, account_for_mw_obs_error=TRUE) )
        )

    #
    # 16pc quantile
    #
    event_rates_16pc = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.16) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.16) )
        )
    event_rates_16pc_mu_vary = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.16, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.16, account_for_mw_obs_error=TRUE) )
        )

    #
    # 84 pc quantile
    #
    event_rates_84pc = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.84) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.84) )
        )
    event_rates_84pc_mu_vary = as.numeric(event_conditional_probabilities * 
        (mw_rate_function(event_table$Mw - dMw/2, quantiles=0.84, account_for_mw_obs_error=TRUE) - 
         mw_rate_function(event_table$Mw + dMw/2, quantiles=0.84, account_for_mw_obs_error=TRUE) )
        )

    gc()


    return(environment())

}



#' Compute source-zone variables efficiently when row_weight=0
#'
#' To allow dropping/replacing source-zones cleanly, we allow row_weight=0. This means
#' that all events on the source-zone will have a rate=0, i.e. it will not contribute
#' to the final hazard. 
#'
#' In that instance the call to source_rate_environment_fun is overly expensive
#' and demanding of memory. However, to run the code without too many changes, we really
#' need some variables to exist even when row_weight=0.
#'
#' This function makes those variables, efficiently, in the 'row_weight=0' case
#'
source_rate_environment_fun_row_weight_zero<-function(
    sourcezone_parameters_row, 
    unsegmented_edge_rate_multiplier=NULL){

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

    # Rake
    target_rake = bird2003_env$unit_source_tables[[source_name]]$rake[1]

    # Edge multiplier
    best_edge_mult = list()
    best_edge_mult$minimum = 0

    #
    # Get the event table. 
    #
    event_table_file = paste0('../SOURCE_ZONES/', source_name, 
        '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_', 
        source_name, '.nc')
    event_table = read_table_from_netcdf(event_table_file)

    # Vector of zeros, reused as empty data elsewhere
    empty_data = rep(0, length(event_table$Mw))

    #
    # 'ZERO' version of the variables that are needed later on,
    # in order to 'cleanly' write zeros to the netcdf files.
    #
    event_rates = empty_data
    event_weight_with_nonzero_rate = empty_data
    event_rates_median = empty_data
    event_rates_16pc = empty_data
    event_rates_84pc = empty_data
    event_rates_upper = empty_data
    event_rates_lower = empty_data

    event_rates_mu_vary = empty_data
    event_weight_with_nonzero_rate_mu_vary = empty_data 
    event_rates_median_mu_vary = empty_data
    event_rates_16pc_mu_vary = empty_data
    event_rates_84pc_mu_vary = empty_data
    event_rates_upper_mu_vary = empty_data
    event_rates_lower_mu_vary = empty_data

    # The above variables should allow the code to execute cleanly in all cases!

    return(environment())
}


#' Suppose a source-zone has some weight on an unsegmented model, and the remaining
#' weight on a collection of segments.
#'
#' How should we compute percentiles (or credible-intervals) of the 'combined' Mw-frequency
#' curve (i.e. to characterise uncertainty)?
#'
#' The approach taken herein is:
#'    A) Firstly take ALL the segmented sources. Their combined Mw-frequency curve
#'       is developed, with percentiles based on assuming uncertainties are co-monotonic. 
#'
#'       For illustration, suppose we are interested in a percentile of the exceedance 
#'       rate of Mw 9 events (e.g. 16th percentile).
#'       On the combined-segments, this is computed assuming the segments are 'co-monotonic'
#'       (i.e. as the sum of the 16th percentile of the exceedance rates on each segment).
#'
#'    B) Once (A) is completed, we have 2 models for the source-zone, one segmented and one unsegmented.
#'       These will have been assigned row_weights which sum to 1.0 (e.g. 0.5 and 0.5 is common).
#'       How to compute the 16th percentile (or any other) of the Mw 9.0 exceedance rate for this combination?
#'       The solution is purely mathmatical -- there is a 50% chance of segmented or unsegmented being true, and
#'       so the 16th percentile should be taken from the 16th percentile of the combined distribution.
#'       In practice, this might mean we evaluate the segmented models at their 32th percentile and ignore
#'       the unsegmented model (this would be correct if the 32th segmented percentile was lower than the 0th 
#'       unsegmented percentile). Or maybe 22th segmented and 10th segmented. It depends on the inputs, but the
#'       answer can be calculated. 
#'        
#' This function does the above calculations, and creates new variables in the environment which
#' can map the 'rate percentiles' onto individual scenarios, which go into the netcdf files.
#' 
#' @param all_sources a list of environments corresponding to all segmented and unsegmented sources on 
#' a single source-zone. Each is an output of compute_rates_all_sources. The first one should be the 
#' unsegmented source, followed by segments.
#' @param percentile_discretization small real number. The function works by numerically discretizing the mw-rate-function
#' at a range of logic-tree percentiles (i.e. inverse quantiles), with spacing = dp which is close to percentile_discretization,
#' while ensuring that 1/dp is an integer. The discretization sequence goes from dp/2 to 1-dp/2.
#' @param quick_exit_if_row_weights_all_zero Logical - Skip the expensive calculations if the row weights are all zero
update_scenario_rate_percentiles_on_source_zones_with_partial_segmentation<-function(all_sources,
    percentile_discretization = 0.0025, quick_exit_if_row_weights_all_zero=TRUE){


    for(i in 1:length(all_sources)){
        #
        # Later, we will have to update the percentiles to get consistent treatment on partly
        # segmented source-zones. These variables will let us do that cleanly.
        #
        all_sources[[i]]$file_event_rates_median = all_sources[[i]]$event_rates_median
        all_sources[[i]]$file_event_rates_16pc = all_sources[[i]]$event_rates_16pc
        all_sources[[i]]$file_event_rates_84pc = all_sources[[i]]$event_rates_84pc
        all_sources[[i]]$file_event_rates_upper = all_sources[[i]]$event_rates_upper
        all_sources[[i]]$file_event_rates_lower = all_sources[[i]]$event_rates_lower
        all_sources[[i]]$file_event_rates_median_mu_vary = all_sources[[i]]$event_rates_median_mu_vary
        all_sources[[i]]$file_event_rates_16pc_mu_vary = all_sources[[i]]$event_rates_16pc_mu_vary
        all_sources[[i]]$file_event_rates_84pc_mu_vary = all_sources[[i]]$event_rates_84pc_mu_vary
        all_sources[[i]]$file_event_rates_upper_mu_vary = all_sources[[i]]$event_rates_upper_mu_vary
        all_sources[[i]]$file_event_rates_lower_mu_vary = all_sources[[i]]$event_rates_lower_mu_vary
    }

    # No need to update unsegmented sources
    if(length(all_sources) == 1){ 
        # Double check the input makes sense
        if(all_sources[[1]]$is_a_segment){
            stop('If a source is segmented, then there should be multiple environments in all_sources')
        }
        # Quick exit without doing anything
        return(invisible(0))
    }
  
    # 
    # Below here, we are definitely on a partially segmented source-zone 
    #
    
    is_a_segment = unlist(lapply(all_sources, f<-function(x) x$is_a_segment))
    if(sum(!is_a_segment) != 1) stop('There should be only one unsegmented source in all_sources')
    if(is_a_segment[1]) stop('The first segment in all_sources should be unsegmented')

    # Weight on segmented and unsegmented models.
    unsegmented_weight = all_sources[[1]]$sourcepar$sourcezone_parameters_row$row_weight
    segmented_weight = all_sources[[2]]$sourcepar$sourcezone_parameters_row$row_weight

    # If the row_weights are zero, we can quick-exit here
    if( (unsegmented_weight + segmented_weight == 0) & quick_exit_if_row_weights_all_zero){
        return(invisible(0))
    }

    stopifnot(abs(unsegmented_weight + segmented_weight - 1) < 1.0e-10)

    #
    # Step 1: for each magnitude bin, we need to determine how to compute the
    # percentiles of interest on the 'unsegmented + sum-of-segments' Mw-rate curve.
    # We do this by numerically evaluating the Mw-rate curves, including
    # uncertainties, on the unsegmented and 'sum-of-segments' separately. Then
    # we determine the appropriate 'segmented' and 'unsegmented' percentiles that
    # lead to the desired 'overall percentile'. 
    #
   
    # Make a sequence of mw values where we will want to evaluate the
    # mw-exceedance rate curve and deal with the percentile estimation issue 
    lower_mw = min(all_sources[[1]]$event_table$Mw) - dMw/2 
    upper_mw = max(all_sources[[1]]$event_table$Mw) + dMw/2
    mw_seq = seq(lower_mw, upper_mw, len=round((upper_mw - lower_mw)/dMw + 1)) #by=dMw)

    # Interpolate credible intervals at the following percentiles for the unsegmented
    # sources, and the segments as a whole 
    # Useful to numerically figure out how to evaluate the percentiles
    dp = 1/round(1/percentile_discretization) # Ensures 1/dp is an integer
    percentiles_to_store = seq(dp/2, 1-dp/2, len=round(1/dp))

    # Make a matrix with the mw_seq vs exceedance-rates for every percentile_to_store
    unsegmented_mw_rate_curves = all_sources[[1]]$mw_rate_function(mw_seq, 
        quantiles=percentiles_to_store)
    unsegmented_mw_rate_curves_mu_vary = all_sources[[1]]$mw_rate_function(mw_seq, 
        quantiles=percentiles_to_store, account_for_mw_obs_error=TRUE)
    gc()

    # As above, for each segment individually
    segmented_mw_rate_curves_list = list()
    segmented_mw_rate_curves_mu_vary_list = list()
    for(i in 2:length(all_sources)){
        segmented_mw_rate_curves_list[[i-1]] = all_sources[[i]]$mw_rate_function(mw_seq, 
            quantiles=percentiles_to_store)
        segmented_mw_rate_curves_mu_vary_list[[i-1]] = all_sources[[i]]$mw_rate_function(mw_seq, 
            quantiles=percentiles_to_store, account_for_mw_obs_error=TRUE)
        gc()
    }

    # Make the 'sum of segments'. Because of our co-monotonic assumption, this
    # can be done with simple summation
    segmented_mw_rate_curves = segmented_mw_rate_curves_list[[1]]
    segmented_mw_rate_curves_mu_vary = segmented_mw_rate_curves_mu_vary_list[[1]]
    if(length(segmented_mw_rate_curves_list) > 1){
        for(i in 2:length(segmented_mw_rate_curves_list)){
            segmented_mw_rate_curves = segmented_mw_rate_curves + segmented_mw_rate_curves_list[[i]]
            segmented_mw_rate_curves_mu_vary = segmented_mw_rate_curves_mu_vary + segmented_mw_rate_curves_mu_vary_list[[i]]
        }
    }

    # For each percentile of interest, we need to compute the percentile at which we evaluate the
    # unsegmented rate model, and the segmented rate model, so that we ultimately get the desired 
    # percentile of the combined model
    desired_inv_quantiles = c(config$lower_ci_inv_quantile, 0.16, 0.5, 0.84, config$upper_ci_inv_quantile)
    names_desired_inv_quantiles = c('lower', '16pc', 'median', '84pc', 'upper')
    # Store the results here
    all_pc_values = list()
    for(i in 1:length(desired_inv_quantiles)){
        all_pc_values[[i]] = list()
        all_pc_values[[i]]$seg = rep(0, length(mw_seq))
        all_pc_values[[i]]$unseg = rep(0, length(mw_seq))
    }
    names(all_pc_values) = names_desired_inv_quantiles
    all_pc_values_mu_vary = all_pc_values
   
    for(i in 1:length(mw_seq)){

        #
        # Get percentiles of the rate curve for magnitude=mw_seq[i]
        # Constant shear modulus
        #
        rate_unseg = unsegmented_mw_rate_curves[i,]
        rate_seg = segmented_mw_rate_curves[i,]
        # Get the 'uncertainty percentile' in the combined unsegmented-segmented curve associated
        # with each rate value in unique_rates. 
        unique_rates = unique(sort(c(0, rate_unseg, rate_seg))) # Ensure we compare with a rate of 0, for numerical reasons
        pc_value = rep(0, length(unique_rates))
        pc_value_unseg = rep(0, length(unique_rates))
        pc_value_seg   = rep(0, length(unique_rates))
        for(j in 1:length(unique_rates)){
            # Percentile of the unseg curve associated with unique_rates[j]
            pc_value_unseg[j] = max( percentiles_to_store * (rate_unseg <= unique_rates[j]) )
            # Percentile of the seg curve associated with unique_rates[j]
            pc_value_seg[j]   = max( percentiles_to_store * (rate_seg   <= unique_rates[j]) )
            # Equation for the 'percentile of the combined distribution'
            pc_value[j] = pc_value_unseg[j]*unsegmented_weight + pc_value_seg[j]*segmented_weight
        }

        #
        # As above, variable shear modulus
        #
        rate_unseg_mu_vary = unsegmented_mw_rate_curves_mu_vary[i,]
        rate_seg_mu_vary = segmented_mw_rate_curves_mu_vary[i,]
        unique_rates = unique(sort(c(0, rate_unseg_mu_vary, rate_seg_mu_vary))) # Ensure we compare with rate = 0 for numerical reasons
        pc_value_mu_vary = rep(0, length(unique_rates))
        pc_value_unseg_mu_vary = rep(0, length(unique_rates))
        pc_value_seg_mu_vary   = rep(0, length(unique_rates))
        for(j in 1:length(unique_rates)){
            # Percentile of the unseg curve associated with unique_rates[j]
            pc_value_unseg_mu_vary[j] = max( percentiles_to_store * (rate_unseg_mu_vary <= unique_rates[j]) )
            # Percentile of the seg curve associated with unique_rates[j]
            pc_value_seg_mu_vary[j]   = max( percentiles_to_store * (rate_seg_mu_vary   <= unique_rates[j]) )
            # Equation for the 'percentile of the combined distribution'
            pc_value_mu_vary[j] = pc_value_unseg_mu_vary[j]*unsegmented_weight + pc_value_seg_mu_vary[j]*segmented_weight
        }


        # Handy function to pick out the segmented/unsegmented percentile values
        # associated with the actual value we want for their combination
        #
        # Note eps is used to deal with floating point imperfections
        get_seg_unseg_pc<-function(desired_inv_quantile, pc_value, pc_value_seg, pc_value_unseg, eps=1e-08){
    
            thresh_index = sum(pc_value <= desired_inv_quantile + eps)
          
            # Later NA will be interpreted as 'set rates to zero' 
            seg_unseg_output = rep(NA, 2)

            # Note if thresh_index = 1, then the rate is 0 (by construction
            # above), so we can 'ignore' these ones
            if(thresh_index > 1){
                # Typical case
                seg_unseg_output[1] = pc_value_seg[thresh_index]
                seg_unseg_output[2] = pc_value_unseg[thresh_index]
            }

            return(seg_unseg_output) 
        }

        # Store the results at the desired percentiles
        eps = 1.0e-8 # Deal with floating point imperfections
        for(j in 1:length(desired_inv_quantiles)){
            # Constant mu
            tmp = get_seg_unseg_pc(desired_inv_quantiles[j], pc_value, pc_value_seg, 
                pc_value_unseg, eps)
            all_pc_values[[j]]$seg[i] = tmp[1]
            all_pc_values[[j]]$unseg[i] = tmp[2]
            # Variable mu
            tmp = get_seg_unseg_pc(desired_inv_quantiles[j], pc_value_mu_vary, 
                pc_value_seg_mu_vary, pc_value_unseg_mu_vary, eps)
            all_pc_values_mu_vary[[j]]$seg[i] = tmp[1]
            all_pc_values_mu_vary[[j]]$unseg[i] = tmp[2]
        }

    }
    # At this stage we know the percentiles at which we should evaluate the
    # segmented and unsegmented distributions, in order to get the 'combined'
    # distribution at the desired percentiles (=desired_inv_quantiles)


    # STEP 1.5
    #
    # Compute the exceedance rate curves for each individual segment (as well as the unsegmented case)
    #
    # We also check that the 'individual segment' rate curves never increase.
    # (which would locally lead to a negative rate increment). This can happen
    # (although rare) because the percentile at which we evaluate the exceedance rate varies with mw. 
    # Furthermore, above we did the analysis based on the 'sum of segments', which does not
    # guarentee monotonicity for an individual segment (although non-monotonicity is rare).
    # We definitely do not want to create negative scenario rate increments anywhere!
    # Better to slightly distort the quantile. 
    #
    ex_rates = vector(mode='list', length=length(all_sources)) # Store exceedance rates at varying mw/percentiles
    ex_rates_mu_vary = vector(mode='list', length=length(all_sources))
    for(j in 1:length(desired_inv_quantiles)){

        nm = names_desired_inv_quantiles[j]

        for(i in 1:length(all_sources)){

            # Get the percentile that matches desired_inv_quantiles[j]
            if(i == 1){
                # Unsegmented
                pc_indices = match(all_pc_values[[nm]]$unseg, percentiles_to_store)
                pc_indices_mu_vary = match(all_pc_values_mu_vary[[nm]]$unseg, percentiles_to_store)
            }else{
                # Segmented
                pc_indices = match(all_pc_values[[nm]]$seg, percentiles_to_store)
                pc_indices_mu_vary = match(all_pc_values_mu_vary[[nm]]$seg, percentiles_to_store)
            }

            # Get the exceedence rates (segment specific)
            ex_rates[[i]][[nm]] = rep(0, len=length(pc_indices))
            ex_rates_mu_vary[[i]][[nm]] = rep(0, len=length(pc_indices))
            for(k in 1:length(pc_indices)){
    
                # Constant mu
                pc_ind = pc_indices[k]
                if(!is.na(pc_ind)){ 

                    if(i == 1){
                        # Unsegmented
                        ex_rates[[i]][[nm]][k] = unsegmented_mw_rate_curves[k, pc_ind]
                    }else{
                        # Segmented
                        ex_rates[[i]][[nm]][k] = segmented_mw_rate_curves_list[[i-1]][k, pc_ind]
                    }
                }

                # Variable mu
                pc_ind = pc_indices_mu_vary[k]
                if(!is.na(pc_ind)){ 

                    if(i == 1){
                        # Unsegmented
                        ex_rates_mu_vary[[i]][[nm]][k] = unsegmented_mw_rate_curves_mu_vary[k, pc_ind]
                    }else{
                        # Segmented
                        ex_rates_mu_vary[[i]][[nm]][k] = segmented_mw_rate_curves_mu_vary_list[[i-1]][k, pc_ind]
                    }
                }

                # Check there is no increase in the rate! This can happen (although rarely in practice).
                # If the rate curve was always evaluated at the same percentile, then this would never happen.
                # However, for this problem, the percentile changes with mw at the individual segment level
                # (such that the percentile is fixed for the combined source exceedance rate curve).
                # An increate in the rate for a single segment would be a
                # problem, because then we would add a 'negative scenario rate' for the segment,
                # and there is no guarentee that will be cancelled by the other sources in every case.
                # To solve it we slightly distort the exceedance rate at the quantile. Considering approximations
                # in these calculations it should not be a big deal.
                # If we get isolated warnings like this there should be no problem, but many warnings suggest trouble.
                if(k > 1){

                    #
                    # Constant shear modulus case
                    #
                    if(ex_rates[[i]][[nm]][k] > ex_rates[[i]][[nm]][k-1]){
            
                        # Try increasing the previous rate
                        test_rate = ex_rates[[i]][[nm]][k]
                        if(k > 2){
                            limit = ex_rates[[i]][[nm]][k-2] 
                        }else{
                            limit = Inf
                        }

                        if(test_rate <= limit){
                            # We can increase the previous rate. This is preferable if it's possible,
                            # because otherwise if the previous rate is zero we
                            # will force zeros for all following k
                            message = paste0('Warning: Adjusting percentile ', nm, ' at Mw=', mw_seq[k-1], 
                                ' on ', all_sources[[i]]$source_segment_name, ' with constant mu \n ', 
                            ' Changing from ', ex_rates[[i]][[nm]][k-1], ' to ', test_rate)
                            print(message)
                            ex_rates[[i]][[nm]][k-1] = test_rate
                        }else{
                            message = paste0('Warning: Adjusting percentile ', nm, ' at Mw=', mw_seq[k], 
                                ' on ', all_sources[[i]]$source_segment_name, ' with constant mu \n ', 
                            ' Changing from ', ex_rates[[i]][[nm]][k], ' to ', ex_rates[[i]][[nm]][k-1])
                            print(message)
                            ex_rates[[i]][[nm]][k] = ex_rates[[i]][[nm]][k-1]
                        }
                    } 

                    #
                    # Same as above
                    # Variable shear modulus case
                    #
                    if(ex_rates_mu_vary[[i]][[nm]][k] > ex_rates_mu_vary[[i]][[nm]][k-1]){
            
                        # Try increasing the previous rate
                        test_rate = ex_rates_mu_vary[[i]][[nm]][k]
                        if(k > 2){
                            limit = ex_rates_mu_vary[[i]][[nm]][k-2] 
                        }else{
                            limit = Inf
                        }

                        if(test_rate <= limit){
                            # We can increase the previous rate. This is preferable if it's possible,
                            # because otherwise if the previous rate is zero we
                            # will force zeros for all following k
                            message = paste0('Warning: Adjusting percentile ', nm, ' at Mw=', mw_seq[k-1], 
                                ' on ', all_sources[[i]]$source_segment_name, ' with variable mu \n ', 
                            ' Changing from ', ex_rates_mu_vary[[i]][[nm]][k-1], ' to ', test_rate)
                            print(message)
                            ex_rates_mu_vary[[i]][[nm]][k-1] = test_rate
                        }else{
                            message = paste0('Warning: Adjusting percentile ', nm, ' at Mw=', mw_seq[k], 
                                ' on ', all_sources[[i]]$source_segment_name, ' with variable mu \n ', 
                            ' Changing from ', ex_rates_mu_vary[[i]][[nm]][k], ' to ', ex_rates_mu_vary[[i]][[nm]][k-1])
                            print(message)
                            ex_rates_mu_vary[[i]][[nm]][k] = ex_rates_mu_vary[[i]][[nm]][k-1]
                        }
                    } 

                }
            }
        }
    }


    
    # STEP 2:
    # Now update the file_event_rates_XX variables, using new percentile values
    # which will ensure that when the rates are summed over unsegmented and segmented models,
    # we will hit the target percentile for the source-zone as a whole.
    #
    for(i in 1:length(all_sources)){

        for(j in 1:(length(mw_seq)-1) ){
    
            # Operate "magnitude by magnitude"
            k = which( (all_sources[[i]]$event_table$Mw > mw_seq[j]  ) & 
                       (all_sources[[i]]$event_table$Mw < mw_seq[j+1]) )
            stopifnot( all( all_sources[[i]]$event_table$Mw[k] == all_sources[[i]]$event_table$Mw[k[1]]) )
 
            # Median, constant mu
            delta_rate = -diff(ex_rates[[i]][['median']][j:(j+1)])  
            all_sources[[i]]$file_event_rates_median[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # Median, variable mu
            delta_rate = -diff(ex_rates_mu_vary[[i]][['median']][j:(j+1)]) 
            all_sources[[i]]$file_event_rates_median_mu_vary[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # Lower, constant mu
            delta_rate = -diff(ex_rates[[i]][['lower']][j:(j+1)])
            all_sources[[i]]$file_event_rates_lower[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # Lower, variable mu
            delta_rate = -diff(ex_rates_mu_vary[[i]][['lower']][j:(j+1)]) 
            all_sources[[i]]$file_event_rates_lower_mu_vary[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # 16pc, constant mu
            delta_rate = -diff(ex_rates[[i]][['16pc']][j:(j+1)])
            all_sources[[i]]$file_event_rates_16pc[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # 16pc, variable mu
            delta_rate = -diff(ex_rates_mu_vary[[i]][['16pc']][j:(j+1)])
            all_sources[[i]]$file_event_rates_16pc_mu_vary[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # 84pc, constant mu
            delta_rate = -diff(ex_rates[[i]][['84pc']][j:(j+1)])
            all_sources[[i]]$file_event_rates_84pc[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # 84pc, variable mu
            delta_rate = -diff(ex_rates_mu_vary[[i]][['84pc']][j:(j+1)]) 
            all_sources[[i]]$file_event_rates_84pc_mu_vary[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # upper, constant mu
            delta_rate = -diff(ex_rates[[i]][['upper']][j:(j+1)])
            all_sources[[i]]$file_event_rates_upper[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  

            # upper, variable mu
            delta_rate = -diff(ex_rates_mu_vary[[i]][['upper']][j:(j+1)]) 
            all_sources[[i]]$file_event_rates_upper_mu_vary[k] = as.numeric(
                all_sources[[i]]$event_conditional_probabilities[k] * delta_rate )  
        }
    }
   
    return(invisible(0)) 
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

        #
        # Put rates onto the event table nc file for uniform slip,
        # deterministic size, fixed mu events.
        # Notice there is no 'bias correction' applied for uniform slip events,
        # because in practice we do not suggest using them for hazard studies,
        # but instead use them as a 'reference point' to weight the stochastic
        # and variable-uniform slip events
        #
        fid = nc_open(event_table_file, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'rate_annual', event_rates)
        # Summation of credible intervals OK for co-monotonic epistemic uncertainties
        ncvar_put_extra(fid, 'rate_annual_upper_ci', file_event_rates_upper)
        ncvar_put_extra(fid, 'rate_annual_lower_ci', file_event_rates_lower)
        ncvar_put_extra(fid, 'weight_with_nonzero_rate', event_weight_with_nonzero_rate)
        ncvar_put_extra(fid, 'rate_annual_median', file_event_rates_median)
        ncvar_put_extra(fid, 'rate_annual_16pc', file_event_rates_16pc)
        ncvar_put_extra(fid, 'rate_annual_84pc', file_event_rates_84pc)

        #
        # As above with variable shear modulus
        #
        ncvar_put_extra(fid, 'variable_mu_rate_annual', event_rates_mu_vary)
        # Summation of credible intervals OK for co-monotonic epistemic uncertainties
        ncvar_put_extra(fid, 'variable_mu_rate_annual_upper_ci', file_event_rates_upper_mu_vary)
        ncvar_put_extra(fid, 'variable_mu_rate_annual_lower_ci', file_event_rates_lower_mu_vary)
        ncvar_put_extra(fid, 'variable_mu_weight_with_nonzero_rate', event_weight_with_nonzero_rate_mu_vary)
        ncvar_put_extra(fid, 'variable_mu_rate_annual_median', file_event_rates_median_mu_vary)
        ncvar_put_extra(fid, 'variable_mu_rate_annual_16pc', file_event_rates_16pc_mu_vary)
        ncvar_put_extra(fid, 'variable_mu_rate_annual_84pc', file_event_rates_84pc_mu_vary)

        nc_close(fid)

        #
        # Put rates onto the uniform slip table nc file that also has tsunami
        # summary stats. We do not store the variable mu information here (easier
        # implementation, plus it is generally more efficient to read from the
        # non-tsunami file). Also, there is no bias correction for uniform-slip-fixed-size
        #
        event_table_fileB = paste0('../SOURCE_ZONES/', source_name, 
            '/TSUNAMI_EVENTS/all_uniform_slip_earthquake_events_tsunami_',
            source_name, '.nc')
        fid = nc_open(event_table_fileB, readunlim=FALSE, write=TRUE)
        ncvar_put_extra(fid, 'event_rate_annual', event_rates)
        # Summation of credible intervals OK for co-monotonic epistemic uncertainties
        ncvar_put_extra(fid, 'event_rate_annual_upper_ci', file_event_rates_upper)
        ncvar_put_extra(fid, 'event_rate_annual_lower_ci', file_event_rates_lower)
        nc_close(fid)

        #
        # Put rates onto the stochastic and variable uniform slip table nc file
        # For these cases we consider bias correction as well as shear modulus variation
        #
        for(slip_type in c('stochastic', 'variable_uniform')){

            # Get the function which applies bias-correction to event weights
            # based on the peak slip quantiles.
            if(target_rake == 90){
                # Bias adjustment depends on slip type and shear modulus treatment
                if(slip_type == 'stochastic'){
                    bias_adjustment_function = config$peak_slip_bias_adjustment_stochastic
                    bias_adjustment_function_variable_mu = config$peak_slip_bias_adjustment_stochastic_mu_vary
                }else if(slip_type == 'variable_uniform'){
                    bias_adjustment_function = config$peak_slip_bias_adjustment_variable_uniform
                    bias_adjustment_function_variable_mu = config$peak_slip_bias_adjustment_variable_uniform_mu_vary
                }else{
                    stop('slip_type not recognized')
                }
            }else{
                # Do not adjust for normal faults, as we don't have enough test
                # data. In particular we should not treat fixed and variable shear modulus
                # differently because for normal faults we always assume mu=60
                bias_adjustment_function<-function(x) 1 + 0*x
                bias_adjustment_function_variable_mu<-function(x) 1 + 0*x
            }

            #
            # Open the _tsunami file for writing
            # As above, this file does not store variable shear modulus info,
            # mainly because it is inconvenient to add the information (given
            # I did not originally define the file to have it). 
            #
            event_table_fileC = paste0('../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_', slip_type, '_slip_earthquake_events_tsunami_',
                source_name, '.nc')

            fid = nc_open(event_table_fileC, readunlim=FALSE, write=TRUE)

            # We need to 'smear' the uniform-slip-fixed-size rates over the events.
            # To do this we need to know the index corresponding to the uniform slip fixed-size row
            event_uniform_event_row = ncvar_get(fid, 'event_uniform_event_row')
            # Number of events corresponding to event row
            nevents = table(event_uniform_event_row)
            names_nevents = as.numeric(names(nevents))
            # Check that it is making sense!
            stopifnot(all(names_nevents == 1:length(event_rates)))

            # Compute an adjustment to the event weights based on the
            # bias-correction from earlier. Use peak-slip to do this
            event_peak_slip = ncvar_get(fid, 'event_slip_string')
            event_peak_slip = sapply(event_peak_slip, f<-function(x) max(as.numeric(strsplit(x, '_')[[1]])))
            event_bias_adjustment = event_peak_slip * 0

            # Apply a peak slip limitation 
            reference_shear_modulus = sourcezone_parameters_row$shear_modulus * 1e+10
            allowed_peak_slip = config$peak_slip_limit_factor * slip_from_Mw(
                event_table$Mw[event_uniform_event_row], # Magnitude at reference shear modulus
                mu=reference_shear_modulus,
                relation=sourcezone_parameters_row$scaling_relation)

            # Compute the weights
            for(euer in names_nevents){

                k = which(event_uniform_event_row == euer)

                bias_adjuster = rep(0, length(k))
                # We must zero rate for events with peak_slip > allowed value
                acceptable_k = which(event_peak_slip[k] <= allowed_peak_slip[k])
                if(length(acceptable_k) > 0){
                    # Weight individual events, based on the 'bias adjustment' functions devised by
                    # comparing DART buoys with family of model events
                    k2 = k[acceptable_k]
                    quantiles_of_peak_slip = rank(event_peak_slip[k2], ties.method='first')/(length(k2)+1)
                    bias_adjuster[acceptable_k] = bias_adjustment_function(quantiles_of_peak_slip)
                }

                # Normalise so it sums to 1
                sba = sum(bias_adjuster)
                if(sba > 0){
                    bias_adjuster = bias_adjuster/sba
                }
                event_bias_adjustment[k] = bias_adjuster

            }

            ncvar_put_extra(fid, 'event_rate_annual', 
                event_rates[event_uniform_event_row] * event_bias_adjustment)
            # Summation of credible intervals OK for co-monotonic epistemic
            # uncertainties
            ncvar_put_extra(fid, 'event_rate_annual_upper_ci', 
                file_event_rates_upper[event_uniform_event_row] * event_bias_adjustment)
            ncvar_put_extra(fid, 'event_rate_annual_lower_ci', 
                file_event_rates_lower[event_uniform_event_row] * event_bias_adjustment)

            nc_close(fid)

            #
            # Now to the same, for the file that contains only the earthquakes.
            # This file stores variable mu information, and has bias adjustment
            # that varies for the fixed-mu and variable-mu cases
            # 
            event_table_fileD = paste0('../SOURCE_ZONES/', source_name, 
                '/TSUNAMI_EVENTS/all_', slip_type, '_slip_earthquake_events_',
                source_name, '.nc')

            fid = nc_open(event_table_fileD, readunlim=FALSE, write=TRUE)

            # Index corresponding to uniform slip row
            event_uniform_event_row = ncvar_get(fid, 'uniform_event_row')
            # Number of events corresponding to event row
            nevents = table(event_uniform_event_row)
            names_nevents = as.numeric(names(nevents))
            stopifnot(all(names_nevents == 1:length(event_rates)))

            # Compute an adjustment to the event weights based on the
            # bias-correction from earlier
            # Use peak-slip to do this
            event_peak_slip = ncvar_get(fid, 'event_slip_string')
            event_peak_slip = sapply(event_peak_slip, f<-function(x) max(as.numeric(strsplit(x, '_')[[1]])))
            event_bias_adjustment = event_peak_slip * 0
            event_bias_adjustment_variable_mu = event_peak_slip * 0

            # Apply a peak slip limitation 
            reference_shear_modulus = sourcezone_parameters_row$shear_modulus * 1e+10
            allowed_peak_slip = config$peak_slip_limit_factor * slip_from_Mw(
                event_table$Mw[event_uniform_event_row], # Magnitude at reference shear modulus
                mu=reference_shear_modulus,
                relation=sourcezone_parameters_row$scaling_relation)

            # Compute the weights
            for(euer in names_nevents){
                k = which(event_uniform_event_row == euer) 
                bias_adjuster = rep(0, length(k))
                # Same peak_slip criterion for both fixed-and-variable shear-modulus
                acceptable_k = which(event_peak_slip[k] <= allowed_peak_slip[k])

                #
                # Fixed mu case
                #
                if(length(acceptable_k) > 0){
                    k2 = k[acceptable_k]
                    quantiles_of_peak_slip = rank(event_peak_slip[k2], ties.method='first')/(length(k2)+1)
                    bias_adjuster[acceptable_k] = bias_adjustment_function(quantiles_of_peak_slip)
                }

                # Normalise so it sums to 1
                sba = sum(bias_adjuster)
                if(sba > 0){
                    bias_adjuster = bias_adjuster/sba
                }
                event_bias_adjustment[k] = bias_adjuster

                #
                # Variable mu case
                #
                bias_adjuster = bias_adjuster * 0
                if(length(acceptable_k) > 0){
                    # quantiles_of_peak_slip already defined above
                    bias_adjuster[acceptable_k] = bias_adjustment_function_variable_mu(quantiles_of_peak_slip)
                }
                # Normalise so it sums to 1
                sba = sum(bias_adjuster)
                if(sba > 0){
                    bias_adjuster = bias_adjuster/sba
                }
                event_bias_adjustment_variable_mu[k] = bias_adjuster

            }

            ncvar_put_extra(fid, 'rate_annual', 
                event_rates[event_uniform_event_row] * event_bias_adjustment)
            # Summation of credible intervals OK for co-monotonic epistemic uncertainties
            ncvar_put_extra(fid, 'rate_annual_upper_ci', 
                file_event_rates_upper[event_uniform_event_row] * event_bias_adjustment)
            ncvar_put_extra(fid, 'rate_annual_lower_ci', 
                file_event_rates_lower[event_uniform_event_row] * event_bias_adjustment)
            ncvar_put_extra(fid, 'weight_with_nonzero_rate', 
                event_weight_with_nonzero_rate[event_uniform_event_row])
            ncvar_put_extra(fid, 'rate_annual_median', file_event_rates_median[event_uniform_event_row] * event_bias_adjustment)
            ncvar_put_extra(fid, 'rate_annual_16pc', file_event_rates_16pc[event_uniform_event_row] * event_bias_adjustment)
            ncvar_put_extra(fid, 'rate_annual_84pc', file_event_rates_84pc[event_uniform_event_row] * event_bias_adjustment)


            #
            # Put the variable_mu rates onto the file. Note they do not
            # go onto the '_tsunami' file (just the 'events only' file)
            #
            ncvar_put_extra(fid, 'variable_mu_rate_annual', 
                event_rates_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)
            # Summation of credible intervals OK for co-monotonic epistemic uncertainties
            ncvar_put_extra(fid, 'variable_mu_rate_annual_upper_ci', 
                file_event_rates_upper_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)
            ncvar_put_extra(fid, 'variable_mu_rate_annual_lower_ci', 
                file_event_rates_lower_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)
            ncvar_put_extra(fid, 'variable_mu_weight_with_nonzero_rate', 
                event_weight_with_nonzero_rate_mu_vary[event_uniform_event_row])
            ncvar_put_extra(fid, 'variable_mu_rate_annual_median', 
                file_event_rates_median_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)
            ncvar_put_extra(fid, 'variable_mu_rate_annual_16pc', 
                file_event_rates_16pc_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)
            ncvar_put_extra(fid, 'variable_mu_rate_annual_84pc', 
                file_event_rates_84pc_mu_vary[event_uniform_event_row] * event_bias_adjustment_variable_mu)

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

        # For better load balancing, split up sources with row_weight>0 and row_weight=0
        unseg_row_weight_zero = which(sourcezone_parameters$segment_name == '' & sourcezone_parameters$row_weight==0)
        unseg_row_weight_nonzero = which(sourcezone_parameters$segment_name == '' & sourcezone_parameters$row_weight!=0)

        # FIXME: Consider replacing the following mclapply with clusterMap,
        # since the 'mc' functions have trouble shutting down workers on NCI

        # The following does the 'heavy computation'
        source_envs[unseg_row_weight_nonzero] = mclapply(
            as.list(1:length(source_segment_names))[unseg_row_weight_nonzero], parfun, 
            mc.cores=config$MC_CORES, mc.cleanup=9L)
        # The following does a quick, shortcut computation
        source_envs[unseg_row_weight_zero] = mclapply(
            as.list(1:length(source_segment_names))[unseg_row_weight_zero], parfun, 
            mc.cores=config$MC_CORES, mc.cleanup=9L)

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

    if(length(seg) > 0){


        # Split the heavy/light computations for better load balance
        seg_zero_row_weight = which(sourcezone_parameters$segment_name != '' & sourcezone_parameters$row_weight == 0)
        seg_nonzero_row_weight = which(sourcezone_parameters$segment_name != '' & sourcezone_parameters$row_weight != 0)

        # FIXME: Consider replacing the mcmapply calls below with clusterMap,
        # since the 'mc' functions have trouble shutting down workers on NCI

        # Heavy computations grouped together for load balance
        source_envs[seg_nonzero_row_weight] = mcmapply(parfun, 
            i = as.list(1:length(source_segment_names))[seg_nonzero_row_weight], 
            unsegmented_edge_rate_multiplier=unsegmented_edge_rate_multiplier[seg_nonzero_row_weight], 
            SIMPLIFY=FALSE,
            mc.cores=config$MC_CORES, mc.cleanup=9L)

        # Quick-exit computations
        source_envs[seg_zero_row_weight] = mcmapply(parfun, 
            i = as.list(1:length(source_segment_names))[seg_zero_row_weight], 
            unsegmented_edge_rate_multiplier=unsegmented_edge_rate_multiplier[seg_zero_row_weight], 
            SIMPLIFY=FALSE,
            mc.cores=config$MC_CORES, mc.cleanup=9L)
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
# On partially segmented sources we need to do some tricky things to the percentile
# calculations to make it consistent. The key reason is that if the 'full-source-zone'
# Mw-frequency curve comes from e.g. 50% unsegmented, and 50% co-monotonic segments,
# then if we are computing percentiles from the combined Mw-frequency curve, we should
# not just evaluate the percentiles at the unsegmented/segemented individually. 
#
unique_source_names = unique(sourcezone_parameters$sourcename)
# Currently only do this in serial. If it is a problem, we could restructure
# update_scenario_rate_per.... so that it does not directly modify the environment
# (because it is hard to get that to work in parallel). Instead it could return
# the required data to the main process, which would then update the environments itself
for(i in 1:length(unique_source_names)){
    inds = which(sourcezone_parameters$sourcename == unique_source_names[i])
    update_scenario_rate_percentiles_on_source_zones_with_partial_segmentation(source_envs[inds])
}

#
# OUTPUT RATES TO NETCDF FILES BELOW HERE
#
if(config$write_to_netcdf){

    # Zero rates in netcdf files by setting scale_rate to zero
    for(i in 1:length(source_segment_names)){
        # No need to repeat this for segments
        if(source_envs[[i]]$is_a_segment) next

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
        # Skip if the row weight is zero
        if(rate_scale == 0) next
        write_rates_to_event_table(source_envs[[i]], scale_rate = rate_scale, add_rate=TRUE)
    }
}

#
# Plot all the rate curves to a single pdf
#

xlim = c(7.0, 9.7)
ylim = c(1.0e-06, 10)
pdf('rate_curves_on_source_zones.pdf', width=9, height=7)
mw_global = seq(xlim[1]-dMw/2, xlim[2]+dMw/2, by=0.1)
global_exceedance_rate_mw_variable_mu = mw_global*0
global_exceedance_rate_mw_fixed_mu = mw_global*0
for(i in 1:length(source_segment_names)){
    
    # Skip ones with a row weight of zero
    rate_scale = as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight)
    if(rate_scale == 0) next
     
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
    event_rates_file_variable_mu = ncvar_get(fid, 'variable_mu_rate_annual')
    event_Mw_file_variable_mu = ncvar_get(fid, 'variable_mu_Mw')
    nc_close(fid)
    empirical_mean_curve = sapply(mw, 
        f<-function(x) sum(event_rates_file * (event_Mw_file >= x)))
    points(mw, empirical_mean_curve, pch=19, cex=0.2, col='pink')
    empirical_mean_curve_variable_mu = sapply(mw, 
        f<-function(x) sum(event_rates_file_variable_mu * (event_Mw_file_variable_mu >= x)))
    points(mw, empirical_mean_curve_variable_mu, pch=19, cex=0.2, col='blue')

    #
    # Compute global empirical rate IN THE FILES for a plot later.
    # To avoid double counting, we should not do this for 'segments',
    # just full sources (segments and full-sources both write to the same file)
    if(!source_envs[[i]]$is_a_segment){
        global_exceedance_rate_mw_variable_mu = global_exceedance_rate_mw_variable_mu + 
            sapply(mw_global, f<-function(x){
                sum(event_rates_file_variable_mu * (event_Mw_file_variable_mu >= x))})
        global_exceedance_rate_mw_fixed_mu = global_exceedance_rate_mw_fixed_mu + 
            sapply(mw_global, f<-function(x){
                sum(event_rates_file * (event_Mw_file >= x))})
    }

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
    abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), 
        col='orange', lty='dotted')

    # Posterior of slip, Mw_max, and b
    num_diff<-function(x, y){
        # Quick numerical derivative
        N = length(x)
        c( (y[2]-y[1])/(x[2]-x[1]), 
           (y[3:N] - y[1:(N-2)])/(x[3:N] - x[1:(N-2)]), 
           (y[N] - y[N-1])/(x[N] - x[N-1]) )
    }

    # Convenience plotting function
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

    # A few more plots showing the posterior distribution of parameters

    par(mfrow=c(2,2))
    
    plot_derivs('slip_rate')
    plot_derivs('Mw_max')
    plot_derivs('b')

    # Add a plot showing the weight that mw is possible
    mw_is_possible = source_envs[[i]]$mw_rate_function(mw_global, epistemic_nonzero_weight=TRUE)
    plot(mw_global, mw_is_possible, title='Weight that Mw is possible', ylim=c(0,1), t='h')
    grid()

}

# Globally integrated rates
mw = seq(MW_MIN, MAXIMUM_ALLOWED_MW_MAX, by=dMw/2) #all_rate_curves$Mw_seq
rate_vals = mw*0
rate_vals_variable_mu = mw*0
gcmt_global = data.frame()
for(i in 1:length(source_segment_names)){

    # Skip sources with row_weight=0
    if( as.numeric(source_envs[[i]]$sourcezone_parameters_row$row_weight) == 0) next

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

    # Add rates derived from files
    points(mw_global, global_exceedance_rate_mw_variable_mu, t='l', col='blue')
    points(mw_global, global_exceedance_rate_mw_fixed_mu, t='l', col='purple')

}
grid(col='orange')
abline(h=c(1,1/10, 1/100, 1/1000, 1/10000, 1/100000, 1/1000000), col='orange', 
    lty='dotted')

dev.off()





#
# Save the image and plot the implied convergence rates
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
    thrust_sources = c('sunda2', 'kermadectonga2', 'puysegur2', 'newhebrides2', 'timortrough', 'solomon2')
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
    # moment rates on outer-rise sources {which used the '0.45%' result in
    # Sleep (2012), and were clearly too low}
    thrust_sources = c('sunda2', 'kermadectonga2', 'puysegur2', 'newhebrides2', 'timortrough', 'solomon2')
    for(i in 1:length(thrust_sources)){
        uss = bird2003_env$unit_source_tables[[thrust_sources[i]]]
        div_vec = pmax(0, -uss$bird_vel_div)
        rl_vec = uss$bird_vel_rl
        deg2rad = pi/180
        allowed_rake_deviation_radians = config$rake_deviation * deg2rad
        # Restrict angle to +- rake_deviation of pure thrust (or pure normal)
        rl_vec = sign(rl_vec) * pmin(abs(rl_vec), div_vec * tan(allowed_rake_deviation_radians))
        convergent_slip = sqrt(rl_vec**2 + div_vec**2)
        slip_x_area = sum(uss$length*uss$width*1e+06*convergent_slip/1000)
        print(c(thrust_sources[i], slip_x_area))
    }
    outerrise_sources = c('outerrisesunda', 'outerrise_kermadectonga', 'outerrise_puysegur', 
        'outerrisenewhebrides', 'outer_rise_timor', 'outerrisesolomon')
    for(i in 1:length(outerrise_sources)){
        uss = bird2003_env$unit_source_tables[[outerrise_sources[i]]]
        div_vec = pmax(0, uss$bird_vel_div)
        rl_vec = uss$bird_vel_rl
        deg2rad = pi/180
        allowed_rake_deviation_radians = config$rake_deviation * deg2rad
        # Restrict angle to +- rake_deviation of pure thrust (or pure normal)
        rl_vec = sign(rl_vec) * pmin(abs(rl_vec), div_vec * tan(allowed_rake_deviation_radians))
        convergent_slip = sqrt(rl_vec**2 + div_vec**2)
        slip_x_area = sum(uss$length*uss$width*1e+06*convergent_slip/1000)
        print(c(outerrise_sources[i], slip_x_area))
    }

}
