#'
#' Given Mw, compute parameters which characterise the 'probability
#' distribution' of slip, supposing area varies as a log-normal distribution
#' according to the Strasser (2010) scaling relation.
#'
#' FIXME: This function is not currently exported -- do we want to keep it?
#' Seems a good idea. \cr
#' We assume mu is constant, and log10(Area) has a normal distribution. Since
#' the seismic moment = mu*slip*Area, this implies that log10(Slip) has a normal
#' distribution as well. So the mean and standard deviation of the log10
#' variable completely describe it. As well as returning those values, we return
#' the mean slip (mean of untransformed slip) since its calculation is slightly
#' inconvenient.
#'
#' @param Mw the moment magnitude
#' @param mu Shear Modulus (Pascals)
#' @param constant Value of constant passed to \code{M0_2_Mw}
#' @param relation scaling relation type, passed to \code{Mw_2_rupture_size}
#' @return A list containing 1) The log10(Slip) mean and stanard deviation and
#' 2) The associated mean slip (mean of the untransformed slip)
#'
compute_slip_density_parameters<-function(Mw, mu=3e+10, constant=9.05, relation='Strasser'){

    M0 = M0_2_Mw(Mw, inverse=TRUE, constant=constant)

    size_scaling = Mw_2_rupture_size(Mw, relation=relation, detailed=TRUE, CI_sd=1)

    # km
    mu_log10_area = log10(size_scaling$values['area'])
    sd_log10_area = size_scaling$log10_sigmas['area']

    # To multiply area by 1e+06 (i.e. km^2 --> m^2) we just add log10(1e+06) to
    # mu_log10_area
    mu_log10_area = mu_log10_area + log10(1.0e+06)

    mu_log10_slip = (-mu_log10_area + log10(M0) - log10(mu))
    sd_log10_slip = sd_log10_area

    log10_slip_par = c(mu_log10_slip, sd_log10_slip)
    names(log10_slip_par) = c('log10_mean', 'log10_sd')

    # If log10(slip) is normally distributed, then this is the mean_slip for
    # ruptures of size Mw:
    mean_slip = exp( mu_log10_slip*log(10) + 0.5*( (sd_log10_slip*log(10))**2))
    # [NOTE: To derive this formula, consider a variable y ~ exp(x) where x is
    # normal with given mean and sigma. We know mean(y) = exp(mean + sigma^2/2).
    # Now consider x' = log10(y). We know this has mean' and sigma' which are the
    # same as x **but divided by ln(10)**
    # We know mean(y) = exp(mean + sigma^2/s), but obviously this is the same as
    #                   exp([mean' x ln(10)] + [sigma' x ln(10)]^2 /2)
    # Hence if y = 10**(x'), then mean(y) can be computed in this way

    return(list(log10_slip_par = log10_slip_par, mean_slip = mean_slip))
}

#'
#' Compute the fraction of seismic moment associated with earthquakes having a
#' magnitude above some threshold
#'
#' For long-term slip rate computations we need to know how much slip is caused
#' by events with Mw > mwmin. This function helps compute that, though it
#' requires appropriately designed inputs
#'
#' @param Mws numeric vector with evenly spaced Mw values, increasing, densely
#' spaced enough for integration
#' @param rate_Mws numeric vector with an occurrence rate for each value in
#' Mws. THIS IS NOT AN EXCEEDANCE RATE, IT IS A RATE FOR EVENTS WITH Mw = Mws[i].
#' The value will thus depend on how closely spaced Mws is.
#' @param moment_Mws vector with the seismic moment for each Mws [follows from
#' standard formulae]
#' @param mwmin numeric We want to know the fraction of the moment associated
#' with events of this size or larger
#' @return number between 0-1, giving the fraction of moment associated with
#' events with Mw > mwmin
compute_moment_fraction_from_events_greater_or_equal_than_mwmin<-function(
    Mws, rate_Mws, moment_Mws, mwmin){

    dMw = Mws[2] - Mws[1]
    if( !isTRUE(all.equal(diff(Mws), rep(dMw,len=length(Mws)-1))) ){
        stop("Mws must be evenly spaced")
    }

    if(mwmin < Mws[1] | mwmin > max(Mws)){
        stop("mwmin not compatible with Mws")
    }

    events_high = which(Mws >= mwmin)

    fraction = sum( rate_Mws[events_high] * moment_Mws[events_high] ) / 
        sum(rate_Mws * moment_Mws)

    return(fraction)
}


#' Assign conditional probabilities to earthquakes with the same Mw
#'
#' Consider all the earthquake events on a discrete source with magnitude Mw.
#' We wish to assign conditional probabilities to these events (conditional on
#' the fact that an event of size Mw occurred). This allows us to compute the rate
#' of any individual event as: \cr
#' 'event_conditional_probability * rate_of_events_of_size_Mw' \cr
#' One way to do this is assume that for a given magnitude, all event
#' conditional probabilities are equal ('all_equal'). \cr
#' This might be reasonable if the events are identical in area and slip, and
#' just differ in location. \cr 
#' Another approach is to assume the conditional probabilities vary
#' inversely with the average slip of the event ('inverse_slip'). This might be
#' reasonable if the individual events are almost the same, but there are slight
#' differences in the event area and slip, and we want each event to make the
#' same contribution to the long-term average slip rate. \cr
#' Other methods could be proposed, and this routine also supports a
#' user-defined function
#'
#' @param event_table A data.frame containing information on all earthquake
#' events on the source-zone. It's the output of
#' \code{get_all_earthquake_events}
#' @param conditional_probability_model character or function. Supported
#' character values are 'all_equal' and 'inverse_slip', corresponding to cases
#' discussed above. \cr
#' If a function, it must take a single argument, which is a subset of the
#' event_table containing all events with fixed Mw, and return the conditional
#' probability of each event.
#' @return A numeric vector with length = nrow(event_table) giving the
#' conditional probabilities. 
#' @export
#'
get_event_probabilities_conditional_on_Mw<-function(
    event_table,
    conditional_probability_model = 'all_equal'){

    # Define functions corresponding to supported conditional_probability
    # models. 
    # These take as input a subset of the event table containing events with
    # identical Mw
    if(is.character(conditional_probability_model) && 
        (conditional_probability_model == 'all_equal')){
        # 'all_equal' probability function
        pfun<-function(event_table_Mw){
            N = length(event_table_Mw[,1])
            return( rep(1, N)/N)
        }

    }else if(is.character(conditional_probability_model) &&
        (conditional_probability_model == 'inverse_slip')){

        # 'inverse_slip' probability function
        pfun<-function(event_table_Mw){
            slip_inv = 1/event_table_Mw$slip
            return(slip_inv/sum(slip_inv))
        }

    }else if( class(conditional_probability_model) == 'function'){

        # Here the user must provide a function which can convert a subset of
        # the event table, containing all events with fixed Mw, to a
        # conditional probability

        pfun = conditional_probability_model

    }else{

        stop('conditional_probability_model value not recognized')

    }

    # Find unique Mw values in event table
    event_table_Mws = sort(unique(event_table$Mw))

    # Compute conditional_probabilities
    event_conditional_probabilities = rep(NA, length(event_table[,1]))
    for(Mw in event_table_Mws){
        ii = which(event_table$Mw == Mw)
        event_conditional_probabilities[ii] = pfun(event_table[ii,])

        # Ensure they sum to 1.0
        if(!isTRUE(all.equal(sum(event_conditional_probabilities[ii]), 1.0))){
            msg = paste0('ERROR: event conditional probabilities with Mw = ',
                Mw, ' do not seem to sum to 1.0. This suggests an error in ',
                'the conditional_probability_model')
            stop(msg)
        }
    }

    # Final check
    if(any(is.na(event_conditional_probabilities))){
        warnings('NA values among event_conditional_probabilities')
    }

    return(event_conditional_probabilities)
}

#'
#' Get a function to compute the rate of earthquakes with magnitude > Mw
#'
#' Create a function f(Mw) which returns the rate of events with 
#' magnitude >= Mw. Here Mw is a vector of earthquake magnitudes. \cr
#' The user provides parameters for a variant of the Gutenberg Richter
#' model, although the parameter 'a' is adjusted to match the long-term
#' seismogenic slip rate (m/year). Details are explained below. \cr
#' To account for uncertainties in the GR model parameters, the user may specify
#' each as a vector of one-or-more parameters with associated probabilities
#' (which must sum to 1, of course). The code then treats all combinations of
#' those parameters as having probability equal to the product of the individual
#' parameter probabilities. It then integrates over the chosen Gutenberg Richter
#' models to produce a single exceedance rate function, which is returned. \cr
#' The output rate function also supports a number of optional arguments, which
#' provide various details on the rate function on each branch of the logic
#' tree. See details below \cr
#'
#' @note The 'a' parameter in the GR model is constrained to match the long-term
#' seismogenic slip rate.  Given particular values of all the parameters except
#' a, and a long-term seisomogenic slip rate (= long_term_slip_rate *
#' coupling_coefficient), the long-term moment rate should equal the rate
#' integrated from the individual events. Mathematically: \cr
#' \cr
#' (definitions) \cr
#' \eqn{ \mu } = Shear modulus
#' \cr
#' S = long term semsmogenic slip rate 
#' \cr
#' A = source-zone area 
#' \cr
#' \eqn{ s_{i} }= slip for event i 
#' \cr
#' \eqn{ a_{i} } = area for event i
#' \cr
#' \eqn{rMw_{i} } = rate of events with the same magnitude Mw as event i
#' \cr
#' \eqn{ Pr(event_{i} | event with size Mw_{i}) } = the conditional probability
#' of each event with the same Mw in the event table [i.e. if we take all events
#' with a given Mw from the table, their conditional probabilities should sum to
#' 1]. 
#' \cr
#' (main equation) \cr
#' \deqn{ \mu  S  A = 
#' \sum_{i \in events} ( \mu s_{i} a_{i} Pr(event_{i} | event with size Mw_{i}) rMw_{i}) }\cr
#' \cr
#' The form of the GR relation implies that a factor of 10^(a) will appear in
#' the RHS term \eqn{rMw_{i}}, and nowhere else. Hence to
#' compute 'a', we just compute the LHS and RHS of the above equations,
#' tentatively assuming a=0 in the latter case, and then it follows that a =
#' log10(LHS/RHS). \cr
#' The code also optionally allows us to account for the fraction of seismic
#' moment released by earthquakes with magnitude smaller than Mw_min. This
#' reduces the LHS term in the equation above.
#'
#' @param slip_rate numeric vector with 1 or more long-term seismogenic slip
#' rates for the source-zone (m/year). This gives the long-term average slip on
#' the fault due to seismic events (creep is ignored, so these rates must
#' account for the rate of coupling)
#' @param slip_rate_prob numeric vector with weights for the slip_rate
#' parameters
#' @param b numeric vector with 1 or more 'b' parameters
#' @param b_prob numeric vector with weights for the b parameters
#' @param Mw_min numeric vector with 1 or more 'Mw_min' parameters
#' @param Mw_min_prob numeric vector with weights for the Mw_min parameters
#' @param Mw_max numeric vector with 1 or more 'Mw_max' parameters
#' @param Mw_max_prob numeric vector with weights for the Mw_max parameters
#' @param sourcezone_total_area The total area of the sourcezone (km^2)
#' @param event_table A data.frame containing information on all events in the
#' source-zone. Output of \code{get_all_earthquake_events}. This is needed to
#' evaluate the RHS of the seismic moment balance equation.
#' @param event_conditional_probabilities A vector of conditional probabilities
#' for all events in the source-zone (conditional on the fact than an event with
#' their magnitude did occur). Output of
#' \code{get_event_probabilities_conditional_on_Mw}. This is needed to evaluate
#' the RHS of the seismic moment balance equation.
#' @param computational_increment For each parameter combination in the logic
#' tree, the rate function is computed at points in a sequence from min(Mw_min)
#' to max(Mw_max), with spacing computational_increment (adjusted if required so
#' that the extremes are included). The function value is then evaluated via
#' lookup. 
#' @param Mw_frequency_distribution character vector giving the variant on the
#' Gutenberg Richter model used (if more than 1 they are passed to the logic
#' tree).  Either 'truncated_gutenberg_richter' or
#' 'characteristic_gutenberg_richter'. 
#' @param Mw_frequency_distribution_prob numeric vector. Logic tree weights for
#' each Mw_frequency_distribution type
#' @param update_logic_tree_weights_with_data logical. If TRUE, then a value for
#' Mw_count_duration must be provided. The weights for each parameter
#' combination in the tree are then updated with Bayes-theorem, based on the
#' probability that they would produce Mw_count_duration[2] events which exceed
#' Mw_count_duration[1] within time duration Mw_count_duration[3]. This assumes
#' earthquake timings follow a Poisson-process.
#' @param Mw_count_duration numeric vector of length 3 providing some data. The
#' first entry is an Mw value, the second is the number of events exceeding Mw
#' at the site, and the third is the duration of observation (in years).  If
#' \code{update_logic_tree_weights_with_data=TRUE} but Mw_obs_data$t is not
#' provided, this is used to re-weight parameter combinations in the logic
#' tree, with more weight given to those which predict similar event
#' frequencies as the data.
#' @param Mw_obs_data optional list containing earthquake event data
#' (corresponding to Mw_count_duration) which is used to re-weight logic tree
#' branches. This allows for more detailed use of data, as compared with the
#' simple Mw_count_duration approach. The use of data is only attempted if
#' "update_logic_tree_weights_with_data=TRUE". \cr
#' Mw_obs_data may contain: (1) a member Mw_obs_data$Mw having moment-magnitude
#' data for the historic events, and/or (2) a member Mw_obs_data$t, giving the
#' time in years of the event **since the observational start time**. \cr 
#' Either of the 'Mw' or 't' vectors vectors may also be NULL, in which case
#' they will be ignored (and in the 't' case, we will just use Mw_count_duration
#' if provided). \cr
#' If Mw_obs_data$Mw is not null, we use the likelihood of the 'Mw' data under
#' each logic tree branch to re-weight the logic tree branches [in addition to
#' use of temporal data]. The data must have "min(Mw_obs_data$Mw) >=
#' Mw_count_duration[1]" and "length(Mw_obs_data$Mw) = Mw_count_duration[2]".
#' The likelihood is evaluated by numerically differentiating each modelled GR
#' curve (normalised to a density for values above Mw_count_duration[1]). The
#' numerical differentiation uses a central difference over eps=1.0e-04
#' magnitude units, which should be fairly accurate although not perfectly
#' exact. \cr
#' If Mw_obs_data$t is not null (NOT RECOMMENDED DUE TO BIAS, see below), then
#' we assume event occurrence times are a realisation of a Poisson process
#' (treating the time before the first event and the time after the last event
#' as censored observations). BEWARE THAT CENSORED MAXIMUM LIKELIHOOD IS BIASED
#' DOWNWARD FOR THIS PROBLEM, SEE
#' https://stats.stackexchange.com/questions/133347/ml-estimate-of-exponential-distribution-with-censored-data
#' In practice it is better to use the 'Poisson counts' approach above, but we
#' include the censored approach because it will generalise to non-Poisson interevent
#' times more easily (not currently implemented). With the censored approach,
#' note the time between the last event and the end of observations is inferred
#' as (Mw_count_duration[3] - t[length(t)]), while the time before the first
#' event is just t[1], so we must have Mw_obs_data$t[1] >= 0 and
#' max(Mw_obs_data$t) <= Mw_count_duration[3]. \cr
#' It is worth noting a situation when the use of detailed temporal data is identical to the
#' simpler approach of providing only Mw_count_duration. Consider a case where
#' Mw_obs_data$t=NULL and Mw_obs_data$Mw=NULL, and Mw_count_duration = c(7.5, 3,
#' 50). This will give the same answer as having Mw_obs_data$Mw=NULL,
#' Mw_obs_data$t=c(0, 10, 30, 50), Mw_count_duration = c(7.5, 4, 50). Notice how
#' A) the observed event times begin at 0 and end at the last observed event
#' time; B) We need one more event than when we just provided the
#' Mw_count_duration data! This is because the 'censored' information provides
#' no information in this case [we already know 'time-between-events' >=0], and
#' we end up with 3 useful 'time-between-event' values.
#' @param Mw_2_M0 function which takes an earthquake magnitude and returns the
#' corresponding seismic moment M0
#' @param account_for_moment_below_mwmin logical. If FALSE assume all
#' seismically coupled slip occurs due to earthquakes with Mw > Mwmin. If TRUE,
#' estimate the fraction of seismically coupled slip that occurs due to
#' earthquakes with Mw > Mwmin, and use that when computing the rates.
#' @return Function f(Mw) of a numeric vector Mw, which by default returns the
#' rate of events with magnitude > Mw. It has a number of optional arguments.
#' \cr
#' If the optional argument 'bounds=TRUE' is passed,
#' then the function returns upper/lower/median rates of the logic tree as well
#' as the probability weighted mean rate. \cr
#' If the optional argument \code{return_all_logic_tree_branches=TRUE} is
#' passed, then the function returns a list containing a data.frame with all
#' combinations of parameters, a vector with the Gutenberg-Richter 'a' parameters
#' (that are derived from the latter in the current function), a vector with the respective
#' probabilities(weights), a vector with the Mw sequence at which the function
#' is tabulated, and a matrix with the tabulated values for every branch of the
#' logic tree. The latter option can be used to scrutinize the results from each
#' branch of the logic tree. \cr
#' If the optional argument \code{quantiles} = a vector of probabilities is
#' passed, then the function evaluates the rate function at each given quantile
#' of the logic tree (e.g. quantiles = c(0.1, 0.5, 0.9) would return the lower
#' 10\%, median, and upper 90\% rates from the logic tree at each Mw) \cr
#' If the optional argument \code{return_random_curve=TRUE} is passed, then the
#' function will randomly pick a parameter combination from the logic tree (with
#' probability equal to the probability of each branch), and evaluate the
#' exceedance rate of Mw for that parameter combination. This can be useful for
#' some monte-carlo computations. \cr
#' Only one of these optional arguments can be passed at a time.
#' @param mw_max_posterior_equals_mw_max_prior logical. If TRUE, then do the Bayesian
#' update on a 'per-mw-max' level, i.e. the posterior of mw_max will be the
#' same as the prior of mw_max.
#' @param mw_observation_error_cdf If not NULL, a function F(x, Mw_true) giving the
#' probability that {(Mw_obs - Mw_true) <= x} if Mw-true takes a known value. Here
#' the observed Mw_obs has some error. In other words, if Mw_true is fixed,
#' then F is a cumulative distribution function for the error in Mw_obs. While this
#' CDF does not have to be a function of Mw_true, it is useful to allow that.
#' @export
#'
rate_of_earthquakes_greater_than_Mw_function<-function(
    slip_rate, 
    slip_rate_prob, 
    b, 
    b_prob, 
    Mw_min, 
    Mw_min_prob, 
    Mw_max, 
    Mw_max_prob, 
    sourcezone_total_area, 
    event_table, 
    event_conditional_probabilities,
    computational_increment=0.01,
    Mw_frequency_distribution='truncated_gutenberg_richter',
    Mw_frequency_distribution_prob = 1,
    update_logic_tree_weights_with_data = FALSE,
    Mw_count_duration = c(NA, NA, NA),
    Mw_obs_data = list(Mw=NULL, t=NULL),
    Mw_2_M0 = function(x) M0_2_Mw(x, inverse=TRUE),
    account_for_moment_below_mwmin=FALSE,
    mw_max_posterior_equals_mw_max_prior = FALSE,
    mw_observation_error_cdf = NULL){

    # Check data
    stopifnot(length(slip_rate) == length(slip_rate_prob))
    stopifnot(length(b) == length(b_prob))
    stopifnot(length(Mw_min) == length(Mw_min_prob))
    stopifnot(length(Mw_max) == length(Mw_max_prob))
    stopifnot(length(Mw_frequency_distribution) == 
        length(Mw_frequency_distribution_prob))

    # Check (conservative) on the input Mw_min/Mw_max
    stopifnot(min(Mw_max) > max(Mw_min))

    # Check probabilities sums to 1 (within floating point roundoff)
    stopifnot(isTRUE(all.equal(sum(slip_rate_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(b_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_min_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_max_prob), 1)))
    stopifnot(isTRUE(all.equal(sum(Mw_frequency_distribution_prob), 1)))

    # # Check that the range of Mw in the table covers the range of
    # # Mw_min/Mw_max)
    # stopifnot(min(event_table$Mw) <= min(Mw_min))
    # stopifnot(max(event_table$Mw) >= max(Mw_max))


    # Get all combinations of the parameters, and all combinations of the
    # probabilities
    all_par_combo = expand.grid(
        slip_rate = slip_rate, 
        b = b, 
        Mw_min = Mw_min, 
        Mw_max = Mw_max, 
        Mw_frequency_distribution=Mw_frequency_distribution)
    all_par_combo_prob = expand.grid(
        slip_rate_prob = slip_rate_prob,
        b_prob = b_prob, 
        Mw_min_prob = Mw_min_prob, 
        Mw_max_prob = Mw_max_prob,
        Mw_frequency_distribution_prob = Mw_frequency_distribution_prob)

    # Collapse the probabilities to a single one for each set of table
    # parameters
    all_par_prob = rep(NA, nrow(all_par_combo))
    for(i in 1:nrow(all_par_combo)){
        all_par_prob[i] = prod(all_par_combo_prob[i,])    
    }
    rm(all_par_combo_prob)
  
    # Check the summed probabilities still = 1
    stopifnot(isTRUE(all.equal(sum(all_par_prob), 1.0)))

    # Make sequence at which function evaluations occur
    min_Mw_min = min(Mw_min)
    max_Mw_max = max(Mw_max)
    npts = round((max_Mw_max - min_Mw_min)/computational_increment) + 1
    Mw_seq = seq(min_Mw_min, max_Mw_max, len=npts)


    # Get all Mw values in the eq table, and find the spacing between unique
    # values
    table_Mw_values = sort(unique(event_table$Mw))
    table_Mw_increment = diff(table_Mw_values)
    range_of_table_increments = diff(range(table_Mw_increment))
    
    if(!isTRUE(all.equal(range_of_table_increments, 0.))){
        print('Earthquake table has these Mw values:')
        print(table_Mw_values)
        stop('The sorted unique Mw values must be an evenly-spaced sequence')
    }
    table_Mw_increment = table_Mw_increment[1]

    # Shorthand notation
    eq_Mw = event_table$Mw
    eq_slip = event_table$slip
    eq_area = event_table$area

    # Check that all event_conditional_probabilities
    # sum to 1.0, for fixed Mw
    for(Mw in table_Mw_values){

        ii = which(eq_Mw == Mw)

        conditional_prob_sum = 
            sum(event_conditional_probabilities[ii])

        if(!isTRUE(all.equal(conditional_prob_sum, 1.0))){
            print(paste0('Error: conditional probabilities for Mw = ', Mw, 
                ' do not sum to 1.0'))
            print(paste0('They sum to: ', conditional_prob_sum))
            print('The values are:')
            print(event_conditional_probabilities[ii])
            stop('Invalid conditional probabilities')
        }

    }
    
    # Compute the rates for each final branch of the logic tree
    all_rate_matrix = matrix(NA, ncol = length(Mw_seq), 
        nrow=nrow(all_par_combo))
    a_parameter = rep(NA, nrow(all_par_combo)) # GR 'a' parameter
    for(i in 1:nrow(all_rate_matrix)){

        # Get parameter vector
        par = all_par_combo[i,] 
        
        if(par$Mw_frequency_distribution == 'truncated_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_truncated_gutenberg_richter
        }else if(par$Mw_frequency_distribution == 'characteristic_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_characteristic_gutenberg_richter
        }else{
            stop(paste0(" Value of Mw_frequency_distribution: ", 
                par$Mw_frequency_distribution, " not recognized"))
        }

        # Ranges of Mw in our event table, with additional dx/2 spacing.
        # These are useful since the rate assigned to events with Mw = X is
        # GR(rate-dx/2) - GR(rate+dx/2)
        lower_Mw = eq_Mw - table_Mw_increment/2
        upper_Mw = eq_Mw + table_Mw_increment/2

        if(account_for_moment_below_mwmin){
            # Compute the fraction of seismic moment associated with
            # earthquakes with Mw ranging over the values in our event table

            # Finely spaced sequence used for numerical integration
            lower = min(par$Mw_min, min(lower_Mw), 0)  - 0.001
            upper = max(par$Mw_max, max(upper_Mw)) + 0.001
            broad_Mw_seq = seq(lower, upper, by=0.001)
            dmw0 = broad_Mw_seq[2] - broad_Mw_seq[1]

            # Get the rate of events with Mw in the finely spaced sequence. 
            broad_rates = Mfd(broad_Mw_seq - dmw0/2, a = 0, b = par$b, 
                              Mw_min=broad_Mw_seq[1], Mw_max = par$Mw_max) -
                          Mfd(broad_Mw_seq + dmw0/2, a = 0, b = par$b, 
                              Mw_min=broad_Mw_seq[1], Mw_max = par$Mw_max) 
            broad_seismic_moment = Mw_2_M0(broad_Mw_seq)

            # Get the fraction of seismic moment between min mw and max mw in
            # our earthquake events table
            # Step 1: above min(eq_Mw)
            moment_fraction = 
                compute_moment_fraction_from_events_greater_or_equal_than_mwmin(
                    broad_Mw_seq, broad_rates, broad_seismic_moment, min(lower_Mw))
            # Step 2: subtract off fraction above max(eq_Mw)
            moment_fraction = moment_fraction - 
                compute_moment_fraction_from_events_greater_or_equal_than_mwmin(
                    broad_Mw_seq, broad_rates, broad_seismic_moment, max(upper_Mw))
           
            stopifnot(moment_fraction >= 0 & moment_fraction <= 1)
        }else{
            moment_fraction = 1
        }

        # Evaluate LHS of seismic moment balance equation
        LHS = (sourcezone_total_area * 1e+06) * par$slip_rate * moment_fraction

        # We need to solve for 'a', given the long-term slip rate.
        # This can be done by initially setting 'a' to zero, then
        # back-calculating the value required to match the long-term slip
        rate_of_events_of_size_Mw_aeq0 = 
            Mfd(lower_Mw, a = 0, b = par$b, Mw_min=par$Mw_min, 
                Mw_max = par$Mw_max) -
            Mfd(upper_Mw, a = 0, b = par$b, Mw_min=par$Mw_min, 
                Mw_max = par$Mw_max)

        RHS = sum(eq_slip * (eq_area*1e+06) * rate_of_events_of_size_Mw_aeq0 * 
            event_conditional_probabilities)

        # Back-calculate 'a'
        # 10^a = LHS/RHS
        # NOTE: This assumes that 'a' appears in Mfd as 10^(a) {only}
        a_parameter[i] = log10(LHS/RHS)

        # Compute rates of exceedance for a truncated Gutenberg Richter model.
        all_rate_vec = Mfd(Mw_seq, a = a_parameter[i], b = par$b, 
            Mw_min = par$Mw_min, Mw_max = par$Mw_max)

        all_rate_matrix[i,] = all_rate_vec

    }
    # Store 'a' in all par combo (allows access from outside the function)
    all_par_combo$a = a_parameter

    # Save original probabilities in case we update them
    all_par_prob_prior = all_par_prob   
 
    if(update_logic_tree_weights_with_data){
        # Update the logic tree weights using provided data sets
        all_par_prob = compute_updated_logic_tree_weights(Mw_seq,
            all_par_combo, all_par_prob_prior, Mw_count_duration,
            Mw_obs_data, mw_max_posterior_equals_mw_max_prior)

        if(is.null(mw_observation_error_cdf)){
            # We did not provide an error model for magnitude observations, so do not
            # change them
            all_par_prob_with_Mw_error = all_par_prob
        }else{
            # Update the logic tree weights using provided data sets
            all_par_prob_with_Mw_error = compute_updated_logic_tree_weights(
                Mw_seq, all_par_combo, all_par_prob_prior, Mw_count_duration,
                Mw_obs_data, mw_max_posterior_equals_mw_max_prior,
                cdf_mw_observation_error=mw_observation_error_cdf)
        }

    }else{
        # Posterior = Prior
        all_par_prob = all_par_prob_prior
        all_par_prob_with_Mw_error = all_par_prob
    }
    gc()

    # Compute the average rate
    final_rates = rep(0, length(Mw_seq))
    final_rates_with_Mw_error = rep(0, length(Mw_seq))
    for(i in 1:nrow(all_rate_matrix)){
        final_rates = final_rates + all_par_prob[i] * all_rate_matrix[i,]
        final_rates_with_Mw_error = final_rates_with_Mw_error + 
            all_par_prob_with_Mw_error[i] * all_rate_matrix[i,]
    }

    # Store the max/min rates from all branches of the logic tree,
    # to help quantify uncertainty later.
    # Consider producing a fuller picture from all_rate_matrix
    upper_rate = apply(all_rate_matrix, 2, max)
    lower_rate = apply(all_rate_matrix, 2, min)

    # Function to compute quantiles (describing epistemic uncertainty)
    # of the logic tree curves at different exceedance rates
    quantile_rate_fun<-function(rates, p, with_mw_error=FALSE){
        sorted_rates = sort(rates, index.return=TRUE)
        if(with_mw_error){
            sorted_prob = all_par_prob_with_Mw_error[sorted_rates$ix]
        }else{
            sorted_prob = all_par_prob[sorted_rates$ix]
        }
        cumulative_sorted_prob = cumsum(sorted_prob)
        # The median occurs when the cumulative probability >= 0.5
        ind = min(which(cumulative_sorted_prob >= p))
        return(sorted_rates$x[ind])
    }


    # Compute the 50th percentile
    median_rate = apply(all_rate_matrix, 2, 
        f<-function(rates) quantile_rate_fun(rates, p=0.5)
        )
    median_rate_with_Mw_error = apply(all_rate_matrix, 2, 
        f<-function(rates) quantile_rate_fun(rates, p=0.5, with_mw_error=TRUE)
        )


    # For the output function, if we exceed the bounds on the LHS we can return
    # the most extreme value, and on the RHS we should return zero.
    # We take care of the LHS in approxfun, and the RHS below.
    output_function = approxfun(Mw_seq, final_rates, rule=2)
    output_function_with_Mw_error = approxfun(Mw_seq, final_rates_with_Mw_error, rule=2)
    upper_function = approxfun(Mw_seq, upper_rate, rule=2)
    lower_function = approxfun(Mw_seq, lower_rate, rule=2)
    median_function = approxfun(Mw_seq, median_rate, rule=2)
    median_function_with_Mw_error = approxfun(Mw_seq, median_rate_with_Mw_error, rule=2)

    # Function to compute the rate, with a number of useful extra options
    #
    # @param Mw the Mw value at which the rate is desired
    # @param bounds logical. If TRUE, return some summary statistics on the
    # range of rates suggested by the logic tree. More detailed information
    # can be obtained using quantiles or return_all_logic_tree_branches.
    # @param return_all_logic_tree_branches logical. If TRUE, return a list
    # with the all_rate_matrix, Mw_seq, all_par_prob, all_par_prob_prior,
    # and all_par. This option is useful for doing non-standard manipulations
    # of the results, and debugging.
    # @param quantiles numeric vector of probabilities or NULL. If not NULL,
    # then return quantiles of the exceedance rate for Mw (based on the
    # logic tree branches and probabilities)
    # @param return_random_curve TRUE/FALSE. If TRUE, randomly pick a parameter
    # combination from the logic tree (with probability equal to the probability
    # of each branch), and evaluate the exceedance rate of Mw for that
    # parameter combination. This can be useful for some monte-carlo
    # computations.
    # @param account_for_mw_obs_error Use curve weights that account for errors
    # in the observations that were used to derive the curve weights
    # @return Depends on the input arguments, see comments above.
    output_function2<-function(
        Mw, 
        bounds=FALSE, 
        return_all_logic_tree_branches=FALSE, 
        quantiles=NULL,
        return_random_curve=FALSE,
        account_for_mw_obs_error=FALSE){

        # Check input arguments for accidental nulls
        input_args = as.list(match.call())

        # If we actually provided the 'quantiles' argument, then it's unlikely
        # we want a NULL value, however, a typo could create this. Check here
        # for that situation.
        if('quantiles' %in% names(input_args)){
            if(is.null(quantiles)){
                print(paste0('Deliberately passed null value to the optional argument "quantiles"',
                        ' -- halting as this suggests a user typo'))
                stop()
            }
        }


        # Check that at most one 'special option' is taken
        quantiles_not_null = !is.null(quantiles)
        options_sum = bounds + return_all_logic_tree_branches + 
            quantiles_not_null + return_random_curve

        if(options_sum > 1){
            msg = paste0('At most one of bounds, ', 
                'return_all_logic_tree_branches, ',
                '!is.null(quantiles), and return_random_curve can be TRUE')
            stop(msg)
        }


        # Option to just dump the key data. This can be useful for
        # debugging / plotting non-standard outputs, etc 
        if(return_all_logic_tree_branches){
            output = list(all_rate_matrix = all_rate_matrix, 
                          Mw_seq = Mw_seq,
                          all_par_prob = all_par_prob,
                          all_par_prob_prior = all_par_prob_prior,
                          all_par = all_par_combo,
                          a_parameter = a_parameter,
                          all_par_prob_with_Mw_error=all_par_prob_with_Mw_error)
            return(output)
        }

        if(return_random_curve){
            # Grab a curve 'at random' with probabilities based on all_par_prob

            stopifnot(bounds==FALSE & return_all_logic_tree_branches==FALSE)

            if(account_for_mw_obs_error){
                par = sample(1:length(all_par_combo[,1]), size=1, 
                    prob=all_par_prob_with_Mw_error)
            }else{
                par = sample(1:length(all_par_combo[,1]), size=1, 
                    prob=all_par_prob)
            }
        
            output = approx(Mw_seq, all_rate_matrix[par,], xout=Mw, rule=2)$y *
                (Mw <= max_Mw_max)
            return(output)
        }

        # Typical case
        if(account_for_mw_obs_error){
            output = (Mw <= (max_Mw_max))*output_function_with_Mw_error(Mw)
        }else{
            output = (Mw <= (max_Mw_max))*output_function(Mw)
        }

        if(bounds){
            # Compute some summary statistics
            output_up = (Mw <= max_Mw_max)*upper_function(Mw)
            output_low = (Mw <= max_Mw_max)*lower_function(Mw)
            if(account_for_mw_obs_error){
                output_med = (Mw <= max_Mw_max)*median_function_with_Mw_error(Mw)
            }else{
                output_med = (Mw <= max_Mw_max)*median_function(Mw)
            }
           
            if(length(Mw) == 1){ 
                output = c(output, output_low, output_up, output_med)
                names(output) = c('rate', 'min', 'max', 'median')
            }else{
                output = cbind(output, output_low, output_up, output_med)
                colnames(output) = c('rate', 'min', 'max', 'median')
            }

        }else if(!is.null(quantiles)){
            # Compute results for a set of quantiles

            if(length(Mw) == 1){
                output = rep(NA, length(quantiles))
                names(output) = quantiles
            }else{
                output = matrix(NA, ncol=length(quantiles), nrow=length(Mw))
                colnames(output) = quantiles
            }
            for(i in 1:length(quantiles)){
                # For each quantile, evaluate the rate curve with interpolation
                quantile_rate = apply(all_rate_matrix, 2, 
                    f<-function(rate) quantile_rate_fun(rate, p=quantiles[i], 
                        account_for_mw_obs_error))
                if(length(Mw) == 1){
                    output[i] = 
                        approx(Mw_seq, quantile_rate, xout=Mw, rule=2)$y * 
                        (Mw <= max_Mw_max)
                }else{
                    output[,i] = 
                        approx(Mw_seq, quantile_rate, xout=Mw, rule=2)$y * 
                        (Mw <= max_Mw_max)
                }
            }
        }

        return(output)
    }

    return(output_function2)

}

#' Evaluate the rate of earthquakes with magnitude >= Mw using a doubly
#' truncated Gutenberg Richter distribution.
#'
#' The Gutenberg Richter distribution gives the number of earthquakes with
#' magnitude > Mw as:\cr
#' N_{GR}(x >= Mw) = 10^(a-bMw) \cr
#' where a and b are constants. Note 10^(a) is the rate of earthquakes with Mw >
#' 0. \cr
#' By differentiating N_{GR} we can estimate the number of events with 
#' (Mw - dx/2 <= magnitude <= Mw +dx/2) as:\cr
#' N(Mw -dx/2 <= x <= Mw + dx/2) ~= dx * n(Mw) = dx * ([10^(a-bMw)] * bln(10)) \cr
#' where n(Mw) is the negative of the derivative of N(x>Mw). \cr
#' n(Mw) is like a scaled probability-density-function of the number of
#' earthquakes. Unlike a pdf, n(Mw) integrates to the rate of earthquakes,
#' instead of 1.\cr
#' For the truncated Gutenberg Richter distribution, n(Mw) is truncated between
#' lower and upper Mw limits (i.e. set to zero outside these limits). \cr
#' We then have the equivalent of N_{GR} for the truncated distribution as: \cr
#' N_{TGR}(x >= Mw) = 10^(-max(Mw, Mw_min)*b + a) - 10^(-Mw_max*b + a) \cr
#' Notice that now, 10^(a) is no-longer the rate of earthquakes with Mw >= 0 --
#' instead that rate is :\cr
#' N_{TGR}(x >= Mw_min) = 10^(a - b*Mw_min) - 10^(a - b*Mw_max) \cr
#' (because there are no events with magnitude <= Mw_min)
#'
#' @param Mw Moment magnitude
#' @param a The a parameter
#' @param b The b parameter
#' @param Mw_min the lower truncated moment magnitude
#' @param Mw_max the upper truncated moment magnitude
#' @return The rate of events with magnitude > Mw
#' @export
Mw_exceedance_rate_truncated_gutenberg_richter<-function(
    Mw, a, b, Mw_min, Mw_max){
   
    N_mag_gt_Mw = 10^(a - b*pmax(Mw, Mw_min)) - 10^(a - b*Mw_max)

    output = N_mag_gt_Mw*(Mw <= Mw_max)
    
    return(output)
}



#' Evaluate the rate of earthquakes with magnitude > Mw using a characteristic
#' variant of the Gutenberg Richter distribution.
#' 
#' The 'pure' Gutenberg Richter distribution gives the number of earthquakes
#' with magnitude >= Mw as:\cr
#' N_{GR}(x >= Mw) = 10^(a-bMw) \cr
#' where a and b are constants. Note 10^(a) is the rate of earthquakes with Mw
#' >= 0. \cr
#' The characteristic formulation used here ajusts this as: \cr
#' N_{GR}(x >= Mw) = 10^(a-bMw) if Mw_min <= Mw <= Mw_max \cr
#' ~~~~~~~~~~~~~  = 10^(a-bMw_min) if Mw < Mw_min \cr
#' ~~~~~~~~~~~~~  = 0 if Mw > Mw_max \cr
#' This means that the rate of earthquakes of size EXACTLY Mw_max is finite.
#' @param Mw Moment magnitude
#' @param a The a parameter
#' @param b The b parameter
#' @param Mw_min the lower truncated moment magnitude
#' @param Mw_max the upper truncated moment magnitude
#' @return The rate of events with magnitude > Mw
#' @export
#' 
Mw_exceedance_rate_characteristic_gutenberg_richter<-function(
    Mw, a, b, Mw_min, Mw_max){

    N_mag_gt_Mw = 10^(a - b*pmax(Mw, Mw_min))

    output = N_mag_gt_Mw*(Mw <= Mw_max)
    return(output)
}

# Update logic tree weights using data
#
# This is a computational-workhorse routine used in
# rate_of_earthquakes_greater_than_Mw_function
#
# 
compute_updated_logic_tree_weights<-function(Mw_seq, all_par_combo, 
    all_par_prob_prior, Mw_count_duration, Mw_obs_data, 
    mw_max_posterior_equals_mw_max_prior, cdf_mw_observation_error=NULL){

    # These numbers are repeatedly used when we treat magnitude errors. 
    # We assume finite support in the mw error, with a max absolute value
    max_obs_mw_error = 0.5
    # and do numerical integration with a given 'dy' increment.
    integration_dy = 0.01

    #
    # Adjust the weights of logic-tree branches based on the data
    #
    if(sum(is.na(Mw_count_duration)) > 0){
        msg = paste0('Must provide Mw_count_duration data if ', 
            'update_logic_tree_weights_with_data=TRUE')
        stop(msg)
    }

    data_thresh = Mw_count_duration[1]
    data_count = Mw_count_duration[2]
    data_observation_duration = Mw_count_duration[3]
    
    if( (data_thresh < min(Mw_seq)) | (data_thresh > max(Mw_seq)) ){
        stop('The Mw threshold of the data must be within [Mw_min, Mw_max]')
    }
   
    model_rates = all_par_prob_prior*0
    for(i in 1:length(model_rates)){

        #
        # Compute the probability of observing the data, given each model
        #

        # GR parameters for this logic tree curve
        a_par = all_par_combo$a[i]
        b_par = all_par_combo$b[i]
        mw_min_par = -Inf # So we can do numerical derivatives
        mw_max_par = all_par_combo$Mw_max[i]
        Mfd_type = all_par_combo$Mw_frequency_distribution[i]

        if(Mfd_type == 'truncated_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_truncated_gutenberg_richter
        }else if(Mfd_type == 'characteristic_gutenberg_richter'){
            Mfd = Mw_exceedance_rate_characteristic_gutenberg_richter
        }else{
            stop(paste0(" Value of Mw_frequency_distribution: ",
                Mfd_type, " not recognized"))
        }


        if(is.null(cdf_mw_observation_error)){
            #model_rates[i] = approx(Mw_seq, all_rate_matrix[i,], 
            #    xout=data_thresh)$y
            model_rates[i] = Mfd(data_thresh, a = a_par, b=b_par, 
                Mw_min=mw_min_par, Mw_max=mw_max_par)
        }else{
            #
            # EXTENSION TO TREAT DATA ERRORS, OR MU VARIATION
            #
            # Below for the integration, we assume the absolute Mw errors are
            # always < 0.5. Good to add some check to ensure that is the case
            #
            tmp = cdf_mw_observation_error(-max_obs_mw_error , 
                data_thresh+max_obs_mw_error)
            if(tmp > 0){
                stop('mw_observation_error is too large')
            }
            tmp = cdf_mw_observation_error( max_obs_mw_error , 
                data_thresh-max_obs_mw_error)
            if(tmp < 1){
                stop('mw_observation_error is too large')
            }
            #
            Mfd_local<-function(y){
                Mfd(y, a=a_par, b=b_par, Mw_min=mw_min_par, Mw_max=mw_max_par)
            }

            int_range = data_thresh + c(-1,1)*max_obs_mw_error 
            model_rates[i] = exceedance_rate_of_observed(Mfd_local, 
                cdf_mw_observation_error, data_thresh,
                true_value_range_with_nontrivial_cdf_value = int_range, 
                integration_dy=integration_dy)
        }
    }

    #
    # Temporal component of likelihood
    #
    if(is.null(Mw_obs_data$t)){
        # Here, detailed temporal data was not provided. 
        # Assume poisson rate of events.

        log_pr_data_given_model = dpois(data_count, 
            lambda = (model_rates*data_observation_duration), 
            log=TRUE)

    }else{
        msg = paste0('Warning: Censored ML gives a biased estimate of the rates for \n',
            'a Poisson process, it is likely better to set Mw_obs_data$t=NULL, \n',
            'while retaining the Mw_obs_data$Mw values')
        print(msg)
        # Detailed temporal data was provided.
        # Times must be sorted
        stopifnot(min(diff(Mw_obs_data$t)) >= 0)
        stopifnot(min(Mw_obs_data$t) >= 0)
        stopifnot(max(Mw_obs_data$t) <= data_observation_duration)

        # Compute times between events
        dts = diff(Mw_obs_data$t)
        dts1_lower_bound = Mw_obs_data$t[1]
        dts_last_lower_bound = data_observation_duration - Mw_obs_data$t[data_count]
        
        log_pr_data_given_model = rep(NA, length=length(model_rates))
        for(i in 1:length(model_rates)){
            ri = model_rates[i] 
            # Likelihood function, exponential model. Account for fact that the 
            # first/last time spacings are not known exactly
            log_pr_data_given_model[i] = 
                # Sum log-likelihood for numerical stability
                sum(dexp(dts, rate=ri, log=TRUE)) + 
                # Integral over upper tails for first/last data points
                pexp(dts1_lower_bound, rate=ri, lower.tail=FALSE, log.p=TRUE) +
                pexp(dts_last_lower_bound, rate=ri, lower.tail=FALSE, log.p=TRUE)
        }

    }

    #
    # Magnitude component of likelihood
    #
    if(!is.null(Mw_obs_data$Mw)){

        # Check consistency between this data and Mw_count_duration
        stopifnot(length(Mw_obs_data$Mw) == data_count)
        stopifnot(all(Mw_obs_data$Mw >= data_thresh))

        # If times were provided, they must be the same length
        if(!is.null(Mw_obs_data$t)){
            stopifnot(length(Mw_obs_data$Mw) == length(Mw_obs_data$t))
        }

        # Update the likelihood
        for(i in 1:length(model_rates)){

            data_thresh = Mw_count_duration[1]
            # GR parameters for this logic tree curve
            a_par = all_par_combo$a[i]
            b_par = all_par_combo$b[i]
            mw_min_par = all_par_combo$Mw_min[i]
            mw_max_par = all_par_combo$Mw_max[i]
            Mfd_type = all_par_combo$Mw_frequency_distribution[i]

            if(Mfd_type == 'truncated_gutenberg_richter'){
                Mfd = Mw_exceedance_rate_truncated_gutenberg_richter
            }else if(Mfd_type == 'characteristic_gutenberg_richter'){
                Mfd = Mw_exceedance_rate_characteristic_gutenberg_richter
            }else{
                stop(paste0(" Value of Mw_frequency_distribution: ",
                    Mfd_type, " not recognized"))
            }       
             

            #
            # Evaluate the 'density' for samples above Mw_observation_threshold 
            #   = -(1/GR(data_thresh))*[ derivative_of_GR_with_respect_to_Mw]
            # (Because the CDF is ( 1  -(1/GR(data_thresh))*GR(Mw) )
            #
            if(is.null(cdf_mw_observation_error)){
                #
                # Ignoring observation errors in Mw data
                #
                gr_mwmin = Mfd(data_thresh, a = a_par, b=b_par, 
                    Mw_min=mw_min_par, Mw_max=mw_max_par)
                eps = 1e-04 # For numerical differentiation
                density_above_data_thresh = -1/(gr_mwmin*2*eps) * (
                    Mfd(Mw_obs_data$Mw+eps, a = a_par, b=b_par, Mw_min=mw_min_par, Mw_max=mw_max_par) - 
                    Mfd(Mw_obs_data$Mw-eps, a = a_par, b=b_par, Mw_min=mw_min_par, Mw_max=mw_max_par) )
            }else{
                #
                # Treating observation errors in Mw data
                #
                Mfd_local<-function(y){
                    Mfd(y, a=a_par, b=b_par, Mw_min=-Inf, Mw_max=mw_max_par)
                }

                int_range = data_thresh + c(-1,1)*max_obs_mw_error
                gr_mwmin = exceedance_rate_of_observed(Mfd_local,
                    cdf_mw_observation_error, data_thresh, 
                    true_value_range_with_nontrivial_cdf_value=int_range,
                    integration_dy = integration_dy)

                eps = 1e-04 # For numerical differentiation
                # Compute the density in a few steps:
                # density = -1/(gr_mwmin*2*eps) ( rate_plus_eps - rate_minus_eps)
                rate_plus_eps = sapply(Mw_obs_data$Mw + eps, f<-function(x){
                    exceedance_rate_of_observed(Mfd_local,
                        cdf_mw_observation_error, x, 
                        true_value_range_with_nontrivial_cdf_value=(x+c(-1, 1)*max_obs_mw_error),
                        integration_dy = integration_dy)
                })
                rate_minus_eps = sapply(Mw_obs_data$Mw - eps, f<-function(x){
                    exceedance_rate_of_observed(Mfd_local,
                        cdf_mw_observation_error, x, 
                        true_value_range_with_nontrivial_cdf_value=(x+c(-1, 1)*max_obs_mw_error),
                        integration_dy = integration_dy)
                })

                density_above_data_thresh = -1/(2*eps*gr_mwmin)*
                    (rate_plus_eps - rate_minus_eps)
            }

            stopifnot(all(density_above_data_thresh >= 0))

            # Add the Mw-related log-likelihood terms to the full log likelihood
            log_dens_values = log(density_above_data_thresh)
            log_pr_data_given_model[i] = log_pr_data_given_model[i] + (sum(log_dens_values))
        }

    }

    # Check that some models have non-zero probability
    
    if( all(!is.finite(log_pr_data_given_model)) ){
        stop('Mw data is impossible under every model')
    }
 
    # For numerical stability, rescale log_pr_data_given_model 
    # This will not effect all_par_prob, but makes it numerically easier to compute
    mx = max(log_pr_data_given_model)

    # Bayes theorem, with likelihood scaled for numerical stability
    all_par_prob =  all_par_prob_prior * exp(log_pr_data_given_model-mx) / 
        sum(all_par_prob_prior * exp(log_pr_data_given_model-mx))


    # Optionally, ENFORCE mw-max-posterior to be mw-max prior.
    # A user might want to do this e.g. if they want a fixed uncertainty in Mw-max, 
    # irrespective of the data. Obviously this option is not a pure Bayesian update.
    if(mw_max_posterior_equals_mw_max_prior){
        if(!is.na(Mw_count_duration[1])){
            if( any((all_par_prob_prior > 0) & (all_par_combo$Mw_max < Mw_count_duration[1])) ){
                stop(paste0('Mw_max prior puts non-zero weight on value that is less ',
                    'than the observations. \n This means some prior values are impossible, \n', 
                    'and we cannot satisfy mw_max_posterior_equals_mw_max_prior'))
            }
        }
        # Change the all_par_prob to ensure that mw_max posterior is the same
        # as mw_max prior. This is like we did the weight update for each mw_max separately.
        mw_max_post_prob = aggregate(all_par_prob, list(all_par_combo$Mw_max), sum)
        #mw_max_prior_prob = aggregate(Mw_max_prob, list(Mw_max), sum)
        mw_max_prior_prob = aggregate(all_par_prob_prior, list(all_par_combo$Mw_max), sum)

        stopifnot(isTRUE(all.equal(sum(mw_max_post_prob[,2]), 1.0)))
        stopifnot(isTRUE(all.equal(sum(mw_max_prior_prob[,2]), 1.0)))

        match_mw_max1 = match(all_par_combo$Mw_max, mw_max_post_prob[,1])
        match_mw_max2 = match(all_par_combo$Mw_max, mw_max_prior_prob[,1])

        k = which(all_par_prob > 0) # Avoid possibility of division by zero if priors included zero
        all_par_prob[k] = all_par_prob[k]/mw_max_post_prob[match_mw_max1[k],2] * mw_max_prior_prob[match_mw_max2[k],2]
        # Check that we equalized over Mw_max 
        mw_max_post_prob = aggregate(all_par_prob, list(all_par_combo$Mw_max), sum)
        stopifnot(isTRUE(all.equal(mw_max_post_prob[,2], mw_max_prior_prob[,2])))
    }

    return(all_par_prob)
}

#' Exceedance rate of observations with random errors
#'
#' Suppose that a process generates values 'y' from a known exceedance rate
#' function (such as a Gutenberg-Richter distribution with fixed parameters),
#' and we observe 'y + obs_error' where obs_error has a known distribution
#' dependent on the value of 'y'. This function can calculate the exceedance
#' rate of the observations. \cr
#' Written another way, suppose: \cr
#' the rate of true values that exceed y is exceedance_true(y); \cr
#' cdf_obs_error(x, y) gives the probability that obs_error > x, given that the true value is 'y'; \cr
#' the observed_value is (y+obs_err); \cr
#' Then this function computes the rate at which observed_value is above data_value, i.e. \cr
#' integral { (1 - cdf_obs_error(data_value - y | y)) * (-d(exceedance_true)/d(true_value) ) } dy \cr
#' This is done with numerical integration. \cr
#' The user should specify the limits of the integration
#' (true_value_range_with_nontrivial_cdf_value). This must be such that 
#' if true_value = true_value_range_with_nontrivial_cdf_value[1], then it 
#' is impossible for an error to make the observed value
#' greater than data_value; and if true_value = true_value_range_with_nontrivial_cdf_value[2], 
#' then it is impossible for an error to make the observed value less than data_value.
#'
#' @param exceedance_true exceedance rate function R(y) giving the rate of
#' events with (true_value > y)
#' @param cdf_obs_error function F(x,y) giving the CDF of the error (x)
#' assuming the true_value (y) is fixed
#' @param data_value the function will return the rate of observations above data_value
#' @param true_value_range_with_nontrivial_cdf_value range of the true_value over 
#' which the 'cdf complement' term \cr
#' (1 - cdf_obs_error(data_value - y | y)) \cr
#' might be in (0,1). Supposing cdf_obs_error has a finite support, then for
#' large enough y, this term will always be 1, and for small enough y it will always
#' be zero. By providing these limits we make the integration more efficient. 
#' If it does not have finite support then the integration will only be approximate
#' (but it may be very good if the probability of large errors is sufficiently small outside
#' the specified range -- this happens in the first example).
#' @param integration_dy numerically integrate with this spacing
#' @return value of the integral, which gives the rate of 'observed' (= true + error) that
#' are greater than data_thresh
#' @importFrom stats ecdf
#' @export
#' @examples
#'
#' # First example: Make a problem with known rate functions and error models,
#' # then compute the answer in multiple ways, and check it works 
#' N = 1e+06 # Also compare with large random sample
#' # Samples of the true value. This has exceedance rate function '(1-pexp)',
#' # supposing the overall rate of events is 1
#' y = rexp(N) 
#' # Samples of the error, with a distribution that depends on the true value.
#' # So this has cdf(x,y)=pgamma(x, shape=1, scale=2+y)
#' delta = rgamma(N, shape=1, scale=2+y) 
#' # The observed value
#' observed = y + delta 
#' data_thresh = 4 # We want the rate of observed>data_thresh
#' 
#' # The example allows us to compute rate(observed>data_thresh) empirically,
#' # and analytically, and check that exceedance_rate_of_observed also gives the right answer
#' 
#' # Integration of this function gives the analytical solution
#' f<-function(y,x) dexp(y) * (1-pgamma(x-y, shape=1, scale=2+y))
#' # These 2 numbers should be equal for a large enough sample
#' pr_observed_gt_4_A = mean(observed > data_thresh)
#' pr_observed_gt_4_B = integrate(f, lower=0, upper=20, x=data_thresh, rel.tol=1e-08)
#' # Further, the integral is also 'the same' (for large enough N) as
#' pr_observed_gt_4_C = mean(1-pgamma(data_thresh - y, shape=1, scale=2+y))
#' 
#' #
#' # Here we evaluate the solution with the current function
#' ex_true<-function(y) (1-pexp(y))
#' cdf_obs_error<-function(x, y) v1 = pgamma(x, shape=1, scale=2+y)
#' # Even though the error density for this problem has infinite support,
#' # the probability of 'large' errors is small enough that we can get a practically
#' # exact answer using a 'large' true_value_range_with_nontrivial_cdf_value.
#' val = exceedance_rate_of_observed(ex_true, cdf_obs_error, data_thresh,
#'    true_value_range_with_nontrivial_cdf_value = c(0, 20))
#' stopifnot(abs(val - pr_observed_gt_4_B$value) < 5.0e-06)
#' 
#' 
#' #
#' # Second example: More earthquake-like. We look at how data errors
#' # change the exceedance rate of observed events, compared to the 
#' # rate without error. For Gutenberg Richter type models, the effect
#' # of symmetric data errors is to increase the rate of 'observed' exceedances, because
#' # there are more smaller events (so it is more common for 'low value + positive error' to
#' # cause an exceedance, than for 'high value + negative error' to cause a non-exceedance)
#' #
#' 
#' # A hypothetical 'true' exceedance rate function of GR type
#' exceedance_true<-function(x){
#'     Mw_exceedance_rate_truncated_gutenberg_richter(x, a=6, b=1, Mw_min=-Inf, Mw_max=9.4)
#' }
#' # Compare results with negligably small and large errors
#' er1<-function(x, mw_true) punif(x, -0.0001, 0.0001)
#' er2<-function(x, mw_true) punif(x, -0.1, 0.1)
#' er3<-function(x, mw_true) punif(x, -0.4, 0.4)
#' 
#' # Compare the 'no error' case to the 'negligable error' case
#' data_value=7.6
#' exceedance_no_error = exceedance_true(data_value)
#' exceedance_small_error = exceedance_rate_of_observed(exceedance_true, er1, data_value=data_value,
#'    true_value_range_with_nontrivial_cdf_value=c(data_value-0.001, data_value+0.001), 
#'    integration_dy=1e-05) # Note true_value_range_with_nontrivial_cdf_value is too large (no problem, just less efficient) 
#' stopifnot(abs(exceedance_no_error - exceedance_small_error) < 1.0e-08)
#' 
#' # Next case has a small but non-negligable error
#' exceedance_mid_error = exceedance_rate_of_observed(exceedance_true, er2, data_value=data_value,
#'    true_value_range_with_nontrivial_cdf_value=c(data_value-0.15, data_value+0.15), 
#'    integration_dy=1e-03)
#' # Should be slightly larger than the 'no error' case, because smaller events
#' # are more frequent, so there are more 'smaller events with positive errors'
#' # than there are 'higher events with negative errors'
#' stopifnot(exceedance_mid_error > exceedance_no_error)
#' # Compare with another numerical integral
#' f<-function(x) -(exceedance_true(x+1.0e-06) - exceedance_true(x-1.0e-06))/(2e-06) * (1 - er2(data_value - x, x))
#' target_value = integrate(f, lower=7.4, upper=11, rel.tol=1e-08)
#' stopifnot(abs(target_value$value - exceedance_mid_error) < 3.0e-08)
#' 
#' # Next case has a larger error
#' exceedance_high_error = exceedance_rate_of_observed(exceedance_true, er3, data_value=data_value,
#'    true_value_range_with_nontrivial_cdf_value=c(data_value-0.45, data_value+0.45), 
#'    integration_dy=1e-03)
#' # Again the rate should have increased
#' stopifnot(exceedance_high_error > exceedance_mid_error)
#' # Compare with another numerical integral
#' f<-function(x) -(exceedance_true(x+1.0e-06) - exceedance_true(x-1.0e-06))/(2e-06) * (1 - er3(data_value - x, x))
#' target_value = integrate(f, lower=7., upper=11, rel.tol=1e-08)
#' stopifnot(abs(target_value$value - exceedance_high_error) < 3.0e-08)
#'  
exceedance_rate_of_observed<-function(exceedance_true, cdf_obs_error, data_value,
    true_value_range_with_nontrivial_cdf_value, integration_dy=0.01){

    stopifnot(length(data_value) == 1)
    # Numerically compute the integral (over y) of 
    #   (1 - cdf_obs_error(data_value - y, y))*(derivative_of_exceedance_true_at_y)
    # Initially try a small-ish range of mw-error to reduce the expense
    # of the integration
    lower= true_value_range_with_nontrivial_cdf_value[1]
    upper = true_value_range_with_nontrivial_cdf_value[2]
    value_true = seq(lower, upper, by=integration_dy)
    cdf_complement = 1 - cdf_obs_error(data_value - value_true, value_true)

    incremental_rate = - (
        exceedance_true(value_true+integration_dy/2) - 
        exceedance_true(value_true-integration_dy/2) )
    # Numerical integration of the 'observed rate' of events > data_value
    # Note we add in the rate of events above 'upper' (which should never have errors making them
    # below data_value)
    output = sum(incremental_rate * cdf_complement) + exceedance_true(upper + integration_dy/2)
    return(output)
}



#' Conditional empirical cumulative distribution function
#'
#' Make a function to which returns the ecdf of 'x' conditional on the values of
#' 'conditional_var'. The input 'conditional_var' must be drawn from a finite set of
#' evenly spaced values.  For instance, the 'conditional_var' might take values
#' in 7.1, 7.2, 7.3, ...  9.3. All values in this set must occur at least once
#' in 'conditional_var'.\cr
#' These are restrictive assumptions, designed a specific application. This
#' function is not really a 'general purpose' conditional ecdf creator. \cr
#' The code works by grouping values of 'x' by their value of 'conditional_var',
#' and computing the ecdf for each group. The returned function
#' F(x,conditional_var) then uses these grouped ecdf's to estimate the ecdf of 'x'
#' given *any* value of 'conditional_var', by interpolation.
#'
#' @param x numeric vector of x values for which we wish to compute the conditional ecdf
#' @param conditional_var numeric vector of conditional_var values satisfying the constraints described above
#' @return a function F(x,y) giving the ecdf of x, conditional on the fact that conditional_var=y. 
#'
#' @export
#' @examples
#' 
#' # Make conditional_vars that satisfies the assumptions of this routine
#' unique_conditional_vars = seq(7.2, 9.3, by=0.1)
#' conditional_vars = rep(unique_conditional_vars, each=50)
#' 
#' # Make 'x' a symmetric perturbation about conditional_vars
#' x = rep(seq(-0.5, 0.5, length=50), length.out=length(conditional_vars)) + conditional_vars
#' 
#' # Make the ecdf
#' ecdf_conditional = make_conditional_ecdf(x, conditional_vars)
#' 
#' #
#' # Below we perform various tests on ecdf_conditional. These give useful
#' # insight into its properties, which users should understand. 
#' #
#' 
#' # By construction in this example, there is a 50% chance of being <= the
#' # conditional value
#' stopifnot(all(ecdf_conditional(unique_conditional_vars, unique_conditional_vars) == 0.5))
#' 
#' # By construction, 40% of values should be beneath conditional_var-0.1
#' m1 = ecdf_conditional(7.6, 7.7)
#' stopifnot(m1 == 0.4)
#' # BUT, if we evaluate the function at conditional values that are not in
#' # unique_conditional_vars, then the relation will only hold approximately, due
#' # to the interpolation involved
#' random_conditional = runif(40, min(unique_conditional_vars), max(unique_conditional_vars))
#' m2 = ecdf_conditional(random_conditional - 0.1, random_conditional)
#' stopifnot(all(abs(m2 - 0.4) < 0.02))
#' #
#' # Confirm that extremes hold at interpolated values. Because of the
#' # interpolation, the following would not work using 0.5 in place of 0.6
#' #
#' m3 = ecdf_conditional(random_conditional-0.6, random_conditional)
#' stopifnot(all(m3 == 0))
#' m4 = ecdf_conditional(random_conditional+0.6, random_conditional)
#' stopifnot(all(m4 == 1))
#' 
#' # If we evaluate the function outside of the unique_conditional_vars, it will use
#' # the nearest value
#' stopifnot(ecdf_conditional(7.4, min(unique_conditional_vars)) == 
#'           ecdf_conditional(7.4, min(unique_conditional_vars) - 10))
#' stopifnot(ecdf_conditional(9.4, max(unique_conditional_vars)) == 
#'           ecdf_conditional(9.4, max(unique_conditional_vars) + 10))
#' 
#' #
#' # Example showing how event rates in a Gutenberg-Richter function are affected
#' # by measurement errors
#' #
#' Mfd_true_mw<-function(x) Mw_exceedance_rate_truncated_gutenberg_richter(x,
#'     a=5.5, b=1, Mw_min=-Inf, Mw_max=9.3)
#' # Error uniformly distributed between +-(0.2 + Mw/50)
#' observation_error_ecdf = make_conditional_ecdf(
#'     rep(seq(-1, 1, length=50), length.out=length(conditional_vars))*(0.2 + conditional_vars/50),
#'     conditional_vars)
#' true_rate_8 = Mfd_true_mw(8.0)
#' observed_rate_8 = exceedance_rate_of_observed(
#'     exceedance_true=Mfd_true_mw, 
#'     cdf_obs_error=observation_error_ecdf, 
#'     data_value = 8.0,
#'     true_value_range_with_nontrivial_cdf_value=c(7.4, 8.6))
#' # We should have more frequent 'observed > 8' than 'true > 8' for
#' # a Gutenberg Richter Mfd with symmetric errors
#' stopifnot(observed_rate_8 > true_rate_8)
#' # Check it, by comparison with an analytical version of the ecdf
#' observation_error_ecdf2<-function(x, y) pmax(0, pmin(1, (x - (-0.2-y/50))/(0.2+y/50 - (-0.2 - y/50))))
#' observed_rate_8b = exceedance_rate_of_observed(
#'     exceedance_true=Mfd_true_mw, 
#'     cdf_obs_error=observation_error_ecdf2, 
#'     data_value = 8.0,
#'     true_value_range_with_nontrivial_cdf_value=c(7.4, 8.6))
#' rate_difference_1 = (observed_rate_8 - true_rate_8)
#' rate_difference_2 = (observed_rate_8b - true_rate_8)
#' # The answers should be close
#' stopifnot(abs(rate_difference_1 - rate_difference_2) < 0.1*rate_difference_1)
#' 
#' # The empirical ecdf should be close to the exact ecdf
#' random_x = runif(1e+05, -0.5, 0.5)
#' random_conditional = runif(1e+05, 7.2, 9.3)
#' p1 = observation_error_ecdf( random_x, random_conditional)
#' p2 = observation_error_ecdf2(random_x, random_conditional)
#' stopifnot(max(abs(p1-p2)) < 0.02)
#' 
#' 
make_conditional_ecdf<-function(x, conditional_var){

    # Figure out the range of 'true' conditional_var values
    conditional_var_values_binned = sort(unique(conditional_var))

    spacing = diff(conditional_var_values_binned)

    if(!isTRUE(all.equal(spacing, rep(spacing[1], length(spacing))))){
        stop('Unique conditional_var values are not evenly spaced, which is required by make_conditional_ecdf')
    }
    spacing = spacing[1]

    # Sort the x to facilitate lookup later
    num_conditional_var = length(conditional_var_values_binned)
    x_ecdfs = vector(mode='list', length=num_conditional_var)
    for(i in 1:num_conditional_var){
        k = which(conditional_var == conditional_var_values_binned[i])
        x_ecdfs[[i]] = ecdf(x[k]) 
    }
    names(x_ecdfs) = as.character(conditional_var_values_binned)

    # Function to do the interpolation
    cdf_x_conditional<-function(xs, cond_var){
    
        # Ensure that input arguments have compatible lengths
        if(length(xs) != length(cond_var)){
            if(length(cond_var) == 1){
                cond_var = xs * 0 + cond_var
            }else if(length(xs) == 1){
                xs = cond_var*0 + xs
            }else{
                stop('Incompatible lengths of xs and cond_var')
            }
        }

        # Find lower/upper
        cond_var_lower_ind = findInterval(cond_var, conditional_var_values_binned)
        cond_var_upper_ind = cond_var_lower_ind + 1

        cond_var_lower_ind[cond_var_lower_ind < 1] = 1
        cond_var_upper_ind[cond_var_upper_ind > num_conditional_var] = num_conditional_var
    
        # Weights for interpolation 
        w_lower = (conditional_var_values_binned[cond_var_upper_ind] - cond_var)/spacing
        w_lower[w_lower < 0] = 0
        w_lower[w_lower > 1] = 1
        w_upper = (1-w_lower)

        output_lower = w_lower * 0 - 1
        output_upper = output_lower 

        # Lookup the values in each true mw
        for(i in 1:num_conditional_var){
            k = which(cond_var_lower_ind == i)
            if(length(k) > 0) output_lower[k] = x_ecdfs[[i]](xs[k])
            k = which(cond_var_upper_ind == i)
            if(length(k) > 0) output_upper[k] = x_ecdfs[[i]](xs[k])
        }

        stopifnot(all(output_lower >= 0) & all(output_upper >= 0))

        output = w_lower*output_lower + w_upper*output_upper

        return(output)
    }

    return(cdf_x_conditional)

}
